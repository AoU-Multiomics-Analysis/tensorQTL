#!/usr/bin/env python3
"""Prepare complete cis windows, map permutations, and calculate global q-values."""
import argparse
from array import array
from bisect import bisect_left, bisect_right
from collections import OrderedDict
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import gzip
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys


def log(message):
    print(f'[{datetime.now(timezone.utc).isoformat()}] {message}', flush=True)


def local_file(value):
    value = os.fspath(value)
    if '://' in value:
        raise ValueError(f'Input localization error: unresolved URI {value}')
    path = Path(value)
    if not path.is_file() or not os.access(path, os.R_OK):
        raise ValueError(f'Input localization error: unreadable file {value}')
    return path.resolve()


def open_text(path):
    path = local_file(path)
    return gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()


def chromosome(value):
    """Normalize standard PLINK chromosome aliases; keep other contigs intact."""
    value = str(value)
    tail = value[3:] if value.lower().startswith('chr') else value
    upper = tail.upper()
    aliases = {'23': 'X', '24': 'Y', '25': 'XY', '26': 'M', 'MT': 'M'}
    upper = aliases.get(upper, upper)
    if upper in {'X', 'Y', 'XY', 'M'}:
        return 'chr' + upper
    if tail.isdigit() and 1 <= int(tail) <= 22:
        return 'chr' + str(int(tail))
    return value


def chrom_key(value):
    tail = value[3:] if value.startswith('chr') else value
    order = {'X': 23, 'Y': 24, 'XY': 25, 'M': 26}
    if tail.isdigit():
        return (0, int(tail))
    return (0, order[tail]) if tail in order else (1, value)


def cis_interval(start, end, window):
    if start < 0 or end <= start or window < 0:
        raise ValueError('BED coordinates require 0 <= start < end and cis_window >= 0')
    return max(1, start + 1 - window), end + window


def merge_intervals(intervals):
    result = []
    for chrom, start, end in sorted(intervals, key=lambda r: (chrom_key(r[0]), r[1], r[2])):
        if result and result[-1][0] == chrom and start <= result[-1][2] + 1:
            result[-1] = (chrom, result[-1][1], max(end, result[-1][2]))
        else:
            result.append((chrom, start, end))
    return result


def pvar_rows(path):
    with open_text(path) as handle:
        header = None
        for line in handle:
            if line.startswith('##') or not line.strip():
                continue
            fields = line.rstrip('\n').split()
            if header is None:
                if fields[:5] != ['#CHROM', 'POS', 'ID', 'REF', 'ALT']:
                    raise ValueError('PVAR must start with #CHROM POS ID REF ALT')
                header = fields
                continue
            if len(fields) < 5 or int(fields[1]) < 1:
                raise ValueError('Invalid PVAR row')
            yield fields
        if header is None:
            raise ValueError('PVAR header is missing')


class VariantIndex:
    def __init__(self, path):
        self.positions = {}
        self.labels = {}
        self.total = 0
        for row in pvar_rows(path):
            label, position = row[0], int(row[1])
            chrom = chromosome(label)
            if chrom in self.labels and self.labels[chrom] != label:
                raise ValueError(f'Ambiguous chromosome aliases: {self.labels[chrom]} and {label}')
            self.labels[chrom] = label
            positions = self.positions.setdefault(chrom, array('q'))
            if positions and position < positions[-1]:
                raise ValueError(f'PVAR positions must be sorted within {label}')
            positions.append(position)
            self.total += 1
        if not self.total:
            raise ValueError('PVAR contains no variants')

    def bounds(self, chrom, start, end):
        positions = self.positions.get(chrom, ())
        return bisect_left(positions, start), bisect_right(positions, end)

    def count(self, intervals):
        return sum(upper - lower for chrom, start, end in merge_intervals(intervals)
                   for lower, upper in [self.bounds(chrom, start, end)])


@dataclass(frozen=True)
class Phenotype:
    chrom: str
    start: int
    end: int
    name: str
    values: tuple
    group: str | None = None


def shard_intervals(rows, window):
    return merge_intervals((row.chrom, *cis_interval(row.start, row.end, window)) for row in rows)


def plan_shards(rows, index, window, max_phenotypes, max_variants):
    if window < 0 or max_phenotypes < 1 or max_variants < 1:
        raise ValueError('Shard limits must be positive; cis_window must be nonnegative')
    units = OrderedDict()
    for row in sorted(rows, key=lambda p: (chrom_key(p.chrom), p.start, p.end, p.name)):
        # Tag the key so a phenotype ID cannot collide with a group ID.
        key = ('group', row.group) if row.group is not None else ('phenotype', row.name)
        units.setdefault(key, []).append(row)
    shards, current = [], []
    for unit in units.values():
        ranges = {(p.chrom, *index.bounds(p.chrom, *cis_interval(p.start, p.end, window))) for p in unit}
        if len(ranges) > 1:
            raise ValueError(f'Group {unit[0].group} must have the same cis variants for all members')
        candidate = current + unit
        if current and (len(candidate) > max_phenotypes or
                        index.count(shard_intervals(candidate, window)) > max_variants):
            shards.append(current)
            current = []
        current.extend(unit)
    if current:
        shards.append(current)
    return [sorted(shard, key=lambda p: (chrom_key(p.chrom), p.start, p.end, p.name)) for shard in shards]


def unique_ids(values, label):
    if not values or any(not v or any(c.isspace() for c in v) for v in values):
        raise ValueError(f'{label} must be nonempty IDs without whitespace')
    if len(set(values)) != len(values):
        raise ValueError(f'Duplicate {label}')


def read_psam(path):
    with open_text(path) as handle:
        reader = (line.split() for line in handle if line.strip())
        header = next(reader, [])
        header = [value.lstrip('#') for value in header]
        if 'IID' not in header:
            raise ValueError('PSAM requires an IID header')
        column = header.index('IID')
        result, records = [], []
        for row in reader:
            if len(row) != len(header):
                raise ValueError('Invalid PSAM row')
            result.append(row[column])
            records.append(row)
    unique_ids(result, 'PSAM sample IDs')
    return header, records


def sample_ids(path):
    header, records = read_psam(path)
    return [row[header.index('IID')] for row in records]


def write_keep_file(psam, samples, output):
    header, records = read_psam(psam)
    columns = [key for key in ['FID', 'IID', 'SID'] if key in header]
    selected = set(samples)
    write_tsv(output, ['#' + columns[0], *columns[1:]],
              ([row[header.index(key)] for key in columns] for row in records
               if row[header.index('IID')] in selected))


def validate_covariates(path, samples):
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader, [])
        unique_ids(header[1:], 'covariate sample IDs')
        missing = set(samples) - set(header[1:])
        if missing:
            raise ValueError(f'Covariates have missing samples: {sorted(missing)[:5]}')
        names = []
        columns = [header.index(sample) for sample in samples]
        for row in reader:
            if len(row) != len(header):
                raise ValueError('Invalid covariate row length')
            names.append(row[0])
            if not all(math.isfinite(float(row[i])) for i in columns):
                raise ValueError('Covariates contain nonfinite values')
        unique_ids(names, 'covariate IDs')
        if len(samples) <= len(names) + 2:
            raise ValueError('Too few samples for the number of covariates')


def read_phenotypes(path, index, window, groups_path=None):
    groups = {}
    if groups_path and local_file(groups_path).stat().st_size:
        with open_text(groups_path) as handle:
            for line in handle:
                fields = line.rstrip('\n').split('\t')
                if len(fields) != 2 or not all(fields) or fields[0] in groups:
                    raise ValueError('Groups require one unique phenotype ID and group ID per row')
                groups[fields[0]] = fields[1]
    rows, skipped, names = [], [], []
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader, [])
        if [v.lower() for v in header[:4]] != ['#chr', 'start', 'end', 'phenotype_id']:
            raise ValueError('BED requires #chr, start, end, phenotype_id, and sample columns')
        samples = header[4:]
        unique_ids(samples, 'phenotype sample IDs')
        for fields in reader:
            if len(fields) != len(header):
                raise ValueError('Invalid phenotype row length')
            chrom, start, end, name = chromosome(fields[0]), int(fields[1]), int(fields[2]), fields[3]
            lower, upper = cis_interval(start, end, window)
            names.append(name)
            group = groups.get(name)
            if groups and group is None:
                raise ValueError(f'Groups file is missing phenotype {name}')
            values = tuple(float(v) for v in fields[4:])
            if not all(math.isfinite(v) for v in values):
                raise ValueError(f'Nonfinite phenotype values: {name}')
            reason = None
            if index.count([(chrom, lower, upper)]) == 0:
                reason = 'no_cis_variants'
            elif min(values) == max(values):
                reason = 'constant_phenotype'
            if reason:
                skipped.append((name, group or '', reason))
            else:
                rows.append(Phenotype(chrom, start, end, name, tuple(fields[4:]), group))
    unique_ids(names, 'phenotype IDs')
    if groups and set(groups) != set(names):
        raise ValueError('Groups must cover exactly the input phenotype IDs')
    return samples, rows, skipped


def write_tsv(path, header, rows):
    with Path(path).open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)


def write_skips(path, rows):
    write_tsv(path, ['phenotype_id', 'group_id', 'reason'], rows)


def prepare(args):
    log('stage=prepare Validate localized source files and samples.')
    for value in [args.plink_pgen, args.plink_pvar, args.plink_psam, args.phenotype_bed, args.covariates]:
        local_file(value)
    index = VariantIndex(args.plink_pvar)
    samples, rows, skipped = read_phenotypes(args.phenotype_bed, index, args.cis_window, args.phenotype_groups)
    missing = set(samples) - set(sample_ids(args.plink_psam))
    if missing:
        raise ValueError(f'Genotypes have missing samples: {sorted(missing)[:5]}')
    validate_covariates(args.covariates, samples)
    write_skips('preparation.skipped.tsv', skipped)
    shards = plan_shards(rows, index, args.cis_window, args.max_phenotypes_per_shard, args.max_variants_per_shard)
    if not shards:
        raise ValueError('No valid phenotypes remain after preparation; see preparation.skipped.tsv')
    Path('shards').mkdir()
    write_keep_file(args.plink_psam, samples, 'samples.keep.tsv')
    outputs, audit = [], []
    source_bytes = sum(local_file(v).stat().st_size for v in [args.plink_pgen, args.plink_pvar, args.plink_psam])
    for i, shard in enumerate(shards):
        name = f'shard_{i:06d}'
        stem = Path('shards') / name
        intervals = shard_intervals(shard, args.cis_window)
        ranges = stem.with_suffix('.ranges.tsv')
        # PLINK bed1 uses one-based inclusive endpoints, like tensorQTL's ranges.
        ranges.write_text(''.join(f'{c}\t{s}\t{e}\tregion_{j}\n' for j, (c, s, e) in enumerate(intervals)))
        log(f'stage=prepare Extract {name}: phenotypes={len(shard)} variants={index.count(intervals)}.')
        subprocess.run(['plink2', '--pgen', str(local_file(args.plink_pgen)),
                        '--pvar', str(local_file(args.plink_pvar)), '--psam', str(local_file(args.plink_psam)),
                        '--extract', 'bed1', str(ranges), '--keep', 'samples.keep.tsv',
                        '--allow-extra-chr', '--output-chr', 'chrM', '--make-pgen',
                        '--threads', str(args.threads), '--memory', str(args.plink_memory_mb),
                        '--out', str(stem)], check=True)
        # tensorQTL expects one text IID column; retain the actual PLINK sample order.
        retained = sample_ids(stem.with_suffix('.psam'))
        if set(retained) != set(samples):
            raise ValueError('PLINK did not retain the requested phenotype samples')
        write_tsv(stem.with_suffix('.psam'), ['#IID'], [(sample,) for sample in retained])
        selected = list(pvar_rows(stem.with_suffix('.pvar')))
        unique_ids([v[2] for v in selected], 'selected variant IDs')
        if any(',' in v[4] for v in selected):
            raise ValueError('Selected multiallelic variants must be split before tensorQTL mapping')
        selected_index = VariantIndex(stem.with_suffix('.pvar'))
        if selected_index.total != index.count(intervals):
            raise ValueError(f'PLINK extraction lost or added cis variants in {name}')
        for row in shard:
            interval = [(row.chrom, *cis_interval(row.start, row.end, args.cis_window))]
            if selected_index.count(interval) != index.count(interval):
                raise ValueError(f'Incomplete cis window for phenotype {row.name}')
        # Nonstandard contigs retain their exact input names; standard aliases are normalized.
        with gzip.open(str(stem) + '.phenotypes.bed.gz', 'wt') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['#chr', 'start', 'end', 'phenotype_id', *samples])
            writer.writerows((p.chrom, p.start, p.end, p.name, *p.values) for p in shard)
        groups_file = Path(str(stem) + '.groups.tsv')
        groups_file.write_text(''.join(f'{p.name}\t{p.group}\n' for p in shard if p.group is not None))
        record = dict(shard_id=name, plink_pgen=str(stem)+'.pgen', plink_pvar=str(stem)+'.pvar',
                      plink_psam=str(stem)+'.psam', phenotype_bed=str(stem)+'.phenotypes.bed.gz',
                      phenotype_groups=str(groups_file))
        outputs.append(record)
        genotype_bytes = sum(Path(record[k]).stat().st_size for k in ['plink_pgen', 'plink_pvar', 'plink_psam'])
        oversized = len(shard) > args.max_phenotypes_per_shard or selected_index.total > args.max_variants_per_shard
        audit.append((name, ','.join(dict.fromkeys(p.chrom for p in shard)), len(shard),
                      len({p.group if p.group is not None else p.name for p in shard}),
                      selected_index.total, len(intervals), len(samples), args.cis_window,
                      genotype_bytes, source_bytes, oversized))
        if oversized:
            log(f'stage=prepare {name} exceeds a size target to retain a complete phenotype or group.')
    # This is a task output manifest. Typed WDL struct fields collect each File.
    Path('shards.json').write_text(json.dumps(outputs, indent=2) + '\n')
    write_tsv('shards.manifest.tsv', ['shard_id', 'chromosomes', 'phenotypes', 'tests', 'variants',
              'intervals', 'samples', 'cis_window', 'genotype_bytes', 'source_genotype_bytes', 'oversized'], audit)
    log(f'stage=prepare Completed {len(shards)} shards; skipped {len(skipped)} phenotypes.')


def empty_results(grouped=False):
    import pandas as pd
    from tensorqtl.core import output_dtype_dict
    columns = dict(output_dtype_dict)
    if grouped:
        columns.update(group_id=str, group_size='int32')
    result = pd.DataFrame({key: pd.Series(dtype=dtype) for key, dtype in columns.items()})
    result.index.name = 'phenotype_id'
    return result


def map_permutations(args):
    log('stage=map Validate localized inputs and align sample IDs.')
    for value in [args.plink_pgen, args.plink_pvar, args.plink_psam, args.phenotype_bed, args.covariates]:
        local_file(value)
    index = VariantIndex(args.plink_pvar)
    samples, rows, skipped = read_phenotypes(args.phenotype_bed, index, args.cis_window, args.phenotype_groups)
    validate_covariates(args.covariates, samples)
    genotype_samples = sample_ids(args.plink_psam)
    if set(samples) - set(genotype_samples):
        raise ValueError('Genotypes have missing phenotype samples')
    # Also enforce group compatibility if this task is called outside the workflow.
    plan_shards(rows, index, args.cis_window, max(1, len(rows)), max(1, index.total))
    import numpy as np
    import pandas as pd
    import pgenlib
    import torch
    from tensorqtl import cis
    torch.set_num_threads(args.threads)
    if args.require_gpu and not torch.cuda.is_available():
        raise ValueError('A GPU was requested, but PyTorch cannot access CUDA')
    log(f'stage=map device={"cuda" if torch.cuda.is_available() else "cpu"}; seed={args.seed}; permutations={args.permutations}.')
    if not rows:
        result = empty_results(bool(args.phenotype_groups and local_file(args.phenotype_groups).stat().st_size))
    else:
        variants = list(pvar_rows(args.plink_pvar))
        unique_ids([v[2] for v in variants], 'selected variant IDs')
        # Read text sample IDs ourselves. tensorQTL 1.0.10's PSAM reader infers
        # numeric IDs before converting to str and can remove leading zeros.
        subset = np.array(sorted(genotype_samples.index(s) for s in samples), dtype=np.uint32)
        retained = [genotype_samples[i] for i in subset]
        with pgenlib.PgenReader(os.fsencode(local_file(args.plink_pgen)), sample_subset=subset) as reader:
            if reader.get_variant_ct() != len(variants) or reader.get_raw_sample_ct() != len(genotype_samples):
                raise ValueError('PGEN, PVAR, and PSAM dimensions do not match')
            matrix = np.empty((len(variants), len(subset)), dtype=np.float32 if args.dosages else np.int8)
            if args.dosages:
                reader.read_dosages_range(0, len(variants), matrix)
            else:
                reader.read_range(0, len(variants), matrix)
        genotype_df = pd.DataFrame(matrix, index=[v[2] for v in variants], columns=retained)
        variant_df = pd.DataFrame({'chrom': [chromosome(v[0]) for v in variants],
                                   'pos': [int(v[1]) for v in variants]}, index=genotype_df.index)
        phenotype_df = pd.DataFrame(np.array([p.values for p in rows], dtype=np.float32),
                                    index=[p.name for p in rows], columns=samples)
        phenotype_pos_df = pd.DataFrame({'chr': [p.chrom for p in rows],
                                        'start': [p.start + 1 for p in rows],
                                        'end': [p.end for p in rows]}, index=phenotype_df.index)
        covariates_df = pd.read_csv(local_file(args.covariates), sep='\t', index_col=0,
                                    dtype=str, keep_default_na=False).loc[:, samples].T.astype(np.float32)
        group_s = pd.Series([p.group for p in rows], index=phenotype_df.index) if rows[0].group is not None else None
        try:
            result = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                                 covariates_df=covariates_df, group_s=group_s, nperm=args.permutations,
                                 window=args.cis_window, maf_threshold=args.maf_threshold,
                                 seed=args.seed, verbose=True)
        except ValueError as error:
            # tensorQTL 1.0.10 calls pd.concat([]) when every window contains
            # only monomorphic variants or variants excluded by MAF. Handle
            # this exact state; other numeric or model errors must fail.
            trace = error.__traceback__
            empty = False
            while trace is not None:
                if (trace.tb_frame.f_code is cis.map_cis.__code__ and
                    str(error) == 'No objects to concatenate' and
                    isinstance(trace.tb_frame.f_locals.get('res_df'), list) and
                    len(trace.tb_frame.f_locals['res_df']) == 0):
                    empty = True
                trace = trace.tb_next
            if not empty:
                raise
            result = empty_results(group_s is not None)
        if group_s is None:
            tested = set(result.index)
            skipped.extend((p.name, '', 'no_valid_variants_after_filters') for p in rows if p.name not in tested)
        else:
            tested = set(result['group_id'])
            skipped.extend((p.name, p.group, 'no_valid_variants_after_filters') for p in rows if p.group not in tested)
    # Keep precision until the global q-value step. No per-shard FDR is run.
    result.to_csv(args.prefix + '.cis_qtl.txt.gz', sep='\t', float_format='%.17g')
    write_skips(args.prefix + '.skipped.tsv', skipped)
    log(f'stage=map Completed {len(result)} tests; skipped {len(skipped)} phenotypes; no shard q-values calculated.')


def load_results(paths):
    import pandas as pd
    frames = []
    schema = None
    for path in paths:
        frame = pd.read_csv(local_file(path), sep='\t', index_col=0, float_precision='round_trip',
                            keep_default_na=False,
                            dtype={'phenotype_id': str, 'variant_id': str, 'group_id': str})
        if 'qval' in frame.columns or not {'pval_perm', 'pval_beta', 'beta_shape1', 'beta_shape2'} <= set(frame.columns):
            raise ValueError('Merge requires raw cis permutation summaries without shard q-values')
        if schema is not None and list(frame.columns) != schema:
            raise ValueError('Shard result schemas do not match')
        schema = list(frame.columns)
        if not frame.empty:
            frames.append(frame)
    if not frames:
        raise ValueError('No valid tests remain after genotype filters')
    result = pd.concat(frames)
    if result.index.has_duplicates or ('group_id' in result and result['group_id'].duplicated().any()):
        raise ValueError('Duplicate phenotype or group tests across shards')
    return result


def merge_results(args):
    log('stage=merge Read all localized permutation summaries.')
    paths = local_file(args.results_list).read_text().splitlines()
    if not paths:
        raise ValueError('No shard results were supplied')
    result = load_results(paths)
    import numpy as np
    import tensorqtl
    from tensorqtl import core
    if len(result) < 2:
        raise ValueError('Global tensorQTL q-value calculation requires at least two valid tests')
    pvalues = result['pval_beta'].to_numpy(dtype=float)
    if not np.isfinite(pvalues).all() or ((pvalues < 0) | (pvalues > 1)).any():
        raise ValueError('Invalid beta-approximated p-values; inspect the mapping logs')
    if not core.has_rpy2:
        raise ValueError('Global q-values require R, the R qvalue package, and rpy2')
    log(f'stage=merge Calculate global Storey q-values for {len(result)} tests with R qvalue.')
    try:
        tensorqtl.calculate_qvalues(result, fdr=args.fdr, qvalue_lambda=args.qvalue_lambda)
    except Exception as error:
        raise ValueError('R qvalue could not estimate global FDR. Inspect the p-values and qvalue_lambda; '
                         'qvalue_lambda=0 can be used for an explicit conservative setting.') from error
    if not np.isfinite(result['qval']).all() or not result['qval'].between(0, 1).all():
        raise ValueError('R qvalue returned invalid q-values')
    result.to_csv(args.prefix + '.cis_qtl.txt.gz', sep='\t', float_format='%.17g')
    log('stage=merge Completed global q-value calculation.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    prep = commands.add_parser('prepare')
    mapping = commands.add_parser('map')
    merge = commands.add_parser('merge')
    for subparser in [prep, mapping]:
        for name in ['plink-pgen', 'plink-pvar', 'plink-psam', 'phenotype-bed', 'covariates']:
            subparser.add_argument('--' + name, required=True)
        subparser.add_argument('--phenotype-groups')
        subparser.add_argument('--cis-window', type=int, default=1000000)
        subparser.add_argument('--threads', type=int, default=4)
    prep.add_argument('--max-phenotypes-per-shard', type=int, default=250)
    prep.add_argument('--max-variants-per-shard', type=int, default=100000)
    prep.add_argument('--plink-memory-mb', type=int, default=16000)
    mapping.add_argument('--prefix', required=True)
    mapping.add_argument('--permutations', type=int, default=10000)
    mapping.add_argument('--maf-threshold', type=float, default=0)
    mapping.add_argument('--seed', type=int, default=1)
    mapping.add_argument('--dosages', action='store_true')
    mapping.add_argument('--require-gpu', action='store_true')
    merge.add_argument('--results-list', required=True)
    merge.add_argument('--prefix', required=True)
    merge.add_argument('--fdr', type=float, default=0.05)
    merge.add_argument('--qvalue-lambda', type=float)
    args = parser.parse_args()
    if hasattr(args, 'prefix') and not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', args.prefix):
        parser.error('prefix must be a filename prefix using letters, digits, dots, underscores, or hyphens')
    if hasattr(args, 'threads') and (args.threads < 1 or args.cis_window < 0):
        parser.error('threads must be positive and cis-window must be nonnegative')
    if args.command == 'map' and (args.permutations < 10 or not 0 <= args.maf_threshold < 0.5 or not 0 <= args.seed < 2**32):
        parser.error('permutations must be >= 10, 0 <= maf-threshold < 0.5, and 0 <= seed < 2^32')
    if args.command == 'merge' and (not 0 < args.fdr < 1 or
                                  args.qvalue_lambda is not None and not 0 <= args.qvalue_lambda < 1):
        parser.error('Require 0 < fdr < 1 and 0 <= qvalue-lambda < 1')
    try:
        {'prepare': prepare, 'map': map_permutations, 'merge': merge_results}[args.command](args)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        log(f'ERROR: {error}')
        print(str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
