"""Run real PLINK2/tensorQTL/R on CPU and compare sharded and unsharded calls."""
import gzip
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from test_sharded_workflow import render_task

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = Path(os.environ.get('SHARDED_CIS_SCRIPT', str(ROOT/'scripts/sharded_cis.py')))


def run_task(name, values, directory, check=True):
    directory.mkdir()
    command = render_task(name, dict(docker_image='smoke', **values), directory)
    command = command.replace('/opt/tensorqtl/sharded_cis.py', shlex.quote(str(SCRIPT)))
    return subprocess.run(['bash', '-c', command], cwd=directory, check=check,
                          capture_output=not check, text=True)


def read_result(path):
    return pd.read_csv(path, sep='\t', index_col=0,
                       keep_default_na=False,
                       dtype={'phenotype_id': str, 'variant_id': str, 'group_id': str})


def combine_results(paths):
    frames = [read_result(path) for path in paths]
    return pd.concat(frame for frame in frames if not frame.empty).sort_index()


def variant_ids(path):
    return [row.split('\t')[2] for row in path.read_text().splitlines() if not row.startswith('#')]


def check_preparation_qc(root):
    """Check observed genotype QC without running tensorQTL permutations."""
    import pgenlib
    root.mkdir()
    n, extra = 20, 20
    samples = [f'{i:06d}' for i in range(n + extra)]
    # The BED sample order differs from the PGEN order. Extra PSAM samples must
    # not affect MAC or missingness, and nonfounders must remain in the count.
    bed_samples = samples[:n][::-1]
    records = [
        ('left_boundary', 90, [1]*4 + [0]*16, [0]*extra),
        ('low_mac', 100, [1]*3 + [0]*17, [0]*extra),
        ('high_missing', 101, [1]*4 + [0]*10 + [-9]*6, [0]*extra),
        ('high_alt', 102, [2]*19 + [1], [0]*extra),
        ('subset_mac', 103, [1] + [0]*19, [2]*extra),
        ('subset_missing', 104, [1]*5 + [0]*15, [-9]*extra),
        # Observed MAC is 3. Mean imputation would increase its allele sum to 4.
        ('before_imputation', 105, [1]*3 + [0]*12 + [-9]*5, [0]*extra),
        ('dosage_only', 106, [0.25]*n, [0]*extra),
        ('hardcall_only', 107, [0.9375]*4 + [0]*16, [0]*extra),
        ('right_boundary', 110, [1]*4 + [0]*11 + [-9]*5, [0]*extra),
        ('no_observed', 200, [-9]*n, [1]*extra),
        ('outside_windows', 1000, [1]*n, [0]*extra),
    ]
    with pgenlib.PgenWriter(os.fsencode(root/'input.pgen'), sample_ct=len(samples),
                           variant_ct=len(records), dosage_present=True) as writer:
        for _, _, analysis, other in records:
            writer.append_dosages(np.array(analysis + other, dtype=np.float32))
    pvar_header = '#CHROM\tPOS\tID\tREF\tALT\n'
    pvar_rows = [f'1\t{pos}\t{name}\tA\tG\n' for name, pos, _, _ in records]
    (root/'input.pvar').write_text(pvar_header + ''.join(pvar_rows))
    (root/'input.psam').write_text('#FID\tIID\tPAT\tMAT\tSEX\n' + ''.join(
        f'FAM\t{sample}\t' + ('0\t0' if i < 2 else f'{samples[0]}\t{samples[1]}') +
        f'\t{2 if i == 1 else 1}\n' for i, sample in enumerate(samples)))
    pd.DataFrame([
        ['chr1', 99, 100, 'QC', *range(n)],
        ['chr1', 199, 200, 'EMPTY_AFTER_QC', *range(n)],
        ['chr1', 399, 400, 'NO_SOURCE_VARIANTS', *range(n)],
    ], columns=['#chr', 'start', 'end', 'phenotype_id', *bed_samples]).to_csv(
        root/'phenotypes.bed.gz', sep='\t', index=False)
    pd.DataFrame([np.linspace(-1, 1, n)], index=['PC1'], columns=bed_samples).to_csv(
        root/'covariates.tsv', sep='\t')
    base = dict(plink_pgen=str(root/'input.pgen'), plink_pvar=str(root/'input.pvar'),
                plink_psam=str(root/'input.psam'), phenotype_bed=str(root/'phenotypes.bed.gz'),
                covariates=str(root/'covariates.tsv'), cis_window=10,
                min_mac=4, max_missing_fraction=0.25)
    for dosages in [False, True]:
        prep = root/('dosages' if dosages else 'hardcalls')
        run_task('prepare_cis', dict(base, dosages=dosages), prep)
        retained = ['left_boundary', 'subset_missing',
                    'dosage_only' if dosages else 'hardcall_only', 'right_boundary']
        assert (prep/'qc.snplist').read_text().splitlines() == retained
        assert variant_ids(prep/'qc.pvar') == retained
        summary = pd.read_csv(prep/'variant_qc.tsv', sep='\t').iloc[0]
        assert summary.source_variants == len(records)
        assert summary.candidate_variants == len(records) - 1
        assert summary.outside_cis_windows == 1
        assert summary.retained_variants == len(retained)
        assert summary.samples == n and summary.min_mac == 4
        assert summary.max_missing_fraction == 0.25
        assert summary.genotype_mode == ('dosages' if dosages else 'hard_calls')
        assert summary.low_mac == 6
        assert summary.high_missing == (2 if dosages else 3)
        assert summary.both_filters == (1 if dosages else 2)
        assert summary.no_observed_calls == (1 if dosages else 2)
        shards = json.loads((prep/'shards.json').read_text())
        assert len(shards) == 1, shards
        shard = shards[0]
        assert variant_ids(prep/shard['plink_pvar']) == retained
        assert len((prep/shard['plink_psam']).read_text().splitlines()) == n + 1
        assert list(pd.read_csv(prep/shard['phenotype_bed'], sep='\t').phenotype_id) == ['QC']
        manifest = pd.read_csv(prep/'shards.manifest.tsv', sep='\t').iloc[0]
        assert manifest.variants == len(retained) and manifest.samples == n
        skipped = pd.read_csv(prep/'preparation.skipped.tsv', sep='\t').set_index('phenotype_id')
        assert set(skipped.index) == {'EMPTY_AFTER_QC', 'NO_SOURCE_VARIANTS'}
        assert skipped.loc['EMPTY_AFTER_QC', 'reason'] == 'no_cis_variants_after_qc'
    # No valid variant may produce a successful preparation with an empty scatter.
    failed_dir = root/'all_filtered'
    failed = run_task('prepare_cis', dict(base, min_mac=n + 1), failed_dir, check=False)
    assert failed.returncode != 0
    assert 'No valid phenotypes remain after variant QC' in failed.stdout + failed.stderr
    assert not list(failed_dir.glob('shards/*.pgen'))
    assert not (failed_dir/'shards.json').exists()
    for label, marker in [('duplicate', 'duplicate'), ('multiallelic', 'multiallelic')]:
        rows = list(pvar_rows)
        if label == 'duplicate':
            rows[0] = rows[0].replace('left_boundary', 'right_boundary')
        else:
            rows[0] = rows[0].replace('\tG\n', '\tG,T\n')
        invalid = root/(label + '.pvar')
        invalid.write_text(pvar_header + ''.join(rows))
        failed = run_task('prepare_cis', dict(base, plink_pvar=str(invalid)), root/label, check=False)
        assert failed.returncode != 0 and marker in (failed.stdout + failed.stderr).lower()


def main():
    import pgenlib
    rng = np.random.default_rng(11)
    n = 96
    samples = [f'{i:06d}' for i in range(n)]
    # Phenotype at BED [99,100), window 10: include 90 and 110, exclude 89 and 111.
    positions = [('1', p) for p in [89, 90, 100, 110, 111, 200, 300, 1000]] + [('2', 100), ('2', 200)]
    genotypes = rng.binomial(2, 0.35, size=(len(positions), n)).astype(np.int8)
    # The remote variant is monomorphic; preparation must omit its empty window.
    genotypes[7, :] = 0
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        check_preparation_qc(root/'preparation_qc')
        with pgenlib.PgenWriter(os.fsencode(root/'input.pgen'), sample_ct=n, variant_ct=len(positions)) as writer:
            for row in genotypes:
                writer.append_biallelic(row)
        dosage_matrix = np.where(genotypes == 2, 1.75, genotypes + 0.25).astype(np.float32)
        dosage_matrix[7, :] = 0
        with pgenlib.PgenWriter(os.fsencode(root/'dosages.pgen'), sample_ct=n,
                               variant_ct=len(positions), dosage_present=True) as writer:
            for row in dosage_matrix:
                writer.append_dosages(row)
        (root/'input.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(
            f'{chrom}\t{pos}\tv{i}\tA\tG\n' for i, (chrom, pos) in enumerate(positions)))
        (root/'input.psam').write_text('#FID\tIID\n' + ''.join(f'F{i}\t{s}\n' for i, s in enumerate(samples)))
        ids = [('chr1', 99, 100, 'A', 2), ('chr1', 199, 200, 'B', 5),
               ('chr1', 299, 300, 'C', 6), ('chr1', 999, 1000, 'MONO', 7),
               ('chr2', 99, 100, 'D', 8), ('chr2', 199, 200, 'E', 9)]
        records = [[c, start, end, pid, *(0.5*genotypes[v] + rng.normal(size=n))]
                   for c, start, end, pid, v in ids]
        records += [['chr3', 9, 10, 'NO_VARIANTS', *rng.normal(size=n)],
                    ['chr1', 99, 100, 'CONSTANT', *np.ones(n)]]
        phenotype = pd.DataFrame(records, columns=['#chr', 'start', 'end', 'phenotype_id', *samples])
        phenotype.to_csv(root/'phenotypes.bed.gz', sep='\t', index=False)
        # Covariates have reversed sample order. The mapper must align them by text ID.
        pd.DataFrame([rng.normal(size=n)], index=['PC1'], columns=samples).loc[:, samples[::-1]].to_csv(
            root/'covariates.tsv', sep='\t')
        base = dict(plink_pgen=str(root/'input.pgen'), plink_pvar=str(root/'input.pvar'),
                    plink_psam=str(root/'input.psam'), phenotype_bed=str(root/'phenotypes.bed.gz'),
                    covariates=str(root/'covariates.tsv'), cis_window=10)
        for grouped, dosages in [(False, False), (True, False), (False, True)]:
            values = dict(base, dosages=dosages)
            if dosages:
                values['plink_pgen'] = str(root/'dosages.pgen')
            if grouped:
                # Two phenotypes in a group must remain together even with a one-row target.
                grouped_bed = pd.concat([phenotype, phenotype.iloc[[0]].assign(phenotype_id='A2')])
                grouped_bed.to_csv(root/'grouped.bed.gz', sep='\t', index=False)
                (root/'groups.tsv').write_text(''.join(f'{pid}\t{"G_A" if pid in ["A", "A2"] else "G_"+pid}\n'
                                                       for pid in grouped_bed.phenotype_id))
                values.update(phenotype_bed=str(root/'grouped.bed.gz'), phenotype_groups=str(root/'groups.tsv'))
            work = root/('dosages' if dosages else 'grouped' if grouped else 'ordinary')
            work.mkdir()
            prep = work/'prepare'
            run_task('prepare_cis', dict(values, max_phenotypes_per_shard=1, max_variants_per_shard=3), prep)
            shards = json.loads((prep/'shards.json').read_text())
            assert len(shards) == 5, shards
            first = shards[0]
            assert [r.split('\t')[2] for r in (prep/first['plink_pvar']).read_text().splitlines()[1:]] == ['v1','v2','v3']
            # Verify all extracted genotypes, not only metadata. PLINK can reorder samples.
            original_lookup = {f'v{i}': row for i, row in enumerate(dosage_matrix if dosages else genotypes)}
            for shard in shards:
                extracted_samples = (prep/shard['plink_psam']).read_text().splitlines()[1:]
                sample_order = [samples.index(s) for s in extracted_samples]
                variant_ids = [r.split('\t')[2] for r in (prep/shard['plink_pvar']).read_text().splitlines()[1:]]
                with pgenlib.PgenReader(os.fsencode(prep/shard['plink_pgen'])) as reader:
                    for i, variant in enumerate(variant_ids):
                        actual = np.empty(n, dtype=np.float32 if dosages else np.int8)
                        if dosages:
                            reader.read_dosages(i, actual)
                        else:
                            reader.read(i, actual)
                        np.testing.assert_allclose(actual, original_lookup[variant][sample_order], atol=1e-4)
            skipped = pd.read_csv(prep/'preparation.skipped.tsv', sep='\t')
            assert set(skipped.phenotype_id) == {'NO_VARIANTS', 'CONSTANT', 'MONO'}
            assert skipped.set_index('phenotype_id').loc['MONO', 'reason'] == 'no_cis_variants_after_qc'
            results = []
            for shard in shards:
                inputs = {k: str(prep/v) for k, v in shard.items() if k != 'shard_id'}
                target = work/shard['shard_id']
                run_task('map_cis', dict(inputs, covariates=base['covariates'], prefix='result',
                                        cis_window=10, permutations=256, seed=17, gpu_count=0, cpu=1,
                                        dosages=dosages), target)
                result = target/'result.cis_qtl.txt.gz'
                assert 'qval' not in read_result(result).columns
                results.append(str(result))
            # One large shard is the unsharded reference, with identical preparation filters.
            whole = work/'whole_prepare'
            run_task('prepare_cis', dict(values, max_phenotypes_per_shard=100, max_variants_per_shard=100), whole)
            ref_shards = json.loads((whole/'shards.json').read_text())
            assert len(ref_shards) == 1
            assert set(pd.read_csv(whole/ref_shards[0]['phenotype_bed'], sep='\t')['#chr']) == {'chr1', 'chr2'}
            ref = work/'reference'
            ref_inputs = {k: str(whole/v) for k, v in ref_shards[0].items() if k != 'shard_id'}
            run_task('map_cis', dict(ref_inputs, covariates=base['covariates'], prefix='reference',
                                    cis_window=10, permutations=256, seed=17, gpu_count=0, cpu=1,
                                    dosages=dosages), ref)
            expected = read_result(ref/'reference.cis_qtl.txt.gz').sort_index()
            observed = combine_results(results)
            pd.testing.assert_frame_equal(observed, expected, rtol=2e-4, atol=1e-6)
            assert set(observed.index) == {'A', 'B', 'C', 'D', 'E'}
            if not grouped and not dosages:
                # Keep the map-only empty-result check after preparation removes
                # monomorphic windows. A stricter mapping MAF removes this shard.
                empty = work/'empty_mapping'
                first_inputs = {k: str(prep/v) for k, v in shards[0].items() if k != 'shard_id'}
                run_task('map_cis', dict(first_inputs, covariates=base['covariates'], prefix='empty',
                                        cis_window=10, permutations=256, seed=17, gpu_count=0, cpu=1,
                                        maf_threshold=0.499999), empty)
                assert read_result(empty/'empty.cis_qtl.txt.gz').empty
                empty_skips = pd.read_csv(empty/'empty.skipped.tsv', sep='\t')
                assert set(empty_skips.reason) == {'no_valid_variants_after_filters'}
            merged = work/'merge'
            run_task('merge_cis', dict(results=results, prefix='merged', qvalue_lambda=0.0), merged)
            reference_merge = work/'reference_merge'
            run_task('merge_cis', dict(results=[str(ref/'reference.cis_qtl.txt.gz')],
                                      prefix='merged', qvalue_lambda=0.0), reference_merge)
            pd.testing.assert_frame_equal(read_result(merged/'merged.cis_qtl.txt.gz').sort_index(),
                                          read_result(reference_merge/'merged.cis_qtl.txt.gz').sort_index(),
                                          rtol=2e-4, atol=1e-6)
            assert read_result(merged/'merged.cis_qtl.txt.gz').qval.between(0, 1).all()
            duplicate_list = merged/'duplicate.list'
            duplicate_list.write_text(results[0] + '\n' + results[0] + '\n')
            duplicate = subprocess.run([sys.executable, str(SCRIPT), 'merge', '--results-list',
                str(duplicate_list), '--prefix', 'bad', '--qvalue-lambda', '0'],
                cwd=merged, text=True, capture_output=True)
            assert duplicate.returncode != 0 and 'Duplicate' in duplicate.stderr, duplicate.stderr
        print('Real PLINK2/tensorQTL/R smoke passed: pre-shard MAC and missingness QC, sample subsets, '
              'exact boundaries, genotype integrity, text sample IDs, complete groups, empty mapping, '
              'mixed chromosomes, fractional dosages, raw permutations, and global Storey q-values.')


if __name__ == '__main__':
    main()
