"""Exercise pre-shard QC with real PGEN calls and analysis sample subsets."""
from contextlib import chdir
import csv
import importlib.util
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/sharded_cis.py'
HAS_PGEN = all(importlib.util.find_spec(name) for name in ['numpy', 'pgenlib'])


@unittest.skipUnless(HAS_PGEN, 'numpy and pgenlib are required for real genotype QC')
class VariantQcTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('qc_sharded_cis', SCRIPT)
        cls.mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = cls.mod
        spec.loader.exec_module(cls.mod)

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.samples = ['003', '001', '000', '002']

    def fixture(self, matrix, dosages=False):
        import numpy as np
        import pgenlib
        pgen, pvar, psam = (self.root / ('input.'+ext) for ext in ['pgen', 'pvar', 'psam'])
        with pgenlib.PgenWriter(os.fsencode(pgen), sample_ct=6, variant_ct=len(matrix),
                               dosage_present=dosages) as writer:
            for row in matrix:
                if dosages:
                    writer.append_dosages(np.array(row, dtype=np.float32))
                else:
                    writer.append_biallelic(np.array(row, dtype=np.int8))
        pvar.write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(
            f'1\t{100+i}\tv{i}\tA\tG\n' for i in range(len(matrix))))
        # Pedigree metadata must not exclude nonfounders from QC counts.
        psam.write_text('#FID\tIID\tPAT\tMAT\n' + ''.join(
            f'F\t{i:03d}\tP\tM\n' for i in range(6)))
        return SimpleNamespace(plink_pgen=pgen, plink_pvar=pvar, plink_psam=psam,
                               min_mac=2, max_missing_fraction=0.25, dosages=dosages, cis_window=0)

    def test_qc_uses_only_analysis_samples_and_reports_overlapping_failures(self):
        args = self.fixture([[0, 1, 1, 0, -9, -9], [0, 0, 0, 0, 2, 2],
                             [-9, 0, 1, 1, 2, 2], [-9, -9, 1, 1, 2, 2],
                             [2, 2, 2, 1, 0, 0], [-9, -9, -9, -9, 0, 1],
                             [0, 1, 1, 0, 0, 0]])
        rows = [self.mod.Phenotype('chr1', 99, 105, 'A', ('0', '1', '2', '3'))]
        with chdir(self.root):
            index = self.mod.filter_variants(args, self.samples, self.mod.VariantIndex(args.plink_pvar), rows)
            self.assertEqual(index.total, 2)
            self.assertEqual(Path('qc.snplist').read_text(), 'v0\nv2\n')
            self.assertEqual([(r[1], r[2]) for r in self.mod.pvar_rows('qc.pvar')], [('100', 'v0'), ('102', 'v2')])
            with Path('variant_qc.tsv').open() as handle:
                summary = next(csv.DictReader(handle, delimiter='\t'))
            for key, expected in dict(source_variants=7, candidate_variants=6, outside_cis_windows=1,
                                      retained_variants=2, low_mac=3, high_missing=2, both_filters=1,
                                      no_observed_calls=1, samples=4).items():
                self.assertEqual(int(summary[key]), expected, key)
            self.assertEqual(summary['genotype_mode'], 'hard_calls')

    def test_qc_reads_dosages_when_mapping_uses_dosages(self):
        args = self.fixture([[0.25]*6], dosages=True)
        args.min_mac = 1
        rows = [self.mod.Phenotype('chr1', 99, 100, 'A', ('0', '1', '2', '3'))]
        with chdir(self.root):
            index = self.mod.filter_variants(args, self.samples, self.mod.VariantIndex(args.plink_pvar), rows)
            self.assertEqual(index.total, 1)
            args.dosages = False
            index = self.mod.filter_variants(args, self.samples, self.mod.VariantIndex(args.plink_pvar), rows)
            self.assertEqual(index.total, 0)

    def test_qc_keeps_variant_metadata_aligned_across_batches(self):
        args = self.fixture([[0, 1, 1, 0, 0, 0]] * 513)
        rows = [self.mod.Phenotype('chr1', 99, 612, 'A', ('0', '1', '2', '3'))]
        with chdir(self.root):
            index = self.mod.filter_variants(args, self.samples, self.mod.VariantIndex(args.plink_pvar), rows)
            self.assertEqual(index.total, 513)
            observed = list(self.mod.pvar_rows('qc.pvar'))
            self.assertEqual([(r[1], r[2]) for r in observed[255:258]],
                             [('355', 'v255'), ('356', 'v256'), ('357', 'v257')])
            self.assertEqual(observed[-1][2], 'v512')

    def test_all_filtered_input_reports_skips_before_creating_shards(self):
        args = self.fixture([[0, 0, 0, 0, 2, 2]])
        args.phenotype_bed = self.root/'phenotypes.bed'
        args.phenotype_bed.write_text('#chr\tstart\tend\tphenotype_id\t003\t001\t000\t002\n'
                                      '1\t99\t100\tA\t0\t1\t2\t3\n')
        args.covariates = self.root/'covariates.tsv'
        args.covariates.write_text('ID\t003\t001\t000\t002\nPC1\t0\t1\t0\t1\n')
        args.phenotype_groups = None
        args.max_phenotypes_per_shard = 1
        args.max_variants_per_shard = 10
        with chdir(self.root):
            with self.assertRaisesRegex(ValueError, 'after variant QC'):
                self.mod.prepare(args)
            self.assertFalse(Path('shards').exists())
            self.assertFalse(Path('shards.json').exists())
            self.assertIn('A\t\tno_cis_variants_after_qc', Path('preparation.skipped.tsv').read_text())
            self.assertTrue(Path('variant_qc.tsv').exists())


if __name__ == '__main__':
    unittest.main()
