"""Check cis boundaries and shard ownership without scientific dependencies."""
import importlib.util
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/sharded_cis.py'


class ShardingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not SCRIPT.exists():
            return
        spec = importlib.util.spec_from_file_location('sharded_cis', SCRIPT)
        cls.mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = cls.mod
        spec.loader.exec_module(cls.mod)

    def setUp(self):
        self.assertTrue(SCRIPT.exists(), 'The sharded cis implementation is missing')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.pvar = self.root / 'input.pvar'
        self.pvar.write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(
            f'{c}\t{p}\t{c}_{p}\tA\tG\n' for c, p in
            [('1', 89), ('1', 90), ('1', 100), ('1', 110), ('1', 111),
             ('1', 150), ('1', 200), ('1', 300), ('2', 100)]))
        self.index = self.mod.VariantIndex(self.pvar)

    def phenotype(self, chrom, start, end, name, group=None):
        return self.mod.Phenotype(chrom, start, end, name, ('1', '2', '3'), group)

    def test_bed_start_conversion_and_inclusive_cis_boundaries(self):
        self.assertEqual(self.mod.cis_interval(99, 100, 10), (90, 110))
        self.assertEqual(self.mod.cis_interval(0, 1, 10), (1, 11))
        self.assertEqual(self.mod.cis_interval(99, 150, 10), (90, 160))
        self.assertEqual(self.index.count([('chr1', 90, 110)]), 3)

    def test_union_omits_gaps_and_counts_overlap_once(self):
        intervals = [('chr1', 90, 110), ('chr1', 100, 111), ('chr1', 290, 310)]
        self.assertEqual(self.mod.merge_intervals(intervals),
                         [('chr1', 90, 111), ('chr1', 290, 310)])
        self.assertEqual(self.index.count(intervals), 5)

    def test_variant_limit_splits_nearby_phenotypes_without_losing_windows(self):
        rows = [self.phenotype('chr1', 99, 100, 'A'),
                self.phenotype('chr1', 109, 110, 'B'),
                self.phenotype('chr1', 299, 300, 'C')]
        shards = self.mod.plan_shards(rows, self.index, 10, 10, 3)
        self.assertEqual([[p.name for p in s] for s in shards], [['A'], ['B'], ['C']])
        self.assertEqual(self.mod.shard_intervals(shards[0], 10), [('chr1', 90, 110)])

    def test_one_shard_can_contain_several_chromosomes(self):
        rows = [self.phenotype('chr2', 99, 100, 'B'),
                self.phenotype('chr1', 99, 100, 'A')]
        shards = self.mod.plan_shards(rows, self.index, 10, 10, 10)
        self.assertEqual([[p.name for p in s] for s in shards], [['A', 'B']])

    def test_group_is_indivisible_when_it_exceeds_a_target(self):
        rows = [self.phenotype('chr1', 99, 100, 'A', 'G'),
                self.phenotype('chr1', 99, 100, 'B', 'G'),
                self.phenotype('chr2', 99, 100, 'C', 'H')]
        shards = self.mod.plan_shards(rows, self.index, 10, 1, 2)
        self.assertEqual([[p.name for p in s] for s in shards], [['A', 'B'], ['C']])

    def test_group_with_different_variant_ranges_is_rejected(self):
        rows = [self.phenotype('chr1', 99, 100, 'A', 'G'),
                self.phenotype('chr1', 299, 300, 'B', 'G')]
        with self.assertRaisesRegex(ValueError, 'same cis variants'):
            self.mod.plan_shards(rows, self.index, 10, 10, 10)

    def test_ambiguous_chromosome_aliases_are_rejected(self):
        with self.pvar.open('a') as handle:
            handle.write('chr1\t400\tother\tA\tG\n')
        with self.assertRaisesRegex(ValueError, 'Ambiguous chromosome'):
            self.mod.VariantIndex(self.pvar)

    def test_invalid_coordinates_and_limits_are_rejected(self):
        for start, end, window in [(-1, 1, 1), (10, 10, 1), (10, 9, 1), (0, 1, -1)]:
            with self.subTest(start=start, end=end, window=window):
                with self.assertRaises(ValueError):
                    self.mod.cis_interval(start, end, window)
        with self.assertRaises(ValueError):
            self.mod.plan_shards([], self.index, 10, 0, 5)

    def test_localization_errors_are_reported_before_open(self):
        for path in ['gs://bucket/file', 's3://bucket/file', str(self.root/'missing')]:
            with self.assertRaisesRegex(ValueError, 'localization'):
                self.mod.local_file(path)

    def write_bed(self, rows):
        path = self.root / 'phenotypes.bed'
        path.write_text('#chr\tstart\tend\tphenotype_id\t001\t002\t003\n' + rows)
        return path

    def test_no_variant_and_constant_phenotypes_are_reported(self):
        bed = self.write_bed('1\t99\t100\tA\t1\t2\t3\n'
                             'chr3\t99\t100\tB\t1\t2\t3\n'
                             'chr1\t99\t100\tC\t1\t1\t1\n')
        samples, rows, skipped = self.mod.read_phenotypes(bed, self.index, 10)
        self.assertEqual(samples, ['001', '002', '003'])
        self.assertEqual([(p.name, p.chrom) for p in rows], [('A', 'chr1')])
        self.assertEqual([(s[0], s[2]) for s in skipped],
                         [('B', 'no_cis_variants'), ('C', 'constant_phenotype')])

    def test_bad_phenotype_values_and_duplicate_ids_fail(self):
        for text in ['1\t99\t100\tA\t1\tNA\t3\n',
                     '1\t99\t100\tA\t1\tNaN\t3\n',
                     '1\t99\t100\tA\t1\t2\t3\n1\t99\t100\tA\t2\t3\t4\n']:
            with self.subTest(text=text):
                with self.assertRaises(ValueError):
                    self.mod.read_phenotypes(self.write_bed(text), self.index, 10)

    def test_groups_must_cover_the_input_and_have_unique_members(self):
        bed = self.write_bed('1\t99\t100\tA\t1\t2\t3\n1\t99\t100\tB\t2\t3\t4\n')
        groups = self.root/'groups.tsv'
        for text in ['A\tG\n', 'A\tG\nA\tH\nB\tG\n']:
            groups.write_text(text)
            with self.assertRaises(ValueError):
                self.mod.read_phenotypes(bed, self.index, 10, groups)
        groups.write_text('B\tG\nA\tG\n')
        _, rows, _ = self.mod.read_phenotypes(bed, self.index, 10, groups)
        self.assertEqual([p.group for p in rows], ['G', 'G'])

    def test_psam_reads_iid_as_text_when_fid_is_present(self):
        psam = self.root/'input.psam'
        psam.write_text('#FID\tIID\tSEX\nF1\t001\t1\nF2\t002\t2\n')
        self.assertEqual(self.mod.sample_ids(psam), ['001', '002'])

    def test_plink_keep_file_preserves_family_and_secondary_ids(self):
        psam = self.root/'input.psam'
        psam.write_text('#FID\tIID\tSID\tSEX\nF1\t001\tS1\t1\nF2\t002\tS2\t2\n')
        keep = self.root/'keep.tsv'
        self.assertTrue(hasattr(self.mod, 'write_keep_file'), 'PLINK keep-file support is missing')
        self.mod.write_keep_file(psam, ['002'], keep)
        self.assertEqual(keep.read_text(), '#FID\tIID\tSID\nF2\t002\tS2\n')
        psam.write_text('#IID\n001\n002\n')
        self.mod.write_keep_file(psam, ['001'], keep)
        self.assertEqual(keep.read_text(), '#IID\n001\n')

    def test_covariate_samples_can_be_reordered_but_not_missing(self):
        cov = self.root/'cov.tsv'
        cov.write_text('ID\t002\t001\t003\t004\nPC1\t2\t1\t3\t4\n')
        self.mod.validate_covariates(cov, ['001', '002', '003', '004'])
        with self.assertRaisesRegex(ValueError, 'missing samples'):
            self.mod.validate_covariates(cov, ['001', '005'])

    def test_cli_exposes_typed_mapping_arguments(self):
        result = subprocess.run([sys.executable, str(SCRIPT), 'map', '--help'], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('--cis-window', result.stdout)
        self.assertIn('--permutations', result.stdout)
        self.assertNotIn('--config', result.stdout)

    @unittest.skipUnless(importlib.util.find_spec('pandas'), 'pandas is required for merge IO checks')
    def test_merge_retains_ids_that_match_pandas_na_tokens(self):
        path = self.root/'raw.tsv'
        path.write_text('phenotype_id\tvariant_id\tgroup_id\tpval_perm\tpval_beta\tbeta_shape1\tbeta_shape2\n'
                        'NA\tNULL\tNone\t0.1\t0.1\t1\t2\n'
                        'NULL\tNA\tNA\t0.2\t0.2\t1\t2\n')
        self.assertTrue(hasattr(self.mod, 'load_results'), 'Raw result loader is missing')
        result = self.mod.load_results([path])
        self.assertEqual(list(result.index), ['NA', 'NULL'])
        self.assertEqual(list(result.variant_id), ['NULL', 'NA'])
        self.assertEqual(list(result.group_id), ['None', 'NA'])

    @unittest.skipUnless(importlib.util.find_spec('pandas'), 'pandas is required for smoke comparison checks')
    def test_smoke_comparison_ignores_empty_shards_without_changing_numeric_types(self):
        import pandas as pd
        spec = importlib.util.spec_from_file_location('smoke_sharded', SCRIPT.parent.parent/'tests/smoke_sharded.py')
        smoke = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(smoke)
        empty, data = self.root/'empty.tsv', self.root/'data.tsv'
        header = 'phenotype_id\tnum_var\tpval_beta\n'
        empty.write_text(header)
        data.write_text(header + 'A\t3\t0.1\nB\t1\t0.2\n')
        observed = smoke.combine_results([empty, data])
        expected = pd.DataFrame({'num_var': [3, 1], 'pval_beta': [0.1, 0.2]},
                                index=pd.Index(['A', 'B'], name='phenotype_id'))
        pd.testing.assert_frame_equal(observed, expected)


if __name__ == '__main__':
    unittest.main()
