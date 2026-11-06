"""Check shard output discovery and File paths across the cloud task boundary."""
import json
from pathlib import Path
import tempfile
import unittest

import WDL

DESCRIPTOR = Path(__file__).resolve().parents[1] / 'workflows/tensorqtl_cis_permutations_sharded.wdl'


class OutputStdLib(WDL.StdLib.TaskOutputs):
    def __init__(self, root):
        super().__init__('1.0', write_dir=str(root))
        self.root = root
        self.glob = WDL.StdLib.StaticFunction(
            'glob', [WDL.Type.String()], WDL.Type.Array(WDL.Type.File()), self.find_files)

    def find_files(self, pattern):
        return WDL.Value.Array(WDL.Type.File(), [
            WDL.Value.File(str(p.relative_to(self.root)))
            for p in sorted(self.root.glob(pattern.value)) if p.is_file()])

    def _devirtualize_filename(self, filename):
        return str(self.root / filename)

    def _virtualize_filename(self, filename):
        return filename


def preparation_fixture(root):
    suffixes = {'plink_pgen': '.pgen', 'plink_pvar': '.pvar', 'plink_psam': '.psam',
                'phenotype_bed': '.phenotypes.bed.gz', 'phenotype_groups': '.groups.tsv'}
    records = []
    (root / 'shards').mkdir()
    # Creation order differs from glob order. Include the shard from the Terra failure.
    for name in ['shard_000162', 'shard_000000', 'shard_000158']:
        record = {'shard_id': name}
        for field, suffix in suffixes.items():
            path = 'shards/' + name + suffix
            (root / path).write_text('' if field == 'phenotype_groups' else name)
            record[field] = path
        records.append(record)
    (root / 'shards.json').write_text(json.dumps(records))
    for name in ['shards.manifest.tsv', 'preparation.skipped.tsv', 'preparation.log']:
        (root / name).write_text('test\n')


def evaluate_preparation_outputs(doc, root):
    task = next(t for t in doc.tasks if t.name == 'prepare_cis')
    stdlib = OutputStdLib(root)
    outputs = WDL.Env.Bindings()
    for decl in task.outputs:
        outputs = outputs.bind(decl.name, decl.expr.eval(outputs, stdlib).coerce(decl.type))
    return task, outputs, stdlib


class ShardOutputCollectionTests(unittest.TestCase):
    def test_preparation_files_are_discoverable_before_the_command_runs(self):
        # Cromwell Batch predicts transfers without IO. read_json exposes only
        # its JSON argument at that point, whereas a literal glob exposes files.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            preparation_fixture(root)
            task, outputs, stdlib = evaluate_preparation_outputs(WDL.load(str(DESCRIPTOR)), root)
            transfer_files = set()

            def inspect(expr):
                if isinstance(expr, WDL.Expr.Apply) and expr.function_name in {'glob', 'read_json'}:
                    argument = expr.arguments[0].eval(WDL.Env.Bindings(), stdlib).value
                    if expr.function_name == 'glob':
                        transfer_files.update(str(p.relative_to(root)) for p in root.glob(argument) if p.is_file())
                    else:
                        transfer_files.add(argument)
                for child in expr.children:
                    inspect(child)

            for decl in task.outputs:
                inspect(decl.expr)
                if isinstance(decl.type, WDL.Type.File) and isinstance(decl.expr, WDL.Expr.String):
                    transfer_files.add(decl.expr.eval(WDL.Env.Bindings(), stdlib).value)
            returned_files = set()
            for binding in outputs:
                def record_file(value):
                    returned_files.add(value.value)
                    return value.value
                WDL.Value.rewrite_paths(binding.value, record_file)
            self.assertFalse(returned_files - transfer_files,
                             'Files absent from the cloud upload plan: ' + str(sorted(returned_files - transfer_files)))

    def test_scatter_keeps_shard_files_together_after_cloud_paths_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            preparation_fixture(root)
            doc = WDL.load(str(DESCRIPTOR))
            _, outputs, stdlib = evaluate_preparation_outputs(doc, root)
            values = WDL.values_from_json(dict(
                plink_pgen='gs://source/input.pgen', plink_pvar='gs://source/input.pvar',
                plink_psam='gs://source/input.psam', phenotype_bed='gs://source/input.bed.gz',
                covariates='gs://source/cov.tsv', prefix='all', docker_image='test:image'),
                doc.workflow.available_inputs)
            for decl in doc.workflow.inputs:
                if decl.name not in values:
                    values = values.bind(decl.name, decl.expr.eval(values, stdlib) if decl.expr else WDL.Value.Null())
            # Each glob family has its own cloud directory. Consume the returned
            # Files rather than reconstructing paths from the JSON manifest.
            for binding in outputs:
                cloud = WDL.Value.rewrite_paths(binding.value, lambda f:
                    'gs://outputs/glob-' + binding.name + '/' + Path(f.value).name)
                values = values.bind('prepare_cis.' + binding.name, cloud)
            scatter = next(n for n in doc.workflow.body if isinstance(n, WDL.Tree.Scatter))
            call = next(n for n in scatter.body if isinstance(n, WDL.Tree.Call))
            mapped = []
            for item in scatter.expr.eval(values, stdlib).value:
                env = values.bind(scatter.variable, item)
                inputs = {key: expr.eval(env, stdlib) for key, expr in call.inputs.items()}
                prefix = inputs['prefix'].value
                mapped.append(prefix)
                for key, suffix in [('plink_pgen', '.pgen'), ('plink_pvar', '.pvar'),
                                    ('plink_psam', '.psam'), ('phenotype_bed', '.phenotypes.bed.gz'),
                                    ('phenotype_groups', '.groups.tsv')]:
                    self.assertIsInstance(inputs[key], WDL.Value.File)
                    self.assertEqual(Path(inputs[key].value).name, prefix + suffix)
                    self.assertIn('gs://outputs/glob-', inputs[key].value)
            self.assertEqual(mapped, ['shard_000000', 'shard_000158', 'shard_000162'])


if __name__ == "__main__":
    unittest.main()
