"""Check Terra-safe file types and execute rendered task commands."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import WDL

ROOT = Path(__file__).resolve().parents[1]
DESCRIPTOR = ROOT/'workflows/tensorqtl_cis_permutations_sharded.wdl'


class CommandStdLib(WDL.StdLib.Base):
    def _virtualize_filename(self, filename):
        return filename

    def _devirtualize_filename(self, filename):
        return filename


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


def render_task(name, inputs, write_dir, local_paths=None):
    task = next(t for t in WDL.load(str(DESCRIPTOR)).tasks if t.name == name)
    values = WDL.values_from_json(inputs, task.available_inputs)
    if local_paths:
        values = WDL.Value.rewrite_env_paths(values, lambda value: local_paths.get(value.value, value.value))
    stdlib = CommandStdLib(task.effective_wdl_version, write_dir=str(write_dir))
    for decl in task.inputs:
        if decl.name not in {b.name for b in values}:
            values = values.bind(decl.name, decl.expr.eval(values, stdlib) if decl.expr else WDL.Value.Null())
    return task.command.eval(values, stdlib).value


class ShardedWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(DESCRIPTOR.exists(), 'The sharded WDL descriptor is missing')

    def test_wdl_types_and_no_workflow_scope_file_writes(self):
        doc = WDL.load(str(DESCRIPTOR))
        self.assertEqual(doc.wdl_version, '1.0')
        def visit(node):
            if isinstance(node, WDL.Expr.Apply):
                self.assertNotIn(node.function_name, {'write_lines', 'write_tsv', 'write_map', 'write_json'})
            for child in node.children:
                visit(child)
        visit(doc.workflow)
        prepare = next(t for t in doc.tasks if t.name == 'prepare_cis')
        outputs = {d.name: d.type for d in prepare.outputs}
        for key in ['plink_pgens', 'plink_pvars', 'plink_psams', 'phenotype_beds', 'phenotype_groups_files']:
            self.assertIsInstance(outputs[key], WDL.Type.Array)
            self.assertIsInstance(outputs[key].item_type, WDL.Type.File)
        merge = next(t for t in doc.tasks if t.name == 'merge_cis')
        types = {b.name: b.value.type for b in merge.available_inputs}
        self.assertIsInstance(types['results'].item_type, WDL.Type.File)

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

    def run_command(self, name, cloud_key=None, localize=True):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_dir = root/'command-files'
            write_dir.mkdir()
            mapping = {}
            values = dict(docker_image='test:image')
            files = (['plink_pgen', 'plink_pvar', 'plink_psam', 'phenotype_bed', 'covariates', 'phenotype_groups']
                     if name != 'merge_cis' else ['results'])
            for key in files:
                folder = root/(key + " space '$literal $(touch injected)")
                folder.mkdir()
                path = folder/'input'
                path.write_text('data\n')
                uri = 'gs://bucket/' + key
                mapping[uri] = str(path)
                value = uri if localize or cloud_key == key else str(path)
                values[key] = [value] if key == 'results' else value
            if name != 'prepare_cis':
                values['prefix'] = 'test'
            command = render_task(name, values, write_dir, mapping if localize else None)
            # The image includes this script. Substitute only its installed location.
            probe = root/'probe.py'
            probe.write_text('import json,pathlib,sys\npathlib.Path("argv.json").write_text(json.dumps(sys.argv[1:]))\n')
            command = command.replace('/opt/tensorqtl/sharded_cis.py', str(probe))
            process = subprocess.run(['bash', '-c', command], cwd=root, text=True, capture_output=True)
            args = json.loads((root/'argv.json').read_text()) if (root/'argv.json').exists() else None
            self.assertFalse((root/'injected').exists())
            list_text = None
            if args and '--results-list' in args:
                list_text = Path(args[args.index('--results-list')+1]).read_text()
            return process, args, mapping, list_text

    def test_cloud_files_localize_in_every_task(self):
        for task in ['prepare_cis', 'map_cis', 'merge_cis']:
            with self.subTest(task=task):
                process, args, mapping, list_text = self.run_command(task)
                self.assertEqual(process.returncode, 0, process.stderr)
                self.assertFalse(any('gs://' in arg for arg in args))
                if task == 'merge_cis':
                    self.assertEqual(list_text, mapping['gs://bucket/results'] + '\n')
                else:
                    for uri, path in mapping.items():
                        self.assertIn(path, args)

    def test_unresolved_uri_fails_before_script_in_every_task(self):
        for task in ['prepare_cis', 'map_cis', 'merge_cis']:
            keys = ['results'] if task == 'merge_cis' else ['plink_pgen', 'plink_pvar', 'plink_psam',
                                                                       'phenotype_bed', 'covariates', 'phenotype_groups']
            for key in keys:
                with self.subTest(task=task, key=key):
                    process, args, _, _ = self.run_command(task, cloud_key=key, localize=False)
                    self.assertNotEqual(process.returncode, 0)
                    self.assertIn('localization', (process.stdout + process.stderr).lower())
                    self.assertIsNone(args)


if __name__ == '__main__':
    unittest.main()
