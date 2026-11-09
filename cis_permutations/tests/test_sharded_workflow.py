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


class CloudCommandStdLib(CommandStdLib):
    """Model the cloud names Cromwell assigns to command-created files."""
    def _virtualize_filename(self, filename):
        return 'gs://cromwell-outputs/call-merge_cis/' + Path(filename).name


def render_task(name, inputs, write_dir, local_paths=None, cloud_command_files=False):
    task = next(t for t in WDL.load(str(DESCRIPTOR)).tasks if t.name == name)
    values = WDL.values_from_json(inputs, task.available_inputs)
    if local_paths:
        values = WDL.Value.rewrite_env_paths(values, lambda value: local_paths.get(value.value, value.value))
    stdlib_class = CloudCommandStdLib if cloud_command_files else CommandStdLib
    stdlib = stdlib_class(task.effective_wdl_version, write_dir=str(write_dir))
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

    def run_command(self, name, cloud_key=None, localize=True, cloud_command_files=False, result_count=1,
                    input_overrides=None):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_dir = root/'command-files'
            write_dir.mkdir()
            mapping = {}
            values = dict(docker_image='test:image')
            files = (['plink_pgen', 'plink_pvar', 'plink_psam', 'phenotype_bed', 'covariates', 'phenotype_groups']
                     if name != 'merge_cis' else ['results'])
            for key in files:
                entries = []
                for index in range(result_count if key == 'results' else 1):
                    label = key + (f'_{index}' if index else '')
                    folder = root/(label + " space '\"$literal $(touch injected) `touch injected_backtick`")
                    folder.mkdir()
                    path = folder/'input'
                    path.write_text('data\n')
                    uri = 'gs://bucket/' + label
                    mapping[uri] = str(path)
                    entries.append(uri if localize or cloud_key == key else str(path))
                values[key] = entries if key == 'results' else entries[0]
            if name != 'prepare_cis':
                values['prefix'] = 'test'
            values.update(input_overrides or {})
            command = render_task(name, values, write_dir, mapping if localize else None,
                                  cloud_command_files=cloud_command_files)
            # The image includes this script. Substitute only its installed location.
            probe = root/'probe.py'
            probe.write_text('import json,pathlib,sys\npathlib.Path("argv.json").write_text(json.dumps(sys.argv[1:]))\n')
            command = command.replace('/opt/tensorqtl/sharded_cis.py', str(probe))
            process = subprocess.run(['bash', '-c', command], cwd=root, text=True, capture_output=True)
            args = json.loads((root/'argv.json').read_text()) if (root/'argv.json').exists() else None
            self.assertFalse((root/'injected').exists())
            self.assertFalse((root/'injected_backtick').exists())
            list_text = None
            if args and '--results-list' in args:
                list_text = (root / args[args.index('--results-list')+1]).read_text()
            return process, args, mapping, list_text

    def test_preparation_passes_default_variant_qc_limits(self):
        process, args, _, _ = self.run_command('prepare_cis')
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertIn('--min-mac', args)
        self.assertEqual(args[args.index('--min-mac') + 1], '10')
        self.assertIn('--max-missing-fraction', args)
        self.assertEqual(float(args[args.index('--max-missing-fraction') + 1]), 0.05)
        self.assertNotIn('--dosages', args)

    def test_preparation_passes_qc_overrides_and_dosage_mode_after_localization(self):
        task = next(t for t in WDL.load(str(DESCRIPTOR)).tasks if t.name == 'prepare_cis')
        for name in ['min_mac', 'max_missing_fraction', 'dosages']:
            self.assertIn(name, {binding.name for binding in task.available_inputs})
        for dosages in [False, True]:
            with self.subTest(dosages=dosages):
                process, args, mapping, _ = self.run_command('prepare_cis', input_overrides={
                    'min_mac': 17, 'max_missing_fraction': 0.02, 'dosages': dosages})
                self.assertEqual(process.returncode, 0, process.stderr)
                self.assertEqual(args[args.index('--min-mac') + 1], '17')
                self.assertEqual(float(args[args.index('--max-missing-fraction') + 1]), 0.02)
                self.assertEqual('--dosages' in args, dosages)
                for path in mapping.values():
                    self.assertIn(path, args)
                self.assertFalse(any('gs://' in arg for arg in args))

    def test_workflow_passes_qc_values_to_preparation_and_exposes_report(self):
        doc = WDL.load(str(DESCRIPTOR))
        prepare = next(t for t in doc.tasks if t.name == 'prepare_cis')
        task_outputs = {decl.name: decl for decl in prepare.outputs}
        workflow_outputs = {decl.name: decl for decl in doc.workflow.outputs}
        self.assertIn('variant_qc', task_outputs)
        self.assertIn('variant_qc', workflow_outputs)
        self.assertIsInstance(task_outputs['variant_qc'].type, WDL.Type.File)
        self.assertIsInstance(workflow_outputs['variant_qc'].type, WDL.Type.File)
        call = next(node for node in doc.workflow.body if isinstance(node, WDL.Tree.Call)
                    and node.name == 'prepare_cis')
        values = WDL.Env.Bindings().bind('min_mac', WDL.Value.Int(17)).bind(
            'max_missing_fraction', WDL.Value.Float(0.02)).bind('dosages', WDL.Value.Boolean(True))
        stdlib = CommandStdLib('1.0')
        for name, expected in [('min_mac', 17), ('max_missing_fraction', 0.02), ('dosages', True)]:
            self.assertIn(name, call.inputs)
            self.assertEqual(call.inputs[name].eval(values, stdlib).value, expected)
        report_values = WDL.Env.Bindings().bind('prepare_cis.variant_qc', WDL.Value.File('/localized/variant_qc.tsv'))
        self.assertEqual(workflow_outputs['variant_qc'].expr.eval(report_values, stdlib).value,
                         '/localized/variant_qc.tsv')

    def test_merge_list_stays_local_when_generated_files_have_cloud_names(self):
        for count in [1, 3]:
            with self.subTest(count=count):
                process, args, mapping, list_text = self.run_command(
                    'merge_cis', cloud_command_files=True, result_count=count)
                self.assertEqual(process.returncode, 0, process.stderr)
                self.assertEqual(list_text, ''.join(path + '\n' for path in mapping.values()))
                self.assertNotIn('gs://', list_text)
                self.assertFalse(any('gs://' in arg for arg in args))

    def test_merge_rejects_an_empty_result_array_before_running_the_script(self):
        process, args, _, _ = self.run_command('merge_cis', cloud_command_files=True, result_count=0)
        self.assertNotEqual(process.returncode, 0)
        self.assertIsNone(args)
        self.assertRegex((process.stdout + process.stderr).lower(), 'no.*result')

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
