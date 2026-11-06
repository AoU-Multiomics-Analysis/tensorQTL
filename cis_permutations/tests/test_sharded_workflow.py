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
        struct = doc.struct_typedefs['CisShard'].members
        for key in ['plink_pgen', 'plink_pvar', 'plink_psam', 'phenotype_bed', 'phenotype_groups']:
            self.assertIsInstance(struct[key], WDL.Type.File)
        merge = next(t for t in doc.tasks if t.name == 'merge_cis')
        types = {b.name: b.value.type for b in merge.available_inputs}
        self.assertIsInstance(types['results'].item_type, WDL.Type.File)

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
