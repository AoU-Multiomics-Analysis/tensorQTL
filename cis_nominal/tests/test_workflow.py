"""Validate runtime and execute the rendered command with localized file paths."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import WDL

ROOT = Path(__file__).resolve().parents[1]
WDL_FILE = 'tensorqtl_cis_nominal.wdl'
TASK_NAME = 'tensorqtl_cis_nominal'

class WorkflowTests(unittest.TestCase):
    def runtime_values(self, **overrides):
        doc = WDL.load(str(ROOT / WDL_FILE))
        task = doc.tasks[0]
        inputs = dict(prefix='test', memory=64, disk_space=100,
                      num_threads=16, num_gpus=1, num_preempt=0,
                      plink_pgen='gs://test/input.pgen',
                      plink_pvar='gs://test/input.pvar',
                      plink_psam='gs://test/input.psam',
                      phenotype_bed='gs://test/phenotypes.bed',
                      covariates='gs://test/covariates.tsv')
        inputs.update(overrides)
        try:
            env = WDL.values_from_json(inputs, task.available_inputs)
        except WDL.Error.InputError as error:
            self.fail(str(error))
        stdlib = WDL.StdLib.Base(task.effective_wdl_version)
        for decl in task.inputs or task.postinputs:
            if decl.name not in [binding.name for binding in env]:
                env = env.bind(decl.name, decl.expr.eval(env, stdlib)
                               if decl.expr else WDL.Value.Null())
        return {name: expr.eval(env, stdlib).json
                for name, expr in task.runtime.items()}

    def test_default_l4_runtime(self):
        doc = WDL.load(str(ROOT / WDL_FILE))
        self.assertEqual(doc.wdl_version, '1.0')
        values = self.runtime_values()
        self.assertEqual(values['predefinedMachineType'], 'g2-standard-16')
        self.assertEqual(values['gpuType'], 'nvidia-l4')
        self.assertEqual(values['gpuCount'], 1)
        self.assertEqual(values['zones'], ['us-central1-c'])

    def test_override_machine_gpu_and_zones(self):
        values = self.runtime_values(machine_type='a2-highgpu-1g',
                                     gpu_type='nvidia-tesla-a100',
                                     zones=['us-central1-a', 'us-central1-f'])
        self.assertEqual(values['predefinedMachineType'], 'a2-highgpu-1g')
        self.assertEqual(values['gpuType'], 'nvidia-tesla-a100')
        self.assertEqual(values['zones'], ['us-central1-a', 'us-central1-f'])

    def test_t4_without_predefined_machine(self):
        values = self.runtime_values(machine_type=None,
                                     gpu_type='nvidia-tesla-t4', num_gpus=2)
        self.assertIsNone(values['predefinedMachineType'])
        self.assertEqual(values['gpuType'], 'nvidia-tesla-t4')
        self.assertEqual(values['gpuCount'], 2)

    def test_runtime_inputs_available_in_workflow(self):
        doc = WDL.load(str(ROOT / WDL_FILE))
        for name in ['machine_type', 'gpu_type', 'zones']:
            key = TASK_NAME + '.' + name
            self.assertIn(key, doc.workflow.available_inputs)
            self.assertNotIn(key, doc.workflow.required_inputs)

    def run_task(self, optional=False, cloud=None, localize_cloud=False):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = dict(prefix='test', memory=64, disk_space=100,
                          num_threads=16, num_gpus=1, num_preempt=0)
            keys = ['plink_pgen', 'plink_pvar', 'plink_psam', 'phenotype_bed', 'covariates']
            if optional:
                keys += ['interaction', 'phenotype_groups']
            for key in keys:
                folder = root / (key + " space '$literal")
                folder.mkdir()
                path = folder / key
                path.write_text(key)
                inputs[key] = str(path)
            if cloud:
                inputs[cloud] = 'gs://bucket/input.tsv'
            task = WDL.load(str(ROOT / 'tensorqtl_cis_nominal.wdl')).tasks[0]
            task_inputs = dict(inputs)
            if localize_cloud:
                localized_paths = {'gs://bucket/' + key: inputs[key] for key in keys}
                task_inputs.update({key: 'gs://bucket/' + key for key in keys})
            env = WDL.values_from_json(task_inputs, task.available_inputs)
            if localize_cloud:
                env = WDL.Value.rewrite_env_paths(
                    env, lambda file: localized_paths.get(file.value, file.value))
            stdlib = WDL.StdLib.Base(task.effective_wdl_version)
            for decl in task.inputs:
                if decl.name not in [b.name for b in env]:
                    env = env.bind(decl.name, decl.expr.eval(env, stdlib) if decl.expr else WDL.Value.Null())
            command = task.command.eval(env, stdlib).value
            package = root / 'tensorqtl'
            package.mkdir()
            (package / '__init__.py').write_text('')
            (package / '__main__.py').write_text(
                'import sys,json,pathlib\n'
                'for ext in ["pgen","pvar","psam"]:\n'
                ' assert pathlib.Path(sys.argv[1]+"."+ext).read_text()=="plink_"+ext\n'
                'pathlib.Path("argv.json").write_text(json.dumps(sys.argv[1:]))\n')
            result = subprocess.run(['bash', '-c', command], cwd=root, text=True,
                capture_output=True, env={**os.environ, 'PYTHONPATH':str(root)})
            args = json.loads((root/'argv.json').read_text()) if (root/'argv.json').exists() else None
            return result, args, inputs

    def test_localized_files_and_absent_options(self):
        result, args, inputs = self.run_task()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(args[1], inputs['phenotype_bed'])
        self.assertIn('cis_nominal', args)
        self.assertNotIn('--interaction', args)
        self.assertNotIn('--phenotype_groups', args)

    def test_optional_localized_files(self):
        result, args, inputs = self.run_task(optional=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        for key in ['interaction', 'phenotype_groups']:
            self.assertEqual(args[args.index('--'+key)+1], inputs[key])

    def test_cloud_files_are_rewritten_before_command_rendering(self):
        result, args, inputs = self.run_task(optional=True, localize_cloud=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(args[1], inputs['phenotype_bed'])
        self.assertEqual(args[args.index('--covariates')+1], inputs['covariates'])
        self.assertEqual(args[args.index('--phenotype_groups')+1], inputs['phenotype_groups'])
        self.assertFalse(any('gs://' in arg for arg in args))

    def test_reject_unlocalized_uri(self):
        for key in ['plink_pgen', 'plink_pvar', 'plink_psam', 'phenotype_bed',
                    'covariates', 'interaction', 'phenotype_groups']:
            with self.subTest(key=key):
                result, args, _ = self.run_task(optional=True, cloud=key)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('localization', result.stderr.lower())
                self.assertIsNone(args)

    def test_no_workflow_scope_file_writes(self):
        doc = WDL.load(str(ROOT / 'tensorqtl_cis_nominal.wdl'))
        def visit(node):
            if isinstance(node, WDL.Expr.Apply):
                self.assertNotIn(node.function_name, {'write_lines','write_tsv','write_map','write_json'})
            for child in node.children:
                visit(child)
        visit(doc.workflow)

if __name__ == '__main__':
    unittest.main()
