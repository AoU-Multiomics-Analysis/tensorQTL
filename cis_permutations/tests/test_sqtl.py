"""Check the sQTL entry point and its required, localized groups file."""
from pathlib import Path
import unittest

import WDL
import test_workflow as helpers

ROOT = Path(__file__).resolve().parents[1]
DESCRIPTOR = 'workflows/tensorqtl_cis_permutations_sQTL.wdl'


class SqtlTests(unittest.TestCase):
    def load_workflow(self):
        try:
            return WDL.load(str(ROOT / DESCRIPTOR))
        except WDL.Error.MultipleDefinitions as error:
            self.fail(f'sQTL descriptor has duplicate inputs: {error}')

    def test_groups_file_is_required_at_workflow_input(self):
        doc = self.load_workflow()
        self.assertEqual(doc.wdl_version, '1.0')
        required = {binding.name: binding.value.type
                    for binding in doc.workflow.required_inputs}
        self.assertIn('phenotype_groups', required)
        self.assertIsInstance(required['phenotype_groups'], WDL.Type.File)
        self.assertFalse(required['phenotype_groups'].optional)
        call = doc.workflow.body[0]
        values = WDL.values_from_json({'phenotype_groups': 'gs://test/groups.tsv'},
                                     doc.workflow.available_inputs)
        group = call.inputs['phenotype_groups'].eval(values, WDL.StdLib.Base('1.0'))
        self.assertIsInstance(group, WDL.Value.File)
        self.assertEqual(group.value, 'gs://test/groups.tsv')

    def test_imported_task_localizes_group_file(self):
        self.load_workflow()
        result, args, inputs = helpers.WorkflowTests().run_task(
            descriptor=DESCRIPTOR, optional=True, localize_cloud=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(args[args.index('--phenotype_groups') + 1],
                         inputs['phenotype_groups'])
        self.assertEqual(args[args.index('--mode') + 1], 'cis')
        self.assertFalse(any('gs://' in arg for arg in args))

    def test_imported_task_rejects_unlocalized_groups(self):
        self.load_workflow()
        result, args, _ = helpers.WorkflowTests().run_task(
            descriptor=DESCRIPTOR, optional=True, cloud='phenotype_groups')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('localization', result.stderr.lower())
        self.assertIsNone(args)

    def test_no_workflow_scope_file_writes(self):
        doc = self.load_workflow()

        def visit(node):
            if isinstance(node, WDL.Expr.Apply):
                self.assertNotIn(node.function_name,
                                 {'write_lines', 'write_tsv', 'write_map', 'write_json'})
            for child in node.children:
                visit(child)

        visit(doc.workflow)


if __name__ == '__main__':
    unittest.main()
