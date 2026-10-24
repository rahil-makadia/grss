"""Offline Python API coverage and helper regression tests.

The real Python modules/package initializers run with a native-module stand-in
and a small observatory-code fixture. No native numerical behavior is mocked
as a successful integration test: these tests validate Python contracts only.
"""
import ast
import builtins
import importlib
import importlib.util
import inspect
import io
import json
import os
from pathlib import Path
import sys
import tarfile
import tempfile
import types
import unittest
from unittest.mock import Mock, patch

import numpy as np

ROOT = Path(__file__).resolve().parents[3]


class ExportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.saved_modules = {k: v for k, v in sys.modules.items()
                             if k == 'grss' or k.startswith('grss.')}
        for name in cls.saved_modules:
            del sys.modules[name]
        package = types.ModuleType('grss')
        package.__path__ = [str(ROOT/'grss')]
        sys.modules['grss'] = package
        native = types.ModuleType('grss.libgrss')
        native.__all__ = []
        package.libgrss = native
        sys.modules[native.__name__] = native
        real_open = builtins.open

        def fixture_open(file, *args, **kwargs):
            if str(file).endswith('/fit/special_codes.json'):
                return io.StringIO(json.dumps(dict(
                    gaia=['258'], occultation=[], spacecraft=[], roving=[])))
            return real_open(file, *args, **kwargs)

        with patch('builtins.open', side_effect=fixture_open), patch(
                'requests.sessions.Session.request',
                side_effect=AssertionError('Unexpected network request during import')):
            cls.fit = importlib.import_module('grss.fit')
            cls.prop = importlib.import_module('grss.prop')
            cls.utils = importlib.import_module('grss.utils')
            cls.kernels = importlib.import_module('grss.kernels')
            cls.debias = importlib.import_module('grss.debias')

    @classmethod
    def tearDownClass(cls):
        for name in list(sys.modules):
            if name == 'grss' or name.startswith('grss.'):
                del sys.modules[name]
        sys.modules.update(cls.saved_modules)

    def test_every_defined_function_is_exported(self):
        for path in (ROOT/'grss').rglob('*.py'):
            if path.name == '__init__.py':
                continue
            module_name = '.'.join(path.relative_to(ROOT).with_suffix('').parts)
            module = importlib.import_module(module_name)
            tree = ast.parse(path.read_text())
            package = sys.modules[module_name.rsplit('.', 1)[0]]
            for node in tree.body:
                if not isinstance(node, (ast.FunctionDef, ast.ClassDef)):
                    continue
                with self.subTest(module=module_name, helper=node.name):
                    self.assertIn(node.name, module.__all__)
                    self.assertTrue(callable(getattr(module, node.name)))
                    if package is not sys.modules['grss']:
                        self.assertIn(node.name, package.__all__)
                        self.assertIs(getattr(module, node.name), getattr(package, node.name))

    def test_real_root_initializer_exports_packages(self):
        spec = importlib.util.spec_from_file_location(
            'grss', ROOT/'grss/__init__.py', submodule_search_locations=[str(ROOT/'grss')])
        package = importlib.util.module_from_spec(spec)
        with patch.object(self.utils, 'initialize') as initialize:
            spec.loader.exec_module(package)
        initialize.assert_called_once_with()
        for name in ('fit', 'prop', 'utils', 'kernels', 'debias'):
            self.assertIn(name, package.__all__)
            self.assertIs(getattr(package, name), sys.modules[f'grss.{name}'])

    def test_initializer_uses_current_python_and_unsplit_paths(self):
        with patch.object(self.utils, '_connected_to_internet', return_value=True), \
                patch.object(self.utils, 'check_and_get_codes_file'), \
                patch.object(self.utils, 'grss_path', '/temporary/path with spaces/grss'), \
                patch.object(self.utils.subprocess, 'run') as run, patch.dict(os.environ):
            self.utils.initialize()
        self.assertEqual(run.call_args_list[0].args[0],
                         [sys.executable, '/temporary/path with spaces/grss/debias/get_debiasing_data.py'])
        self.assertEqual(run.call_args_list[1].args[0],
                         [sys.executable, '/temporary/path with spaces/grss/kernels/get_kernels.py'])

    def test_private_methods_have_public_callable_aliases(self):
        for cls in (self.fit.FitSimulation, self.fit.IterationParams, self.prop.SigmaPoints):
            for name, method in vars(cls).items():
                if name.startswith('_') and not name.startswith('__'):
                    with self.subTest(cls=cls.__name__, method=name):
                        public = [key for key, value in vars(cls).items()
                                  if not key.startswith('_') and value is method]
                        self.assertTrue(public, f'{cls.__name__}.{name} has no public alias')

    def test_callable_catalog_covers_every_export_and_method(self):
        catalog = (ROOT/'docs/source/python_helpers.rst').read_text()
        for path in (ROOT/'grss').rglob('*.py'):
            if path.name == '__init__.py':
                continue
            module_name = '.'.join(path.relative_to(ROOT).with_suffix('').parts)
            module = importlib.import_module(module_name)
            for name in module.__all__:
                helper = getattr(module, name)
                if callable(helper):
                    self.assertIn(f'{module_name}.{name}', catalog)
                if inspect.isclass(helper):
                    for method_name, method in vars(helper).items():
                        if inspect.isfunction(method):
                            self.assertIn(f'{module_name}.{name}.{method_name}', catalog)

    def test_module_aliases_are_identical_and_exported(self):
        for package in (self.fit, self.prop, self.utils):
            for name in package.__all__:
                helper = getattr(package, name)
                if name.startswith('_') and callable(helper):
                    with self.subTest(package=package.__name__, helper=name):
                        public = [key for key in package.__all__ if not key.startswith('_')
                                  and getattr(package, key) is helper]
                        self.assertTrue(public)

    def test_all_exports_exist_without_duplicates_or_dependency_leaks(self):
        for module in (self.fit, self.prop, self.utils, self.kernels, self.debias):
            self.assertEqual(len(module.__all__), len(set(module.__all__)))
            for name in module.__all__:
                self.assertTrue(hasattr(module, name), (module.__name__, name))
            for dependency in ('np', 'pd', 'plt', 'requests', 'os', 'warnings'):
                self.assertNotIn(dependency, module.__all__)
            namespace = {}
            exec(f'from {module.__name__} import *', namespace)
            self.assertTrue(set(module.__all__).issubset(namespace))

    def test_legacy_fitting_names_and_signatures(self):
        aliases = {
            '_check_priors': 'prepare_priors',
            '_x_dict_to_state': 'solution_to_state',
            '_x_dict_to_nongrav_params': 'solution_to_nongrav_params',
            '_x_dict_to_events': 'solution_to_events',
            '_get_perturbed_state': 'get_perturbed_state',
            '_get_perturbation_info': 'get_perturbation_info',
            '_get_residuals_and_partials': 'compute_residuals_and_partials',
            '_get_rms_and_reject_outliers': 'compute_fit_statistics',
            '_get_lsq_state_correction': 'solve_state_correction',
        }
        for old, public in aliases.items():
            self.assertIs(getattr(self.fit.FitSimulation, old),
                          getattr(self.fit.FitSimulation, public))
        fitter = object.__new__(self.fit.FitSimulation)
        values = [[1., np.nan], [2., 3.]]
        np.testing.assert_array_equal(fitter._flatten_and_clean(values), [1, 2, 3])
        np.testing.assert_array_equal(
            self.fit.FitSimulation._flatten_and_clean(fitter, values), [1, 2, 3])
        self.assertEqual(list(inspect.signature(fitter.flatten_and_clean).parameters), ['arr'])

    def test_sigma_dictionaries_preserve_parameter_order(self):
        for state in ({'x': 2., 'y': 4.}, {'t': 60000., 'x': 2., 'y': 4.},
                      {'x': 2., 't': 60000., 'y': 4.}):
            points = self.prop.SigmaPoints(state, np.diag([1., 2.]), 'julier', kappa=1.)
            for i, converted in enumerate(points.get_points_dict()):
                for column, key in enumerate(key for key in state if key != 't'):
                    self.assertEqual(converted[key], points.sigma_points[i, column])
                if 't' in state:
                    self.assertEqual(converted['t'], state['t'])
            mean, covariance = points.reconstruct(points.sigma_points)
            np.testing.assert_allclose(mean, [2, 4])
            np.testing.assert_allclose(covariance, np.diag([1, 2]))
            np.testing.assert_allclose(points.get_weights(), points._get_weights())

    def test_merwe_weight_api(self):
        points = self.prop.SigmaPoints({'x': 2.}, np.eye(1), 'merwe',
                                       alpha=.5, beta=2., kappa=0.)
        np.testing.assert_allclose(points.merwe_weights(), points.get_weights())
        mean, covariance = points.reconstruct(points.sigma_points)
        np.testing.assert_allclose(mean, [2])
        np.testing.assert_allclose(covariance, [[1]])

    def test_log_parser_accepts_marker_on_first_line(self):
        lines = ['$$CA_START', 't | flybyBody', '60000 | asteroid', '$$CA_END']
        parallel = importlib.import_module('grss.prop.prop_parallel')
        with patch.object(parallel, '_reconstruct_ca_params', side_effect=lambda row: row):
            approaches, impacts = self.prop.reconstruct_ca_and_impact(lines)
        self.assertEqual(approaches, [{'t': 60000., 'flybyBody': 'asteroid'}])
        self.assertEqual(impacts, [])

    def test_download_helpers_do_not_write_http_errors(self):
        import requests
        for module_name in ('grss.kernels.get_kernels', 'grss.debias.get_debiasing_data'):
            module = importlib.import_module(module_name)
            with tempfile.TemporaryDirectory() as directory:
                output = Path(directory)/'data'
                output.write_bytes(b'original')
                kwargs = {'force': True} if module is sys.modules['grss.kernels.get_kernels'] else {}
                with patch.object(module.requests, 'head', return_value=Mock(status_code=200)), \
                        patch.object(module.requests, 'get', return_value=Mock(
                            raise_for_status=Mock(side_effect=requests.HTTPError('404')))):
                    with self.assertRaises(ConnectionError):
                        module.check_and_download('https://example.invalid/data', output, 'test', **kwargs)
                self.assertEqual(output.read_bytes(), b'original')

    def test_existing_kernel_is_skipped(self):
        module = importlib.import_module('grss.kernels.get_kernels')
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)/'existing'
            output.write_bytes(b'cached')
            with patch.object(module.requests, 'head') as head:
                module.check_and_download('https://example.invalid/data', output, 'test')
                head.assert_not_called()

    def test_safe_archive_extracts_paths_with_spaces(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive_path = root/'tables.tgz'
            with tarfile.open(archive_path, 'w:gz') as archive:
                info = tarfile.TarInfo('nested/table file.txt')
                data = b'bias values'
                info.size = len(data)
                archive.addfile(info, io.BytesIO(data))
            destination = root/'tables with spaces'
            destination.mkdir()
            (destination/'unrelated').write_text('keep')
            self.debias.extract_debiasing_archive(archive_path, destination)
            self.assertEqual((destination/'nested/table file.txt').read_bytes(), data)
            self.assertTrue((destination/'unrelated').exists())

    def test_archive_rejects_traversal_links_and_symlink_destinations(self):
        for name, entry_type in [('../escaped', tarfile.REGTYPE),
                                 ('/absolute', tarfile.REGTYPE),
                                 ('linked', tarfile.SYMTYPE)]:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                archive_path = root/'tables.tgz'
                with tarfile.open(archive_path, 'w:gz') as archive:
                    info = tarfile.TarInfo(name)
                    info.type = entry_type
                    info.linkname = '../escaped'
                    archive.addfile(info, io.BytesIO(b''))
                with self.assertRaises(ValueError):
                    self.debias.extract_debiasing_archive(archive_path, root/'tables')
                self.assertFalse((root/'escaped').exists())
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            outside, destination = root/'outside', root/'tables'
            outside.mkdir()
            destination.mkdir()
            (destination/'nested').symlink_to(outside, target_is_directory=True)
            archive_path = root/'tables.tgz'
            with tarfile.open(archive_path, 'w:gz') as archive:
                info = tarfile.TarInfo('nested/file')
                archive.addfile(info, io.BytesIO(b''))
            with self.assertRaises(ValueError):
                self.debias.extract_debiasing_archive(archive_path, destination)
            self.assertEqual(list(outside.iterdir()), [])


if __name__ == '__main__':
    unittest.main()
