"""Offline IOD tests; analytic Kepler orbits supply independent reference data.

Run with ``python -m unittest discover -s tests/python/iod -v``. Loading the
numerical module directly avoids GRSS's package-level kernel downloads.
"""
import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[3] / 'grss' / 'fit' / 'fit_iod.py'
SPEC = importlib.util.spec_from_file_location('fit_iod', MODULE_PATH)
iod = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = iod
SPEC.loader.exec_module(iod)


def ellipse(times, a, eccentricity, mean0, mu, rotation):
    """Independent elliptic-anomaly solution (does not use the IOD propagator)."""
    n = np.sqrt(mu/a**3)
    anomaly = mean0+n*np.asarray(times)
    for _ in range(20):
        anomaly -= (anomaly-eccentricity*np.sin(anomaly)-mean0-n*np.asarray(times)) / (
            1-eccentricity*np.cos(anomaly))
    c, s = np.cos(anomaly), np.sin(anomaly)
    q = np.sqrt(1-eccentricity**2)
    positions = a*np.column_stack((c-eccentricity, q*s, np.zeros_like(c)))
    velocities = a*n*np.column_stack((-s, q*c, np.zeros_like(c))) / (
        1-eccentricity*c)[:, None]
    return positions @ rotation.T, velocities @ rotation.T


class GaussTests(unittest.TestCase):
    def setUp(self):
        self.mu = 0.01720209895**2
        self.times = np.array([-4., 0., 5.])
        i, w = 0.3, 0.7
        self.rotation = np.array([[np.cos(w), -np.sin(w)*np.cos(i), np.sin(w)*np.sin(i)],
                                  [np.sin(w), np.cos(w)*np.cos(i), -np.cos(w)*np.sin(i)],
                                  [0, np.sin(i), np.cos(i)]])
        self.observers = ellipse(self.times, 1., 0.0167, 0.2, self.mu, np.eye(3))[0]
        self.positions, self.velocities = ellipse(
            self.times, 2.2, 0.25, 0.5, self.mu, self.rotation)
        directions = self.positions-self.observers
        directions /= np.linalg.norm(directions, axis=1)[:, None]
        self.ra = np.arctan2(directions[:, 1], directions[:, 0])
        self.dec = np.arcsin(directions[:, 2])

    def solve(self, **kwargs):
        args = dict(times=self.times+60000, ra=self.ra, dec=self.dec,
                    observer_positions=self.observers, mu=self.mu)
        args.update(kwargs)
        return iod.gauss_iod(**args)

    def best(self, solutions):
        truth = np.r_[self.positions[1], self.velocities[1]]
        return min(solutions, key=lambda sol: np.linalg.norm(sol.state-truth))

    def test_heliocentric_recovery(self):
        solution = self.best(self.solve())
        np.testing.assert_allclose(solution.state[:3], self.positions[1], atol=2e-7)
        np.testing.assert_allclose(solution.state[3:], self.velocities[1], atol=2e-9)
        np.testing.assert_allclose(solution.ranges, np.linalg.norm(
            self.positions-self.observers, axis=1), atol=2e-7)
        self.assertLessEqual(solution.residual, 1e-10)
        self.assertGreater(solution.iterations, 0)
        self.assertEqual(solution.epoch, 60000)

    def test_unit_invariance(self):
        au, day = 149597870.7, 86400.
        au_solution = self.best(self.solve())
        km_solutions = self.solve(times=self.times*day,
                                 observer_positions=self.observers*au,
                                 mu=self.mu*au**3/day**2)
        truth = np.r_[self.positions[1]*au, self.velocities[1]*au/day]
        solution = min(km_solutions, key=lambda sol: np.linalg.norm(sol.state-truth))
        np.testing.assert_allclose(solution.state[:3]/au, au_solution.state[:3], atol=1e-8)
        np.testing.assert_allclose(solution.state[3:]*day/au, au_solution.state[3:], atol=1e-10)

    def test_ra_wrap_and_input_not_modified(self):
        original = [x.copy() for x in (self.times, self.ra, self.dec, self.observers)]
        a, b = self.best(self.solve()), self.best(self.solve(ra=self.ra+2*np.pi))
        np.testing.assert_allclose(a.state, b.state, atol=1e-8)
        for before, after in zip(original, (self.times, self.ra, self.dec, self.observers)):
            np.testing.assert_array_equal(before, after)

    def test_barycentric_fit_conversion(self):
        solution = self.best(self.solve())
        origin = np.array([0.01, 0.02, 0.03, 0.0001, 0.0002, 0.0003])
        converted = solution.to_fit_dict(origin)
        self.assertEqual(list(converted), ['t', 'x', 'y', 'z', 'vx', 'vy', 'vz'])
        np.testing.assert_allclose(list(converted.values())[1:], solution.state+origin)
        self.assertEqual(converted['t'], 60000)
        for invalid in ([0, 0, 0], [np.nan]*6):
            with self.assertRaises(ValueError):
                solution.to_fit_dict(invalid)

    def test_invalid_inputs(self):
        for kwargs in [dict(times=[0, 0, 1]), dict(times=[2, 1, 0]),
                       dict(times=[0, 1]), dict(ra=[0, 1]),
                       dict(dec=[0, 0, 2]), dict(ra=[0, np.nan, 1]),
                       dict(observer_positions=np.zeros((3, 2))),
                       dict(observer_positions=np.zeros((3, 3))),
                       dict(observer_positions=np.full((3, 3), np.inf)),
                       dict(mu=0), dict(mu=-1), dict(mu=np.nan), dict(mu=[1]),
                       dict(tol=0), dict(tol=np.inf), dict(tol=2),
                       dict(max_iter=-1), dict(max_iter=1.5), dict(max_iter=True),
                       dict(ra=[0, 0, 0], dec=[0, 0, 0]),
                       dict(ra=[0, 0.1, 0.2], dec=[0, 0, 0])]:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                self.solve(**kwargs)

    def test_iteration_limit(self):
        with self.assertRaisesRegex(ValueError, 'No physical Gauss candidate'):
            self.solve(max_iter=0)

    def test_multiple_physical_candidates(self):
        i, w = 0.1717079510417686, -0.2976843726265972
        rotation = np.array([[np.cos(w), -np.sin(w)*np.cos(i), np.sin(w)*np.sin(i)],
                             [np.sin(w), np.cos(w)*np.cos(i), -np.cos(w)*np.sin(i)],
                             [0, np.sin(i), np.cos(i)]])
        pos, vel = ellipse(self.times, 3.929304819746592, 0.5066838211942117,
                           1.7163858316617229, self.mu, rotation)
        directions = pos-self.observers
        directions /= np.linalg.norm(directions, axis=1)[:, None]
        solutions = self.solve(ra=np.arctan2(directions[:, 1], directions[:, 0]),
                               dec=np.arcsin(directions[:, 2]))
        self.assertEqual(len(solutions), 2)
        self.assertLess(solutions[0].root, solutions[1].root)
        self.assertGreater(np.linalg.norm(solutions[0].state-solutions[1].state), 1.)
        self.assertTrue(all(sol.residual <= 1e-10 for sol in solutions))
        np.testing.assert_allclose(solutions[-1].state, np.r_[pos[1], vel[1]], atol=1e-7)

    def test_geocentric_satellite(self):
        mu = 398600.4418
        times = np.array([-60., 0., 80.])
        theta = times*7.292115e-5
        observers = 6378.137*np.column_stack((np.cos(theta), np.sin(theta), np.zeros(3)))
        rotation = np.array([[1, 0, 0], [0, np.cos(.5), -np.sin(.5)],
                             [0, np.sin(.5), np.cos(.5)]])
        pos, vel = ellipse(times, 12000., .2, .4, mu, rotation)
        directions = pos-observers
        directions /= np.linalg.norm(directions, axis=1)[:, None]
        solutions = iod.gauss_iod(times, np.arctan2(directions[:, 1], directions[:, 0]),
                                  np.arcsin(directions[:, 2]), observers, mu)
        solution = min(solutions, key=lambda s: np.linalg.norm(s.state[:3]-pos[1]))
        np.testing.assert_allclose(solution.state[:3], pos[1], atol=1e-4, rtol=0)
        np.testing.assert_allclose(solution.state[3:], vel[1], atol=1e-7, rtol=0)

    def test_hyperbolic_orbit_recovery(self):
        a, e, h0 = 2., 1.4, .4
        mean = e*np.sinh(h0)-h0+np.sqrt(self.mu/a**3)*self.times
        h = np.full(3, h0)
        for _ in range(20):
            h -= (e*np.sinh(h)-h-mean)/(e*np.cosh(h)-1)
        rotation = np.array([[1, 0, 0], [0, np.cos(.3), -np.sin(.3)],
                             [0, np.sin(.3), np.cos(.3)]])
        pos = a*np.column_stack((e-np.cosh(h), np.sqrt(e*e-1)*np.sinh(h), np.zeros(3)))
        vel = np.sqrt(self.mu/a)*np.column_stack(
            (-np.sinh(h), np.sqrt(e*e-1)*np.cosh(h), np.zeros(3))) / (e*np.cosh(h)-1)[:, None]
        pos, vel = pos @ rotation.T, vel @ rotation.T
        directions = pos-self.observers
        directions /= np.linalg.norm(directions, axis=1)[:, None]
        solutions = self.solve(ra=np.arctan2(directions[:, 1], directions[:, 0]),
                               dec=np.arcsin(directions[:, 2]))
        solution = min(solutions, key=lambda s: np.linalg.norm(s.state[:3]-pos[1]))
        np.testing.assert_allclose(solution.state, np.r_[pos[1], vel[1]], atol=2e-8, rtol=0)

    def test_universal_propagator_elliptic_and_near_parabolic(self):
        for eccentricity in [0., 0.6, 0.9999]:
            times = np.array([-100., 0., 100.])
            pos, vel = ellipse(times, 2.2, eccentricity, 0.5, self.mu, self.rotation)
            for j in (0, 2):
                result = iod._propagate(np.r_[pos[1], vel[1]], times[j], self.mu)
                np.testing.assert_allclose(result, pos[j], atol=1e-11)

    def test_universal_propagator_hyperbolic(self):
        # Analytic hyperbolic anomaly at periapsis; a denotes |semi-major axis|.
        a, e, h = 2., 1.5, 0.6
        t = np.sqrt(a**3/self.mu)*(e*np.sinh(h)-h)
        state = np.array([a*(e-1), 0, 0, 0, np.sqrt(self.mu/a*(e+1)/(e-1)), 0])
        for sign in [-1, 1]:
            expected = [a*(e-np.cosh(h)), sign*a*np.sqrt(e*e-1)*np.sinh(h), 0]
            np.testing.assert_allclose(iod._propagate(state, sign*t, self.mu), expected, atol=1e-11)


if __name__ == '__main__':
    unittest.main()
