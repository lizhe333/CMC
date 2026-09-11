"""Named numerical contracts, using synthetic verification cases (not competition data)."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from openpyxl import Workbook
from scipy.optimize import brentq
from scipy.special import j0, j1

from q1_solver import ALPHA, H, K, R, read_boundary, solve, solve_bdf, stability


class Q1Contracts(unittest.TestCase):
    def test_uniform_equilibrium_and_axis(self):
        """No gradients or boundary driving implies exactly unchanged fields."""
        boundary = np.array([[0., 28., 2.55], [100., 28., 2.55]])
        t, c, balance, _ = solve(boundary, end=100)
        np.testing.assert_array_equal(t, 28.)
        np.testing.assert_array_equal(c, 2.55)
        self.assertLess(abs(balance[:, 2]).max(), 1e-13)

    def test_variable_diffusivity_mass_balance_and_independent_time_solvers(self):
        """Conservative nonlinear diffusion and the BDF implementation agree with explicit."""
        boundary = np.array([[0., 28., .02], [60., 32., .03], [180., 40., .04]])
        reference = solve_bdf(boundary, n=20, rtol=1e-11, end=180)
        explicit = solve(boundary, n=20, substeps=128, end=180)
        self.assertLess(np.max(abs(reference[0] - explicit[0])), 5e-5)
        self.assertLess(np.max(abs(reference[1] - explicit[1])), 3e-5)
        self.assertLess(abs(explicit[2][:, 2]).max(), 1e-12)
        self.assertLess(abs(reference[2][:, 2]).max(), 1e-12)

    def test_robin_cylinder_against_bessel_series(self):
        """Independent analytic heat solution verifies curvature, center and surface signs."""
        boundary = np.array([[0., 50., 2.55], [1000., 50., 2.55]])
        field = solve_bdf(boundary, n=160, end=1000)[0]
        bi = H * R / K
        equation = lambda x: x * j1(x) - bi * j0(x)
        grid = np.linspace(.001, 100., 2001)
        roots = np.array([brentq(equation, a, b) for a, b in zip(grid[:-1], grid[1:])
                          if equation(a) * equation(b) < 0])
        coefficients = 2 * j1(roots) / (roots * (j0(roots)**2 + j1(roots)**2))
        radius = np.linspace(0, R, 21)
        for second in (100, 300, 600, 1000):
            theta = np.sum(coefficients[:, None] * j0(roots[:, None] * radius / R) *
                           np.exp(-roots[:, None]**2 * ALPHA * second / R**2), axis=0)
            exact = 50. - 22. * theta
            self.assertLess(np.max(abs(field[second] - exact)), 4e-4)

    def test_reject_unstable_explicit_step(self):
        with self.assertRaisesRegex(ValueError, "nonnegative weights"):
            stability(40, 1)

    def test_reject_duplicate_times_and_missing_coverage(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.xlsx"
            for rows in (
                [(0, 28, .02), (0, 29, .03), (1800, 40, .04)],
                [(0, 28, .02), (60, 29, .03)],
            ):
                wb = Workbook()
                ws = wb.active
                ws.append(["时间", "温度", "水分浓度"])
                for row in rows:
                    ws.append(row)
                wb.save(path)
                wb.close()
                with self.assertRaises(ValueError):
                    read_boundary(path)


if __name__ == "__main__":
    unittest.main(verbosity=2)
