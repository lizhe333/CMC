"""Targeted contracts for the question 3 endpoint pipeline.

These tests cover the boundary decision, event/output definitions, and the
read-back contract.  They do not substitute for the Q3 spatial/time
refinement runs.
"""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

CODE_DIR = Path(__file__).resolve().parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

import q2_solver as q2
import q3_solver as q3


class Q3Contracts(unittest.TestCase):
    def test_final_hour_platform_diagnostic(self):
        boundary, _ = q2.read_boundary(q3.DEFAULT_BOUNDARY, end=14400)
        diagnostic = q3.platform_diagnostic(boundary)
        self.assertTrue(diagnostic["pass"])
        self.assertEqual(diagnostic["sample_count"], 61)
        self.assertAlmostEqual(diagnostic["temperature_C"]["mean"], 49.998934, places=5)
        self.assertAlmostEqual(diagnostic["moisture_kg_per_kg"]["mean"], 0.04998754, places=7)
        plateau = q3.make_plateau_boundary(diagnostic, 3600 * 24)
        self.assertEqual(plateau.shape, (2, 3))
        self.assertAlmostEqual(float(plateau[0, 1]), 50.0)
        self.assertAlmostEqual(float(plateau[0, 2]), 0.05)

    def test_platform_failure_blocks_extension(self):
        # Keep the required 61 one-minute samples, but introduce a clear
        # temperature trend so the test exercises the platform gate rather
        # than the input-length guard.
        times = np.linspace(10800.0, 14400.0, 61)
        boundary = np.column_stack([times, np.linspace(49.0, 55.0, 61), np.full(61, 0.05)])
        diagnostic = q3.platform_diagnostic(boundary)
        self.assertFalse(diagnostic["pass"])
        with self.assertRaises(ValueError):
            q3.make_plateau_boundary(diagnostic, 20000)

    def test_full_grid_monotonicity_metrics(self):
        n = 20
        _, nodes, _, volumes = q2.radial_geometry(n)
        z = np.empty((2 * (n + 1), 2))
        z[0::2, :] = 28.0
        z[1::2, :] = np.log(np.column_stack([np.linspace(2.0, 1.0, n + 1), np.linspace(1.5, 0.5, n + 1)]))
        sampled_t, sampled_c, audit, physical = q3._sample_state(z, n, volumes, np.arange(0, n + 1))
        self.assertEqual(sampled_c.shape, (2, n + 1))
        self.assertLessEqual(float(np.max(audit[:, 3])), 0.0)
        self.assertLessEqual(float(np.max(audit[:, 2])), 0.0)
        self.assertAlmostEqual(float(np.max(audit[:, 1])), 2.0)

    def test_workbook_is_exact_minute_schema(self):
        times = np.array([60, 120, 180], dtype=int)
        moisture = np.arange(63, dtype=float).reshape(3, 21) / 100.0
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result3.xlsx"
            q3.write_result3_workbook(path, times, moisture, template=q3.DEFAULT_TEMPLATE)
            report = q3.validate_result3_workbook(path, times, moisture, template=q3.DEFAULT_TEMPLATE)
            self.assertTrue(report["four_decimal_match"])
            self.assertEqual(report["rows"], 3)
            with self.assertRaises(ValueError):
                q3.write_result3_workbook(Path(directory) / "bad.xlsx", np.array([60, 90]), np.zeros((2, 21)), template=q3.DEFAULT_TEMPLATE)

    def test_event_definition_uses_first_full_minute_only_for_workbook(self):
        result = {
            "event": {"time_s": 3612.4, "state": np.zeros(42), "center_C": 0.15, "max_C": 0.15, "center_minus_max_C": 0.0},
            "second_segment": {"end_s_requested": 10000},
        }
        meta = q3._event_metadata(result)
        self.assertEqual(meta["status"], "EVENT_LOCATED")
        self.assertAlmostEqual(meta["time_h"], 3612.4 / 3600.0)

    def test_formal_publish_gate_evidence(self):
        validation_path = q3.DEFAULT_OUTPUT / "q3_validation.json"
        self.assertTrue(validation_path.exists())
        validation = json.loads(validation_path.read_text(encoding="utf-8"))
        self.assertEqual(validation["status"], "NUMERICAL_CHECKS_PASSED")
        self.assertTrue(all(validation["gate_summary"].values()))
        self.assertLessEqual(
            abs(float(validation["event"]["event_residual_C"])),
            float(validation["strict_event_evidence"]["event_residual_tolerance_C"]),
        )
        self.assertTrue(validation["event"]["strict_before_pass"])
        self.assertTrue(validation["event"]["strict_after_pass"])
        selected = validation["convergence"]["selected_final_settings"]
        self.assertEqual(int(selected["n_intervals"]), int(validation["n_intervals"]))
        self.assertEqual(float(selected["rtol"]), float(validation["rtol"]))
        self.assertEqual(float(selected["atol"]), float(validation["atol"]))
        self.assertTrue(validation["official_workbook_validation"]["four_decimal_match"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
