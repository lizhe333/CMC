"""Named contracts for Q2 supplementary constitutive and observer interfaces."""
from __future__ import annotations

import sys
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

CODE_DIR = Path(__file__).resolve().parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

import q2_solver as q2
import q2_supplementary_validation as supplementary


class Q2SupplementaryContracts(unittest.TestCase):
    @staticmethod
    def _boundary(end: float = 30.0) -> np.ndarray:
        return np.array(
            [[0.0, 28.0, 2.55], [end / 2.0, 31.0, 2.20], [end, 34.0, 1.90]],
            dtype=float,
        )

    def test_property_model_default_and_q1_constitutive_derivatives(self):
        rho_default = q2.material_properties(2.55, 28.0)
        rho_explicit = q2.material_properties(2.55, 28.0, property_model="q2")
        for first, second in zip(rho_default, rho_explicit):
            np.testing.assert_array_equal(first, second)

        rho, cp, k, diffusivity = q2.material_properties(2.55, 28.0, property_model="q1")
        self.assertAlmostEqual(float(rho), 820.0)
        self.assertAlmostEqual(float(cp), 2600.0)
        self.assertAlmostEqual(float(k), 0.36)
        self.assertAlmostEqual(float(diffusivity), 4.937655094173937e-9, places=18)
        rho_c, cp_c, k_c, D_c, D_t = q2.property_derivatives(
            2.55, 28.0, property_model="q1"
        )
        self.assertEqual(float(rho_c), 0.0)
        self.assertEqual(float(cp_c), 0.0)
        self.assertEqual(float(k_c), 0.0)
        self.assertGreater(float(D_c), 0.0)
        self.assertEqual(float(D_t), 0.0)
        with self.assertRaises(ValueError):
            q2.material_properties(2.55, 28.0, property_model="unknown")

    def test_q1_degradation_removes_q2_cross_dependencies(self):
        n = 20
        boundary = self._boundary()
        state = q2.initial_state(n, property_model="q1")
        state[0::2] = np.linspace(28.0, 35.0, n + 1)
        state[1::2] = np.linspace(2.5, 1.4, n + 1)
        base = q2.assemble_rhs(state, 3.0, boundary, n, property_model="q1")
        temperature_perturbed = state.copy()
        temperature_perturbed[2] += 1.0
        moisture_perturbed = state.copy()
        moisture_perturbed[1] += 0.1
        moisture_change_from_temperature = q2.assemble_rhs(
            temperature_perturbed, 3.0, boundary, n, property_model="q1"
        )[1::2] - base[1::2]
        heat_change_from_moisture = q2.assemble_rhs(
            moisture_perturbed, 3.0, boundary, n, property_model="q1"
        )[0::2] - base[0::2]
        self.assertLess(float(np.max(np.abs(moisture_change_from_temperature))), 1e-18)
        self.assertLess(float(np.max(np.abs(heat_change_from_moisture))), 1e-18)

    def test_segment_observer_dense_nodes_and_weights_are_read_only(self):
        boundary = np.array([[0.0, 28.0, 2.55], [12.0, 28.0, 2.55]], dtype=float)
        plain = q2.integrate_case(boundary, n=20, end=12, chunk_seconds=6)
        observer = q2.SegmentObserver(retain=True)
        observed = q2.integrate_case(
            boundary,
            n=20,
            end=12,
            chunk_seconds=6,
            segment_observer=observer,
        )
        np.testing.assert_array_equal(plain["temperature_C"], observed["temperature_C"])
        np.testing.assert_array_equal(plain["moisture_kg_per_kg"], observed["moisture_kg_per_kg"])
        self.assertEqual(observer.total_segments, 2)
        self.assertEqual(len(observer.segments), 2)
        for segment in observer.segments:
            nodes = segment["accepted_step_nodes"]
            self.assertGreaterEqual(nodes.size, 2)
            self.assertTrue(np.all(np.diff(nodes) > 0.0))
            self.assertEqual(segment["accepted_step_trapezoid_weights"].shape, nodes.shape)
            self.assertAlmostEqual(
                float(segment["accepted_step_trapezoid_weights"].sum()), 6.0, places=10
            )
            self.assertEqual(segment["control_volume_weights_m2"].shape, (21,))
            self.assertFalse(segment["control_volume_weights_m2"].flags.writeable)
            log_state = segment["dense_log_solution"](nodes[0])
            physical_state = segment["dense_solution"](nodes[0])
            np.testing.assert_allclose(physical_state[1::2], np.exp(log_state[1::2]))
            self.assertTrue(np.all(np.isfinite(segment["dense_solution"](nodes[0]))))
            self.assertTrue(
                np.all(np.isfinite(segment["moisture_surface_log_solution"](nodes)))
            )
            self.assertGreater(float(segment["moisture_surface_solution"](nodes[-1])), 0.0)

    def test_q1_degraded_integrate_case_keeps_compatible_result_keys(self):
        result = q2.integrate_case(
            self._boundary(), n=20, end=10, chunk_seconds=5, property_model="q1"
        )
        self.assertEqual(result["property_model"], "q1")
        for key in (
            "times_s",
            "temperature_C",
            "moisture_kg_per_kg",
            "mass_balance_residual",
            "temperature_final_C",
            "moisture_final_kg_per_kg",
        ):
            self.assertIn(key, result)
        self.assertEqual(result["temperature_C"].shape, (11, 21))
        self.assertEqual(result["moisture_kg_per_kg"].shape, (11, 21))
        self.assertTrue(np.isfinite(result["temperature_C"]).all())
        self.assertTrue(np.isfinite(result["moisture_kg_per_kg"]).all())

    def test_balance_union_gl_quadrature_and_signed_B_on_short_segment(self):
        n = 20
        _, radius_nodes, _, control_volume_weights = q2.radial_geometry(n)
        end = 10
        boundary = np.array([[0.0, 28.0, 2.0], [float(end), 28.0, 2.0]])
        log_c = np.log(2.55)

        def dense_log_solution(query):
            values = np.asarray(query, dtype=float)
            if values.ndim == 0:
                state = np.empty(2 * (n + 1), dtype=float)
                state[0::2] = 28.0
                state[1::2] = log_c
                return state
            state = np.empty((2 * (n + 1), values.size), dtype=float)
            state[0::2, :] = 28.0
            state[1::2, :] = log_c
            return state

        payload = {
            "start_s": 0.0,
            "end_s": float(end),
            "accepted_step_nodes": np.array([0.0, 3.0, 7.0, 10.0]),
            "accepted_step_trapezoid_weights": np.array([1.5, 3.5, 3.5, 1.5]),
            "boundary_nodes": np.array([0.0, 5.0, 10.0]),
            "radius_nodes_m": radius_nodes.copy(),
            "control_volume_weights_m2": control_volume_weights.copy(),
            "hm_m_s": 1e-6,
            "property_model": "q2",
            "dense_log_solution": dense_log_solution,
            "moisture_surface_log_solution": lambda query: np.full(
                np.asarray(query, dtype=float).shape, log_c
            ),
        }
        accumulator = supplementary._BalanceAccumulator(boundary, end, orders=(4, 8))
        union = accumulator._union_nodes(payload)
        for value in (3.0, 5.0, 7.0):
            self.assertIn(value, union)
        accumulator.consume(payload)
        result = accumulator.finalize()
        expected_loss = q2.RADIUS_M * 1e-6 * (2.55 - 2.0) * np.arange(end + 1)
        expected_B = expected_loss / supplementary.TOTAL_VOLUME
        np.testing.assert_allclose(result["cumulative_loss_gl4"], expected_loss, atol=1e-14)
        np.testing.assert_allclose(result["cumulative_loss_gl8"], expected_loss, atol=1e-14)
        np.testing.assert_allclose(result["B_gl4_signed"], expected_B, atol=1e-12)
        np.testing.assert_allclose(result["B_gl8_signed"], expected_B, atol=1e-12)
        self.assertLess(
            float(np.max(np.abs(result["B_gl8_signed"] - result["B_gl4_signed"]))),
            1e-13,
        )
        self.assertEqual(result["accepted_time_s"].shape, (4,))
        self.assertEqual(result["accepted_B_gl4_abs"].shape, (4,))
        self.assertAlmostEqual(float(result["accepted_time_s"][-1]), 10.0)

    def test_quadrature_selection_uses_signed_B_floor_and_named_status(self):
        balance = {
            "B_gl4_signed": np.zeros(3),
            "B_gl8_signed": np.full(3, 5e-10),
            "B_gl16_signed": np.full(3, 2e-9),
            "B_gl32_signed": np.full(3, 4e-9),
            "accepted_B_gl4_signed": np.zeros(2),
            "accepted_B_gl8_signed": np.full(2, 4e-10),
            "accepted_B_gl16_signed": np.full(2, 2e-9),
            "accepted_B_gl32_signed": np.full(2, 4e-9),
        }
        selection = supplementary._select_quadrature_orders(balance)
        self.assertEqual(selection["selected_orders"], [4, 8])
        self.assertTrue(selection["pass"])
        self.assertLessEqual(selection["comparisons"][0]["max_abs_signed_B_difference"], 1e-9)
        self.assertAlmostEqual(
            selection["comparisons"][0]["accepted_max_abs_signed_B_difference"], 4e-10
        )

    def test_time_gate_has_independent_temperature_and_moisture_thresholds(self):
        times = np.array([0, 1800], dtype=int)
        radius = np.array([0.0, 2.0])
        first_t = np.zeros((2, 2))
        second_t = first_t.copy()
        second_t[1, 1] = 5e-8
        first_c = np.zeros((2, 2))
        second_c = first_c.copy()
        second_c[1, 1] = 9e-9
        comparison = supplementary._field_comparison(
            first_t,
            second_t,
            first_c,
            second_c,
            times,
            radius,
            temperature_threshold=1e-7,
            moisture_threshold=1e-8,
        )
        self.assertTrue(comparison["pass"])
        self.assertEqual(comparison["temperature_threshold_C"], 1e-7)
        self.assertEqual(comparison["moisture_threshold_kg_per_kg"], 1e-8)
        self.assertNotEqual(supplementary.DEFAULT_OUTPUT.resolve(), q2.DEFAULT_OUTPUT.resolve())

    def test_plot_stage_is_separate_and_does_not_invoke_solver(self):
        plot_script = CODE_DIR / "plot_q2_supplementary_validation.py"
        self.assertTrue(plot_script.exists())
        source = plot_script.read_text(encoding="utf-8")
        self.assertNotIn("integrate_case(", source)
        self.assertNotIn("solve_bdf(", source)
        self.assertIn("q2_balance.npz", source)
        self.assertIn("q2_balance_time_refined.npz", source)
        self.assertIn("np.abs", source)
        self.assertIn("time_s <= 60.0", source)

    def test_degradation_scope_and_stable_reproduction_recipe(self):
        checks = supplementary._annotate_q1_checks(
            {"method": "BDF with analytic sparse Jacobian, same spatial discretization"}
        )
        self.assertIn("internal", checks["method_scope"])
        self.assertIn("cross-program", checks["method_scope"])
        self.assertIn("thermal surface", checks["cross_program_discretization_scope"])
        commands = supplementary._stable_reproduction_commands()
        self.assertEqual(len(commands), 4)
        self.assertIn("q2_supplementary_validation.py degradation", commands[0])
        self.assertIn("q2_supplementary_validation.py balance", commands[1])
        self.assertIn("q2_supplementary_validation.py all", commands[2])
        self.assertIn("plot_q2_supplementary_validation.py all", commands[3])
        validation_source = (
            CODE_DIR / "q2_supplementary_validation.py"
        ).read_text(encoding="utf-8")
        self.assertIn("plot_q2_supplementary_validation.py", validation_source)
        self.assertNotIn("q2_supplementary_plots.py", validation_source)

    def test_combined_summary_refresh_reads_both_completed_experiments(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            (output / "q2_degradation.json").write_text(
                json.dumps({"status": "PASS", "experiment": "degradation"}),
                encoding="utf-8",
            )
            (output / "q2_balance.json").write_text(
                json.dumps({"status": "PASS", "experiment": "balance"}),
                encoding="utf-8",
            )
            combined = supplementary.refresh_combined_summary(output)
            self.assertEqual(set(combined), {"degradation", "balance"})
            self.assertEqual(combined["degradation"]["status"], "PASS")
            self.assertEqual(combined["balance"]["status"], "PASS")
            stored = json.loads(
                (output / "q2_supplementary_validation.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(set(stored), {"degradation", "balance"})

    def test_degradation_finest_gate_and_one_level_expansion(self):
        coarse_fail = [{"n_intervals": 320, "comparison": {"pass": False}}]
        self.assertEqual(
            supplementary._append_next_degradation_grid([320, 640, 1280], 1280, False),
            [320, 640, 1280, 2560],
        )
        self.assertEqual(
            supplementary._append_next_degradation_grid([320, 640, 1280], 1280, True),
            [320, 640, 1280],
        )
        self.assertEqual(
            supplementary._append_next_degradation_grid([10240], 10240, False),
            [10240],
        )
        final_pass = [
            {"n_intervals": 320, "comparison": {"pass": False}},
            {"n_intervals": 1280, "comparison": {"pass": True}},
        ]
        self.assertFalse(supplementary._degradation_gate(coarse_fail, True, True))
        self.assertTrue(supplementary._degradation_gate(final_pass, True, True))

    @staticmethod
    def _selected_summary(integer_max, accepted_max, floor, passed=True):
        return {
            "quadrature_floor_max_abs_B": floor,
            "integer_seconds": {"max_abs_B": integer_max, "pass": passed},
            "accepted_nodes": {"max_abs_B": accepted_max, "pass": passed},
        }

    def test_balance_platform_exception_strict_drop_and_accepted_gate(self):
        baseline_selection = {"pass": True}
        refined_selection = {"pass": True}
        platform_base = self._selected_summary(2e-10, 3e-10, 1e-10)
        platform_refined = self._selected_summary(4e-10, 5e-10, 1e-10)
        platform = supplementary._balance_time_gate(platform_base, platform_refined)
        self.assertTrue(platform["platform_pass"])
        self.assertTrue(platform["pass"])

        no_drop_base = self._selected_summary(5e-6, 4e-6, 1e-7)
        no_drop_refined = self._selected_summary(5e-6, 4e-6, 1e-7)
        no_drop = supplementary._balance_time_gate(no_drop_base, no_drop_refined)
        self.assertFalse(no_drop["platform_pass"])
        self.assertFalse(no_drop["pass"])

        reduced_refined = self._selected_summary(3e-6, 2e-6, 1e-7)
        reduced = supplementary._balance_time_gate(no_drop_base, reduced_refined)
        self.assertTrue(reduced["strict_reduction_pass"])
        self.assertTrue(reduced["pass"])

        accepted_fail = self._selected_summary(1e-8, 2e-7, 1e-10, passed=False)
        self.assertFalse(
            supplementary._balance_gate(
                baseline_selection,
                refined_selection,
                accepted_fail,
                platform_refined,
                True,
                True,
                True,
            )
        )
        self.assertFalse(
            supplementary._balance_gate(
                baseline_selection,
                {"pass": False},
                platform_base,
                platform_refined,
                True,
                True,
                True,
            )
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
