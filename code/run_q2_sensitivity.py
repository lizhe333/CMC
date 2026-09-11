"""Reproduce Q2 +/-20% h/hm sensitivity from the saved formal baseline."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

CODE_DIR = Path(__file__).resolve().parent
PROJECT = CODE_DIR.parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

import q2_solver as q2


def main() -> None:
    output = PROJECT / "results" / "q2"
    fields = np.load(output / "q2_fields.npz", allow_pickle=False)
    restart = np.load(output / "q2_restart.npz", allow_pickle=False)
    validation_path = output / "q2_validation.json"
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    spatial_records = validation.get("spatial_records", [])
    baseline_runtime = float(spatial_records[-1].get("runtime_s", 0.0)) if spatial_records else 0.0
    result = {
        "n_intervals": int(restart["n_intervals"]),
        "times_s": fields["time_s"],
        "radius_output_cm": fields["radius_cm"],
        "temperature_C": fields["temperature_C"],
        "moisture_kg_per_kg": fields["moisture_kg_per_kg"],
        "mean_moisture": fields["mean_moisture"],
        "cumulative_loss": fields["cumulative_loss"],
        "mass_balance_residual": fields["mass_balance_residual"],
        "final_time_s": float(restart["final_time_s"]),
        "radius_nodes_m": restart["radius_nodes_m"],
        "control_volume_weights_m2": restart["control_volume_weights_m2"],
        "temperature_final_C": restart["temperature_final_C"],
        "moisture_final_kg_per_kg": restart["moisture_final_kg_per_kg"],
        "runtime_s": baseline_runtime,
        "h_W_m2K": float(restart["h_W_m2K"]),
        "hm_m_s": float(restart["hm_m_s"]),
        "rtol": float(restart["rtol"]),
        "atol": float(restart["atol"]),
        "max_step_s": float(restart["max_step_s"]),
    }
    end = int(result["final_time_s"])
    boundary, _ = q2.read_boundary(q2.DEFAULT_BOUNDARY, end=end)
    sensitivity = q2.sensitivity_runs(boundary, result, output, end)
    validation["sensitivity"] = sensitivity
    validation["sensitivity_grid"] = {
        "n_intervals": int(result["n_intervals"]),
        "end_s": end,
        "one_factor_levels": "baseline and +/-20 percent",
    }
    validation_path.write_text(
        json.dumps(q2._json_safe(validation), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (output / "q2_parameters_environment.json").write_text(
        json.dumps(
            {
                "parameters": validation.get("parameters", {}),
                "environment": validation.get("environment", {}),
                "boundary_source": validation.get("boundary_source"),
                "boundary_sheet": validation.get("boundary_sheet"),
                "boundary_interpolation": validation.get("boundary_interpolation"),
                "model": validation.get("model"),
                "moisture_parameterization": validation.get("moisture_parameterization"),
                "requested_end_s": validation.get("requested_end_s"),
                "sensitivity_grid": validation["sensitivity_grid"],
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print(json.dumps(q2._json_safe(sensitivity), ensure_ascii=False))


if __name__ == "__main__":
    main()
