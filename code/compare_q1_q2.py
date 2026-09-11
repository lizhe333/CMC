"""Compare accepted Q1/Q2 fields on their common domain; never rerun solvers."""
from pathlib import Path
import json
import numpy as np


def main():
    root = Path(__file__).resolve().parents[1]
    q1 = np.load(root / 'results/q1/q1_fields.npz', allow_pickle=False)
    q2 = np.load(root / 'results/q2/q2_fields.npz', allow_pickle=False)
    count = len(q1['time_s'])
    np.testing.assert_array_equal(q1['time_s'], q2['time_s'][:count])
    np.testing.assert_allclose(q1['radius_cm'], q2['radius_cm'], atol=1e-12, rtol=0)
    for name in ('temperature_C', 'moisture_kg_per_kg'):
        np.testing.assert_allclose(q1[name][0], q2[name][0], atol=1e-12, rtol=0)
    b1 = np.loadtxt(root / 'results/q1/boundary_used.csv', delimiter=',', skiprows=1, encoding='utf-8-sig')
    b2 = np.loadtxt(root / 'results/q2/boundary_used.csv', delimiter=',', skiprows=1, encoding='utf-8-sig')
    np.testing.assert_allclose(b1, b2[:len(b1)], atol=1e-12, rtol=0)
    c0, t0 = 2.55, 28.0
    rho, cp, k = 650+128*c0, 1450+2736*c0/(1+c0), .21+.38*c0/(1+c0)
    report = {
        'source_q1': 'results/q1/q1_fields.npz',
        'source_q2': 'results/q2/q2_fields.npz',
        'common_time_s': [0, int(q1['time_s'][-1])],
        'initial_and_boundary_match': True,
        'common_radius_m': .02, 'common_h': 25.0, 'common_hm': 8e-7,
        'initial_properties': {
            'q1': {'rho': 820.0, 'cp': 2600.0, 'k': .36,
                   'volumetric_heat_capacity': 820*2600, 'alpha': .36/(820*2600),
                   'D': 7e-9*np.exp(-.89/c0)},
            'q2': {'rho': rho, 'cp': cp, 'k': k,
                   'volumetric_heat_capacity': rho*cp, 'alpha': k/(rho*cp),
                   'D': 2.4e-3*np.exp(-.45/c0)*np.exp(-3850/(t0+273.15))},
        },
        'at_1800s': {}, 'maximum_difference_common_domain': {},
        'interpretation': 'Both constitutive relations and coupling differ. Differences are model comparisons, not experimental errors or isolated causal effects of coupling.',
    }
    for name in ('temperature_C', 'moisture_kg_per_kg'):
        delta = q2[name][:count] - q1[name]
        loc = np.unravel_index(np.argmax(np.abs(delta)), delta.shape)
        report['maximum_difference_common_domain'][name] = {
            'absolute_difference': float(abs(delta[loc])),
            'time_s': int(q1['time_s'][loc[0]]), 'radius_cm': float(q1['radius_cm'][loc[1]])}
        report['at_1800s'][name] = [
            {'radius_cm': float(q1['radius_cm'][j]), 'q1': float(q1[name][-1,j]),
             'q2': float(q2[name][count-1,j]), 'q2_minus_q1': float(delta[-1,j])}
            for j in (0,5,10,15,20)]
    restart = np.load(root / 'results/q2/q2_restart.npz', allow_pickle=False)
    report['q2_published_numerics'] = {key: float(restart[key]) for key in ('rtol','atol','max_step_s','n_intervals')}
    target = root / 'results/q2/q1_q2_comparison.json'
    target.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
