"""Check the concrete table and unified-figure contracts of the Q2 revision."""
from pathlib import Path
import csv
import json
import re
import numpy as np


def main():
    root = Path(__file__).resolve().parents[1]
    result = root / 'results/q2'
    paper = root / 'paper'
    times = [1800, 3600, 5400, 7200, 9000, 10800]
    columns = [0, 5, 10, 15, 20]
    checks = {}
    with np.load(result / 'q2_fields.npz') as original:
        for name, key in [('temperature', 'temperature_C'),
                          ('moisture', 'moisture_kg_per_kg')]:
            tex = (paper / f'tables/q2_{name}.tex').read_text(encoding='utf-8')
            rows = re.findall(r'^\s*(?:0\.5|[1-3]\.[05]) & ([0-9. &]+)', tex, re.M)
            assert len(rows) == 6, (name, 'six time rows required')
            with (result / f'q2_{name}_paper.csv').open(encoding='utf-8-sig') as stream:
                csv_rows = list(csv.reader(stream))[1:]
            for i, row in enumerate(rows):
                values = [v.strip() for v in row.split('&')]
                assert values == [f'{v:.4f}' for v in original[key][times[i], columns]]
                assert values == [f'{float(v):.4f}' for v in csv_rows[i][1:]]
            for token in [r'\multirow{2}{*}{时间/h}',
                          r'\multicolumn{5}{|c|}{到药材中心的距离/cm}',
                          r'\cline{2-6}', r'\fontsize{10.5pt}',
                          r'\begin{tabularx}{0.92\textwidth}']:
                assert token in tex, (name, token)
            assert tex.count(r'\hline') == 8
            checks[name + '_table'] = 'PASS: 30 values, original grid and merged headers'
        metadata = json.loads((result / 'q2_combined_fields_metadata.json').read_text(encoding='utf-8'))
        with np.load(result / 'q2_combined_fields_data.npz') as combined:
            indices = combined['surface_indices']
            assert np.array_equal(combined['time_s'], original['time_s'])
            assert np.array_equal(combined['radius_cm'], original['radius_cm'])
            assert np.array_equal(combined['temperature_C'], original['temperature_C'])
            assert np.array_equal(combined['moisture_kg_per_kg'], original['moisture_kg_per_kg'])
        assert np.array_equal(indices, np.unique(np.r_[np.arange(61), np.arange(60, 10801, 60)]))
        assert metadata['source'] == 'q2_fields.npz'
        assert metadata['smoothing'] is False and metadata['extrapolation'] is False
        checks['combined_figure_data'] = 'PASS: all 10801 x 21 states preserved; surface uses exact source samples'
    name = 'q2_combined_fields.pdf'
    assert (paper / 'figures' / name).read_bytes() == (result / name).read_bytes()
    checks['paper_figure_source'] = 'PASS: embedded combined PDF equals the result PDF'
    source = (paper / 'sections/6_q2_full_drying.tex').read_text(encoding='utf-8')
    assert r'\exp' not in source
    assert r'\farc' not in source
    assert 'fig:q2-combined-fields' in source
    assert 'fig:q2-spacetime' not in source
    assert 'fig:q2-fields' not in source
    assert 'q2_combined_fields.pdf' in source
    assert 'q2_fields.pdf' not in source
    assert 'q2_spacetime.pdf' not in source
    checks['exponent_notation_and_figures'] = 'PASS'
    report = {'scope': 'Q2 user revision; no PDE rerun or official workbook modification',
              'status': 'PASS', 'checks': checks}
    (result / 'q2_paper_revision_validation.json').write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
