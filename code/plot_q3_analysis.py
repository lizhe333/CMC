"""Paper-only Q3 visual analysis from accepted minute samples and event state.

No PDE solve, workbook write, smoothing, or extrapolation. Local crossing times
use piecewise-linear time interpolation; only the center time is the validated
solver event. Surface and contour plots depict the saved 21 radial positions.
"""
from pathlib import Path
import csv
import json
import shutil

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors, font_manager, ticker


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/q3'
THRESHOLD = 0.15


def load_data():
    validation = json.loads((OUT / 'q3_validation.json').read_text(encoding='utf-8'))
    if validation['status'] != 'NUMERICAL_CHECKS_PASSED':
        raise ValueError('Accepted Q3 results are required')
    with np.load(OUT / 'q3_fields.npz', allow_pickle=False) as data:
        t = data['time_s'].astype(float)
        r = data['radius_cm'].copy()
        c = data['moisture_kg_per_kg'].copy()
    with np.load(OUT / 'q3_event.npz', allow_pickle=False) as event:
        end = float(event['time_s'])
        er = event['radius_nodes_m'].copy() * 100
        ec = event['moisture_kg_per_kg'].copy()
    assert c.shape == (len(t), len(r)) and len(r) == 21
    assert np.all(np.diff(t) == 60) and t[0] == 0
    assert t[-1] < end < t[-1] + 60
    assert abs(end - validation['event']['time_s']) < 1e-7
    assert r[0] >= er[0] and r[-1] <= er[-1]
    # Display uses linear spatial reconstruction; endpoint is never rounded.
    cend = np.interp(r, er, ec)
    t = np.r_[t, end]
    c = np.vstack([c, cend])
    assert np.isfinite(c).all() and np.min(c) > 0
    assert np.max(np.diff(c, axis=1)) < 1e-8
    assert np.max(np.diff(c, axis=0)) < 1e-8
    assert abs(c[-1, 0] - THRESHOLD) < 1e-12
    return t, r, c, validation


def first_crossings(t, r, c):
    times, brackets, residuals = [], [], []
    for j in range(len(r)):
        crossing = np.flatnonzero(c[:, j] <= THRESHOLD)
        assert len(crossing) > 0 and crossing[0] > 0
        b = int(crossing[0])
        a = b - 1
        assert np.all(c[:b, j] > THRESHOLD)
        fraction = (c[a, j] - THRESHOLD) / (c[a, j] - c[b, j])
        assert 0 <= fraction <= 1
        time = t[a] + fraction * (t[b] - t[a])
        value = c[a, j] + fraction * (c[b, j] - c[a, j])
        times.append(time)
        brackets.append((t[a], t[b]))
        residuals.append(abs(value - THRESHOLD))
    times = np.asarray(times)
    assert np.max(np.diff(times)) <= 1e-7
    assert abs(times[0] - t[-1]) < 1e-7
    assert max(residuals) < 1e-12
    return times, np.asarray(brackets), residuals


def save(fig, name):
    # Keep text readable after reducing a two-panel figure to paper width.
    from matplotlib.text import Text
    fig.canvas.draw()
    for item in fig.findobj(match=Text):
        item.set_fontsize(item.get_fontsize() * 1.4)
    fig.savefig(OUT / f'{name}.pdf', bbox_inches='tight', pad_inches=.12)
    fig.savefig(OUT / f'{name}.png', dpi=300, bbox_inches='tight', pad_inches=.12)
    plt.close(fig)
    for ext in ['pdf', 'png']:
        shutil.copy2(OUT / f'{name}.{ext}', ROOT / f'paper/figures/{name}.{ext}')


def main():
    t, r, c, validation = load_data()
    crossing, brackets, residuals = first_crossings(t, r, c)
    hours = t / 3600
    available = {f.name for f in font_manager.fontManager.ttflist}
    font = next(n for n in ('Microsoft YaHei', 'SimHei', 'SimSun') if n in available)
    plt.rcParams.update({'font.family': font, 'font.size': 10,
                         'axes.unicode_minus': False, 'pdf.fonttype': 42,
                         'axes.linewidth': .7, 'legend.frameon': False,
                         'legend.fontsize': 9})
    blue = '#0072B2'
    red = '#B64536'
    water = colors.LinearSegmentedColormap.from_list(
        'q3_water', ['#edf6fb', '#56B4E9', '#0072B2', '#12334b'])

    # Keep every minute for the first hour, then every ten minutes for the
    # overview surface. Contours retain every saved minute and the endpoint.
    indices = np.unique(np.r_[np.arange(min(61, len(t))),
                              np.arange(0, len(t), 10), len(t)-1])
    xx, yy = np.meshgrid(hours[indices], r, indexing='ij')
    fig = plt.figure(figsize=(11.2, 4.8), facecolor='white')
    ax = fig.add_subplot(121, projection='3d')
    surface = ax.plot_surface(xx, yy, c[indices], rstride=1, cstride=1,
                              cmap=water, vmin=.05, vmax=2.55, linewidth=0,
                              shade=False, antialiased=False, rasterized=True)
    ax.set(xlabel='时间/h', ylabel='径向位置/cm',
           xlim=(0, hours[-1]), ylim=(2, 0), zlim=(0, 2.6))
    ax.set_xticks([0, 12, 24, 36, 57.47])
    ax.set_xticklabels(['0', '12', '24', '36', '57.47'])
    ax.set_yticks([0, .5, 1, 1.5, 2])
    ax.set_zticks([0, .5, 1, 1.5, 2, 2.5])
    ax.tick_params(labelsize=8, pad=0)
    ax.set_proj_type('ortho')
    ax.set_box_aspect((1.35, 1, .85))
    ax.view_init(elev=26, azim=-125)
    ax.set_title('(a) 全程含水率曲面', fontsize=12, fontweight='bold', pad=8)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((.96, .97, .98, .45))
        axis._axinfo['grid'].update(color='#d9e0e5', linewidth=.45)
    bar = fig.colorbar(surface, ax=ax, orientation='horizontal',
                       shrink=.75, pad=.08, fraction=.045, aspect=26)
    bar.set_label('干基含水率/(kg/kg)', fontsize=10)
    bar.set_ticks([.05, .5, 1, 1.5, 2, 2.55])
    bar.ax.tick_params(labelsize=8)
    bar.outline.set_visible(False)

    ax = fig.add_subplot(122)
    late = hours >= 12
    levels = np.linspace(.05, .46, 42)
    fill = ax.contourf(hours[late], r, c[late].T, levels=levels, cmap=water, zorder=-1)
    # Rasterize the filled layer only, avoiding PDF polygon seam artifacts;
    # the threshold line, labels and axes remain vector objects.
    ax.set_rasterization_zorder(0)
    ax.contour(hours[late], r, c[late].T, levels=[THRESHOLD], colors=[red],
               linewidths=1.8)
    # The center endpoint is supplied by the solver, not extrapolated contouring.
    ax.scatter([hours[-1]], [0], s=38, color=red, clip_on=False, zorder=5)
    ax.annotate('中心最后达标', xy=(hours[-1], 0), xytext=(39, .26),
                 fontsize=9, color=red,
                 arrowprops={'arrowstyle': '->', 'color': red, 'lw': .8})
    ax.text(21, .43, '未达标区', color='#173347', fontsize=11,
            bbox={'facecolor': 'white', 'alpha': .75, 'edgecolor': 'none'})
    ax.text(40, 1.74, '已达标区', color='#173347', fontsize=11)
    ax.set(xlabel='时间/h', ylabel='径向位置/cm', xlim=(12, hours[-1]), ylim=(0, 2))
    ax.set_xticks([12, 24, 36, 48, 57.47])
    ax.set_xticklabels(['12', '24', '36', '48', '57.47'])
    ax.set_yticks([0, .5, 1, 1.5, 2])
    ax.set_title('(b) 后期达标分界', fontsize=12, fontweight='bold', pad=12)
    bar = fig.colorbar(fill, ax=ax, orientation='horizontal',
                       shrink=.95, pad=.20, fraction=.045, aspect=26)
    bar.set_ticks([.05, .15, .25, .35, .45])
    bar.set_label('干基含水率/(kg/kg)', fontsize=10)
    bar.outline.set_visible(False)
    bar.solids.set_edgecolor('face')
    fig.subplots_adjust(left=.01, right=.98, bottom=.18, top=.90, wspace=.18)
    save(fig, 'q3_spacetime')

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.1), layout='constrained')
    profile_h = [12, 24, 36, 48, hours[-1]]
    palette = [blue, '#56B4E9', '#009E73', '#E69F00', '#8B5190']
    markers = ['o', 's', '^', 'D', 'v']
    for hour, color, marker in zip(profile_h, palette, markers):
        i = int(np.flatnonzero(np.isclose(hours, hour, atol=1e-10, rtol=0))[0])
        label = f'{hour:.0f} h' if hour != hours[-1] else '终点 57.47 h'
        axes[0].plot(r, c[i], color=color, lw=1.5, marker=marker,
                      markevery=5, markersize=3.8, label=label)
    axes[0].axhline(THRESHOLD, color=red, ls='--', lw=1.2)
    axes[0].text(.06, .105, '达标阈值 0.15', color=red, fontsize=9)
    axes[0].set(xlabel='径向位置/cm', ylabel='干基含水率/(kg/kg)',
                xlim=(0, 2), ylim=(.04, .57))
    axes[0].set_title('(a) 代表时刻的径向分布', fontweight='bold', pad=12)
    axes[0].legend(loc='upper center', ncol=3, fontsize=8,
                   columnspacing=.7, handlelength=1.5, labelspacing=.25)
    axes[1].plot(r, crossing / 3600, '-o', color=blue, lw=1.6, markersize=3)
    for j in [0, 5, 10, 15, 20]:
        axes[1].scatter([r[j]], [crossing[j]/3600], color=red, s=23, zorder=4)
        axes[1].annotate(f'{crossing[j]/3600:.2f} h', (r[j], crossing[j]/3600),
                          xytext=((8 if j == 0 else -7 if j == 20 else 0), 9),
                          textcoords='offset points', fontsize=9,
                          ha=('left' if j == 0 else 'right' if j == 20 else 'center'))
    axes[1].set(xlabel='径向位置/cm', ylabel='首次达到阈值的时间/h',
                xlim=(-.04, 2.04), ylim=(9, 64))
    axes[1].set_yticks([12, 24, 36, 48, 60])
    axes[1].set_title('(b) 各位置的临界时间', fontweight='bold', pad=12)
    for ax in axes:
        ax.set_xticks([0, .5, 1, 1.5, 2])
        ax.grid(alpha=.2, linewidth=.5)
        ax.spines[['top', 'right']].set_visible(False)
    save(fig, 'q3_radial_completion')

    with (OUT / 'q3_local_crossing_times.csv').open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(['radius_cm', 'time_s', 'time_h', 'bracket_left_s', 'bracket_right_s'])
        for radius, time, bracket in zip(r, crossing, brackets):
            writer.writerow([radius, time, time/3600, *bracket])
    np.savez_compressed(OUT / 'q3_analysis_figure_data.npz', time_s=t, radius_cm=r,
                        moisture_kg_per_kg=c, surface_time_indices=indices,
                        crossing_times_s=crossing, crossing_brackets_s=brackets)
    metadata = {
        'sources': ['q3_fields.npz', 'q3_event.npz', 'q3_validation.json'],
        'solver_rerun': False, 'smoothing': False, 'extrapolation': False,
        'threshold_kg_per_kg': THRESHOLD, 'late_window_start_h': 12,
        'display_radial_points': len(r), 'display_time_points_with_endpoint': len(t),
        'surface_sampling': 'Every minute in first hour, every ten minutes later, exact endpoint appended',
        'time_interpolation': 'Piecewise linear between saved 60 s samples; exact center solver event',
        'spatial_reconstruction': 'Piecewise linear on 21 saved radial locations; no subgrid precision claim',
        'max_crossing_residual': max(residuals),
        'max_crossing_bracket_width_s': float(np.max(np.diff(brackets, axis=1))),
        'center_event_time_h': validation['event']['time_h'],
        'surface_crossing_time_h': float(crossing[-1]/3600),
        'center_surface_delay_h': float((crossing[0]-crossing[-1])/3600),
        'representative_crossings_h': {str(r[j]):float(crossing[j]/3600) for j in [0,5,10,15,20]},
        'checks': {'finite_positive': True, 'radial_monotonicity': True,
                   'temporal_monotonicity': True, 'first_crossing_brackets': True,
                   'center_matches_solver_event': True},
    }
    (OUT / 'q3_analysis_figures.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
