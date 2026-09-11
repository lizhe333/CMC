"""Render Q2 space-time surfaces from accepted data, without a new solve."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager, colors, ticker


def main():
    root = Path(__file__).resolve().parents[1]
    output = root / 'results/q2'
    with np.load(output / 'q2_fields.npz', allow_pickle=False) as fields:
        time = fields['time_s'].copy()
        radius = fields['radius_cm'].copy()
        temperature = fields['temperature_C'].copy()
        moisture = fields['moisture_kg_per_kg'].copy()
    assert temperature.shape == moisture.shape == (10801, 21)
    assert np.isfinite(temperature).all() and np.isfinite(moisture).all()
    # Retain early one-second samples and every minute thereafter, as well
    # as all 21 original spatial positions. No fitted/smoothed data.
    indices = np.unique(np.r_[np.arange(61), np.arange(60, 10801, 60)])
    x, y = np.meshgrid(time[indices] / 3600, radius, indexing='ij')
    available = {f.name for f in font_manager.fontManager.ttflist}
    family = next(n for n in ('Microsoft YaHei', 'SimHei', 'SimSun') if n in available)
    plt.rcParams.update({'font.family': family, 'font.size': 10,
                         'axes.unicode_minus': False, 'pdf.fonttype': 42,
                         'savefig.dpi': 300, 'axes.linewidth': .6})
    # Sequential maps use Okabe-Ito blue and vermillion as the endpoints;
    # redundant height and labelled bars keep the values readable in print.
    heat = colors.LinearSegmentedColormap.from_list('heat', ['#fff4d5', '#F0E442', '#E69F00', '#D55E00'])
    water = colors.LinearSegmentedColormap.from_list('water', ['#edf6fb', '#56B4E9', '#0072B2'])
    fig = plt.figure(figsize=(10.4, 5.2), facecolor='white')
    specifications = [
        (temperature, heat, (28, 51), [28, 34, 40, 46, 50], '(a) 温度场', '温度/°C'),
        (moisture, water, (1, 2.6), [1, 1.4, 1.8, 2.2, 2.55], '(b) 含水率场', '干基含水率/(kg/kg)'),
    ]
    for i, (data, cmap, limits, ticks, title, barlabel) in enumerate(specifications):
        ax = fig.add_subplot(1, 2, i+1, projection='3d')
        z = data[indices]
        surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cmap,
                               vmin=limits[0], vmax=limits[1], linewidth=0,
                               antialiased=False, shade=False, rasterized=True)
        # Six exact paper-time profiles guide depth perception without a
        # dense wire mesh. These curves use the original nodal samples.
        for second in (1800, 3600, 5400, 7200, 9000, 10800):
            ax.plot(np.full(21, second/3600), radius, data[second],
                    color='#263c4a', alpha=.50, linewidth=.55)
        ax.set(xlim=(0, 3), ylim=(2, 0), zlim=limits,
               xlabel='时间/h', ylabel='径向位置/cm')
        ax.set_xticks([0, 1, 2, 3])
        ax.set_yticks([0, .5, 1, 1.5, 2])
        ax.set_zticks(ticks)
        ax.zaxis.set_major_formatter(ticker.FuncFormatter(lambda value, pos: f'{value:g}'))
        ax.tick_params(labelsize=9, pad=0)
        ax.xaxis.labelpad = 6
        ax.yaxis.labelpad = 7
        ax.set_title(title, fontsize=12, fontweight='bold', pad=5)
        ax.set_proj_type('ortho')
        ax.set_box_aspect((1.25, 1, .82))
        ax.view_init(elev=27, azim=-125)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.set_facecolor((.96, .97, .98, .45))
            axis.pane.set_edgecolor('#dbe1e5')
            axis._axinfo['grid'].update(color='#d9e0e5', linewidth=.45)
        bar = fig.colorbar(surf, ax=ax, orientation='horizontal',
                           shrink=.68, pad=.045, fraction=.035, aspect=25)
        bar.set_label(barlabel, labelpad=4, fontsize=10)
        bar.set_ticks([28, 34, 40, 46, 50] if i == 0 else [1, 1.4, 1.8, 2.2, 2.55])
        bar.ax.tick_params(labelsize=9, length=2)
        bar.outline.set_visible(False)
    fig.subplots_adjust(left=.015, right=.985, bottom=.12, top=.92, wspace=.06)
    fig.savefig(output / 'q2_spacetime.pdf', bbox_inches='tight', pad_inches=.08)
    fig.savefig(output / 'q2_spacetime.png', dpi=300, bbox_inches='tight', pad_inches=.08)
    plt.close(fig)
    np.savez_compressed(output / 'q2_spacetime_data.npz', time_s=time[indices],
                        radius_cm=radius, temperature_C=temperature[indices],
                        moisture_kg_per_kg=moisture[indices])
    metadata = {'source': 'q2_fields.npz', 'time_indices': indices.tolist(),
                'sampling': '1 s for first minute; 60 s afterwards; original 21 radial nodes',
                'smoothing': False, 'extrapolation': False,
                'temperature_range_C': [float(temperature.min()), float(temperature.max())],
                'moisture_range_kg_per_kg': [float(moisture.min()), float(moisture.max())],
                'shape': list(temperature[indices].shape),
                'paper_time_profiles_s': [1800, 3600, 5400, 7200, 9000, 10800]}
    (output / 'q2_spacetime_metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in metadata.items() if k != 'time_indices'}, ensure_ascii=False))


if __name__ == '__main__':
    main()
