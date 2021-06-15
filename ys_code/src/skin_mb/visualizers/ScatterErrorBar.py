

def single_series(single_data, ax, entry_list=None, pos_mapper=None, shift=0,
                  marker=None, marker_color=None, marker_size=None, vertical=True,
                  sig_flag=None, **kwargs):
    import numpy as np
    import pandas as pd

    if not isinstance(single_data, pd.DataFrame):
        data = single_data.get_data(entry_list).sort_values('mid', ascending=False) #assume a handle to fetch data
    else:
        data = single_data.copy()

    if pos_mapper is None:
        data['pos'] = np.linspace(0, data.shape[0] - 1, data.shape[0], dtype=int)
    else:
        data['pos'] = data.index.map(pos_mapper)

    if marker is None:
        data['marker'] = np.repeat('o', data.shape[0])
    elif isinstance(marker, str):
        data['marker'] = np.repeat(marker, data.shape[0])
    elif isinstance(marker, dict):
        data['marker'] = data.index.map(marker)
    else:
        raise ValueError('Customized marker should be dict: {otu_id -> marker}')

    if marker_color is None:
        marker_color = ['#B2112A', '#AEAEAE', '#2C73B4'] # Default show significance, tri-color [pos, neu, neg]

    if isinstance(marker_color, list):
        if len(marker_color) == 3:
            def colorer(row):
                if sig_flag is not None:
                    if np.isnan(row['mid']):
                        return marker_color[1]
                    else:
                        return marker_color[1-row[sig_flag]] # first use flag col: -1 --> 2, 0 --> 1, 1 --> 0
                else:
                    # then use mid, upper, lower
                    if row['lower'] > 0:
                        return marker_color[0]
                    elif row['upper'] > 0:
                        return marker_color[1]
                    else:
                        return marker_color[2]
            data['color'] = data.apply(colorer, axis=1)
        elif len(marker_color) == data.shape[0]:
            data['color'] = marker_color
        else:
            raise ValueError('List of color needs to be 3 or the size of data')
    elif isinstance(marker_color, dict):
        data['color'] = data.index.map(marker_color)
    elif isinstance(marker_color, str):
        data['color'] = np.repeat(marker_color, data.shape[0])

    if marker_size is None:
        data['size'] = np.repeat(3, data.shape[0])
    elif isinstance(marker_size, float):
        data['size'] = np.repeat(marker_size, data.shape[0])
    elif isinstance(marker_size, dict):
        data['size'] = data.index.map(marker_size)
    else:
        raise ValueError('Customized size should be dict: {otu_id -> size}')

    if 'ls' not in data.columns:
        data['ls'] = np.repeat('-', data.shape[0])

    def plot_single_marker(row):
        if np.isnan(row['mid']):
            # If there is no record
            if vertical:
                ax.scatter(x=[0], y=[row['pos'] + shift], marker='o', s=row['size'], color=row['color'], zorder=3)
            else:
                ax.scatter(x=[row['pos'] + shift], y=[0], marker='o', s=row['size'], color=row['color'], zorder=3)
        else:
            if vertical:
                eb = ax.errorbar(x=[row['mid']], y=[row['pos'] + shift],
                                 xerr=[[row['mid'] - row['lower']], [row['upper'] - row['mid']]],
                                 marker=row['marker'], capsize=2, markersize=row['size'], lw=1,
                                 color=row['color'], zorder=2, **kwargs)
            else:
                eb = ax.errorbar(y=[row['mid']], x=[row['pos'] + shift],
                                 yerr=[[row['mid'] - row['lower']], [row['upper'] - row['mid']]],
                                 marker=row['marker'], capsize=2, markersize=row['size'], lw=1,
                                 color=row['color'], zorder=2, **kwargs)
            eb[-1][0].set_linestyle(row['ls'])

    data.apply(plot_single_marker, axis=1)
    return data


def multi_series(series_to_plot, ax=None, plot_config=None, entry_list=None, plot_scope=None, pos_mapper=None,
                 figsize=None, vertical=True, label_mapper=None, value_label=None, marker_colors=None, entry_label_off=False,
                 legend_config=None, **kwargs):
    import pandas as pd
    import numpy as np

    if isinstance(series_to_plot, list):
        if plot_config is None:
            series_to_plot = {'Series {}'.format(ix): {'data': data, 'plot_config': {}}
                              for ix, data in enumerate(series_to_plot)}
        elif isinstance(plot_config, list) and len(plot_config) == len(series_to_plot):
            series_to_plot = {'Series {}'.format(ix): {'data': data, 'plot_config': plot_config[ix]}
                              for ix, data in enumerate(series_to_plot)}
        else:
            raise ValueError('plot_config does not match series_to_plot')
    elif isinstance(series_to_plot, dict):
        if plot_config is None:
            plot_config = {}
        series_to_plot = {ix:{'data': data, 'plot_config': plot_config[ix] if ix in plot_config.keys() else {}}
                          for ix,data in series_to_plot.items()}

    if entry_list is None:
        # By default search for all data that are detected significant in any selected dataset
        entry_list = set()
        for series in series_to_plot.values():
            if isinstance(series['data'], pd.DataFrame):
                entry_list.update(series['data'].index.values)
            else:
                # assume a single data handler if not pd.DataFrame
                entry_list.update(list(series['data'].get_data().index.values))

        def nested_sort_fn(entry):
            def get_value(entry, series):
                if isinstance(series['data'], pd.DataFrame):
                    try:
                        return series['data'].loc[entry]['mid']
                    except KeyError:
                        return np.nan
                else:
                    try:
                        return series['data'].get_data(entry).loc[entry]['mid']
                    except KeyError:
                        return np.nan

            values = np.array([get_value(entry, series) for series in series_to_plot.values()])
            values = [np.sum(~np.isnan(values))] + list(values)
            return values

        # By default sort by point estimation in each data series
        entry_list = sorted(list(entry_list), key=nested_sort_fn, reverse=False)

    if plot_scope is not None:
        entry_list = [entry for entry in entry_list if entry in plot_scope]

    if pos_mapper is None:
        pos_mapper = {entry: ix for ix, entry in enumerate(entry_list)}

    if ax is None:
        import matplotlib.pyplot as plt
        if figsize is None:
            figsize = (6, len(entry_list) / 2) if vertical else (len(entry_list) / 2, 6)
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)

    plot_log = {}
    for ix, series in enumerate(series_to_plot.values()):
        if len(series_to_plot) > 1:
            shift = -0.2 + ix * (0.4 / (len(series_to_plot) - 1))
        else:
            shift = 0
        plot_log[ix] = single_series(series['data'], ax, entry_list=entry_list, pos_mapper=pos_mapper,
                                     shift=shift, vertical=vertical, sig_flag='sig', marker_color=marker_colors,
                                     **series['plot_config'], **kwargs)

    if label_mapper is not None:
        labels = [label_mapper[otu] for otu in pos_mapper.keys()]
    else:
        labels = [otu for otu in pos_mapper.keys()]

    if vertical:
        if not entry_label_off:
            ax.set_yticks(list(pos_mapper.values()))
            ax.set_yticklabels(position=list(pos_mapper.values()), labels=labels, fontsize=12)
            ax.tick_params(axis='x', labelsize=12)
        ylim = ax.get_ylim()
        ax.plot([0, 0], ylim, ls='--', color='#AEAEAE', zorder=1)
        ax.set_ylim(ylim)
        xlim = ax.get_xlim()
        xlim = [abs(xlim[0]), abs(xlim[1])]
        xlim = [-max(xlim), max(xlim)]
        for pos in pos_mapper.values():
            ax.plot(xlim, [pos - 0.5, pos - 0.5], ls='-', lw=1, color='#F8F7ED')
        ax.set_xlim(xlim)
        ax.grid(False, axis='y')
        ax.set_ylim(-0.5, len(pos_mapper) - 0.5)
        ax.set_xlabel(value_label, fontsize=12)
    else:
        if not entry_label_off:
            ax.set_xticks(list(pos_mapper.values()))
            ax.set_xticklabels(position=list(pos_mapper.values()), labels=labels, fontsize=12)
            ax.tick_params(axis='y', labelsize=12)
        xlim = ax.get_xlim()
        ax.plot(xlim, [0, 0], ls='--', color='#AEAEAE', zorder=1)
        ax.set_xlim(xlim)
        ylim = ax.get_ylim()
        ylim = [abs(ylim[0]), abs(ylim[1])]
        ylim = [-max(ylim), max(ylim)]
        for pos in pos_mapper.values():
            ax.plot([pos - 0.5, pos - 0.5], ylim, ls='-', lw=1, color='#F8F7ED')
        ax.set_ylim(ylim)
        ax.grid(False, axis='x')
        ax.set_ylim(-0.5, len(pos_mapper) - 0.5)
        ax.set_ylabel(value_label, fontsize=12)

    if legend_config is not None:
        from ..visualizers.PlotTools import make_legends
        handles, labels = make_legends(legend_config)
        ax.legend(handles=handles, labels=labels, bbox_to_anchor=(0.1, 1, 0.8, 0.2), loc='lower center',
                  ncol=min(len(handles), 4), frameon=False)

    return plot_log