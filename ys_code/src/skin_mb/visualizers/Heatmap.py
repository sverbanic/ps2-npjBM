def _log_zero_resolver(values, method='min'):
    """
    return a log transformed values with zero reolved
    Args:
        values:
        method:

    Returns:

    """
    import numpy as np
    import pandas as pd

    if isinstance(values, pd.DataFrame):
        return_np = False
    else:
        values = pd.DataFrame(data=values)
        return_np = True

    if method == 'min':
        replace_to = values.min().min()
    elif method == 'max':
        replace_to = values.max().max()
    elif method == '1':
        replace_to = 1
    elif isinstance(method, float):
        replace_to = method
    else:
        replace_to = np.nan

    values.replace(to_replace=0, value=replace_to, inplace=True)
    values = np.log10(values[values > 0])

    if return_np:
        return values.values
    else:
        return values


def heatmap(values, ax, origin='lower', bg_color='#FFFFFF',
            z_log=False, log_zero_resolve='min', nan_color='#151515',
            y_label_off=False, x_label_off=False,
            y_labels=None, x_labels=None,
            y_label_colors=None, x_label_colors=None,
            cmap=None, show_cbar=True, cbar_label=None, cbar_pos='top right',
            midpoint=None, **kwargs):

    import numpy as np
    import pandas as pd
    import matplotlib.colors as colors

    if isinstance(values, list):
        values = np.array(values, dtype=float)
    if isinstance(values, np.ndarray):
        values = pd.DataFrame(data=values, dtype=float)
    values = values.astype(dtype=float)

    if cbar_label is None:
        cbar_label = 'z'

    if z_log:
        if np.min(np.min(values)) < 0:
            kwargs['norm'] = colors.SymLogNorm(linthresh=1e-6, linscale=1e-6)
            log_ticks = np.linspace(-1, -5, 5)
            ticks = [-10 ** value for value in log_ticks] + [0] + [10 ** value for value in log_ticks]
        else:
            kwargs['norm'] = colors.LogNorm()
        # values = _log_zero_nan_resolver(values, method=log_zero_resolve)
        cbar_label = r'$\log_{10}$(' + cbar_label + ')'
    elif midpoint is not None:
        from .PlotTools import MidpointFixNorm
        kwargs['norm'] = MidpointFixNorm(vmin=values.min().min(),
                                         vmax=values.max().max(),
                                         midpoint=midpoint)
    if cmap is None:
        if midpoint is None:
            import matplotlib as mpl
            cmap = mpl.cm.get_cmap('Blues')
        else:
            from .PlotTools import CMapsDi
            cmap = CMapsDi.BluWhtRed
    if nan_color.lower() in ['bg', 'background']:
        cmap.set_bad(color=bg_color)
    else:
        cmap.set_bad(color=nan_color)
    ax.set_facecolor(bg_color)
    im = ax.imshow(values, cmap=cmap, origin=origin, aspect='auto', **kwargs)

    if show_cbar:
        from .PlotTools import add_color_bar
        if values.shape[0] < 10 and cbar_pos == 'top right':
            import matplotlib.pyplot as plt
            fig = plt.gca()
            ax_box = ax.get_position().bounds
            height = 1
            width = 0.05 * ax_box[2]
            x = ax_box[0] + ax_box[2] + 0.01
            y = ax_box[1] + ax_box[3] - height
            cbar_pos = (x, y, width, height)
            cax = fig.add_axes(cbar_pos)
            cbar = add_color_bar(target=im, ax=ax, cax=cax, cbar_label=cbar_label, orientation='vertical',
                                 ticks=ticks, **kwargs)
        else:
            cbar = add_color_bar(target=im, ax=ax, cbar_label=cbar_label, cbar_pos=cbar_pos, ticks=ticks, **kwargs)
        if 'norm' in kwargs:
            def remove_center_ticks(ticks):
                if len(ticks) % 2 == 0:
                    return ticks[:int((len(ticks) - 2)/2)] + ticks[int((len(ticks) + 2)/2):]
                else:
                    return ticks[:int((len(ticks) - 1)/2)] + ticks[int((len(ticks) + 1)/2):]

    if not y_label_off:
        ax.set_yticks(np.linspace(0, values.shape[0] - 1, values.shape[0]))
        if y_labels is None:
            y_labels = values.index
        ax.set_yticklabels(y_labels, fontdict={'fontsize': 12})
        if y_label_colors is not None:
            for ytick, color in zip(ax.get_yticklabels(), y_label_colors):
                ytick.set_color(color)
    if not x_label_off:
        ax.set_xticks(np.linspace(0, values.shape[1] - 1, values.shape[1]))
        if x_labels is None:
            x_labels = values.columns
        ax.set_xticklabels(x_labels, fontdict={'fontsize': 12, 'rotation': 90})
        if x_label_colors is not None:
            for xtick, color in zip(ax.get_xticklabels(), x_label_colors):
                xtick.set_color(color)
    return im