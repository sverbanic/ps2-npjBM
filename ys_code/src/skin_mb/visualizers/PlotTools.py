import matplotlib as mpl


class _CMapsDi:

    def __init__(self):
        pass

    @property
    def BluWhtOrg(self, reverse=False):
        return mpl.colors.LinearSegmentedColormap.from_list(name='BluWhtOrg', colors=['#F39730', '#FFFFFF', '#2C73B4'])

    @property
    def BluBlkOrg(self):
        return mpl.colors.LinearSegmentedColormap.from_list(name='BluBlkOrg', colors=['#F39730', '#000000', '#2C73B4'])

    @property
    def BluWhtYlw(self):
        return mpl.colors.LinearSegmentedColormap.from_list(name='BluWhtYlw', colors=['#F8DB36', '#FFFFFF', '#2C73B4'])

    @property
    def BluBlkYlw(self):
        return mpl.colors.LinearSegmentedColormap.from_list(name='BluBlkYlw', colors=['#F8DB36', '#000000', '#2C73B4'])

    def BluWhtRed(self, reverse=False):
        if reverse:
            return mpl.colors.LinearSegmentedColormap.from_list(name='BluWhtRed',
                                                                colors=['#2C73B4', '#FFFFFF', '#B2112A'])
        else:
            return mpl.colors.LinearSegmentedColormap.from_list(name='BluWhtRed',
                                                                colors=['#B2112A', '#FFFFFF', '#2C73B4'])

    @property
    def BluBlkRed(self):
        return mpl.colors.LinearSegmentedColormap.from_list(name='BluBlkRed', colors=['#B2112A', '#000000', '#2C73B4'])


CMapsDi = _CMapsDi()


class MidpointFixNorm(mpl.colors.Normalize):
    def __init__(self, vmin, vmax, midpoint=0, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        import scipy as sp
        normalized_min = max(0, 1 / 2 * (1 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
        normalized_max = min(1, 1 / 2 * (1 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return sp.ma.masked_array(sp.interp(value, x, y))


def add_title_bar(ax, title, pos='top', background='#AEAEAE', height=None, fontsize=14, *kwargs):
    """
    Add a text box with fix grey background to match the aethetics in ggplot
    Args:
        ax:
        title:
        pos:
        height:
        fontsize:
        *kwargs:

    Returns:

    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    if height is None:
        fig = plt.gcf()
        fig_size = fig.get_size_inches()
        ax_box = ax.get_position().bounds
        ax_box_ratio = np.array([ax_box[2], ax_box[3]]) * fig_size
        ax_box_ratio = ax_box_ratio[0]/ax_box_ratio[1]
        height = fontsize * 0.006 * (ylim[1] - ylim[0]) * ax_box_ratio
    bx = mpl.patches.Rectangle(xy=(xlim[0], ylim[1]),
                               width=(xlim[1] - xlim[0]), height=height,
                               clip_on=False, facecolor=background, edgecolor='w', linewidth=3, alpha=0.5, zorder=0)
    ax.text(s=title, fontsize=fontsize,
            x=(xlim[0] + xlim[1]) / 2,
            y=ylim[1] + height/2,
            va='center', ha='center',
            clip_on=False)
    ax.add_patch(bx)


def make_legends(legend_configs):
    """
    Args:
        data_configs: name: {marker, color, lt}

    Returns: a list of errbar handles
    """
    import matplotlib.pyplot as plt
    ax = plt.gca()

    def get_fake_series(kwargs):
        plot_type = kwargs.pop('type', None)
        if plot_type == 'scatter':
            return ax.scatter(y=[], x=[], s=20, **kwargs)
        else:
            return ax.errorbar(y=[], x=[], xerr=[[], []], ls='', capsize=3, markersize=3, **kwargs)

    return [get_fake_series(kwargs) for name, kwargs in legend_configs.items()], legend_configs.keys()


def add_color_bar(target, fig=None, ax=None, cax=None, orientation='vertical', ticks=None, cbar_pos='right top', cbar_label=None, **kwargs):
    import matplotlib.pyplot as plt
    if fig is None:
        fig = plt.gcf()
    if ax is None:
        ax = plt.gca()
    if cax is None:
        # If no specified Axes to put the color bar
        if cbar_pos == 'right top':
            ax_box = ax.get_position().bounds
            height = max(0.2 * ax_box[3], 0.3)
            width = 0.05 * ax_box[2]
            x = ax_box[0] + ax_box[2] + 0.01
            y = ax_box[1] + ax_box[3] - height
            cbar_pos = (x, y, width, height)
            orientation = 'vertical'
        elif cbar_pos == 'top right':
            ax_box = ax.get_position().bounds
            height = 0.05 * ax_box[3]
            width = 0.2 * ax_box[2]
            x = ax_box[0] + ax_box[2] - width
            y = ax_box[1] + ax_box[3] + 0.03
            cbar_pos = (x, y, width, height)
            orientation = 'horizontal'
        cax = fig.add_axes(cbar_pos)   # This position is relative to figure (0, 1)

    cb = fig.colorbar(target, cax=cax, ticks=ticks, orientation=orientation)
    if cbar_label:
        if orientation == 'vertical':
            cb.ax.set_ylabel(cbar_label, fontsize=12)
        else:
            cb.ax.set_xlabel(cbar_label, fontsize=12)
            cb.ax.yaxis.tick_right()
            ticks = cb.ax.get_xticklabels()
            cb.ax.set_xticklabels(labels=ticks, rotation=90)

    return cb


