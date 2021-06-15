def three_panel_barplot(data_to_plot, percent=True, ylabel=None, ncol=3, bbox_to_anchor=None, **kwargs):
    import matplotlib.pyplot as plt
    from .PlotTools import add_title_bar

    fig, axes = plt.subplots(1, 3, figsize=[12, 6], sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0)
    sample_type_map = {
        'A': 'Pre-debridement',
        'B': 'Post-debridement',
        'C': 'Skin'
    }
    for ix, type in enumerate(['A', 'B', 'C']):
        sample_names = [sample for sample in data_to_plot.index if type in sample]
        data_to_plot.loc[sample_names].plot.bar(ax=axes[ix], legend=[], **kwargs)
        axes[ix].set_xticks([ix for ix, _ in enumerate(sample_names)])
        axes[ix].set_xticklabels([name[:-1] for name in sample_names], fontsize=12)
        axes[ix].set_xlabel('')
        axes[ix].grid(False)

    if percent:
        axes[0].set_ylim([0, 100])
        axes[0].tick_params(axis='y', direction='in', labelsize=12)
        axes[0].set_ylabel('Relative abundance (%)', fontsize=14)
    else:
        if ylabel is None:
            ylabel = 'Number of reads'
            axes[0].tick_params(axis='y', direction='in', labelsize=12)
            axes[0].set_ylabel(ylabel, fontsize=14)

    if bbox_to_anchor is None:
        bbox_to_anchor = (0, -0.18, 1, 0.2)
    axes[1].legend(loc='lower center', bbox_to_anchor=bbox_to_anchor, ncol=ncol, fontsize=14, frameon=False)
    for ix, type in enumerate(['A', 'B', 'C']):
        add_title_bar(ax=axes[ix], title=sample_type_map[type])

