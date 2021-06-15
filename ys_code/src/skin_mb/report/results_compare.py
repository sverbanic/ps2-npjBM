

def tri_panel_plot(dataset, compare,
                   otu_list=None, otu_scope=None, add_otu=None, otu_sort_by='BGLMM',
                   otu_group=None, otu_group_bg_color=None,
                   da_percent=False, transform_pipe=None, z_log=False, cbar_label='diff. rel. abun.',
                   subject_to_show=None,
                   label_mapper='genus name', otu_marker_colors=None,
                   panel_order=('bayes', 'deseq', 'diff_abun'),
                   bglmm_label=None,
                   dry_run=False, **kwargs):
    from ..data import OtuTable, BayesianRes, DeseqResult, DiffAbunRes

    otu_table = OtuTable.load_default(dataset)
    if transform_pipe is None:
        transform_pipe = ['rel abun', 'diff']
    diff_abun_res = DiffAbunRes.load_default(dataset=dataset, transform_pipe=transform_pipe, percent=da_percent)
    bayes_res = BayesianRes.load_default(dataset)
    deseq_res = DeseqResult.load_default(dataset)

    if otu_list is None:
        # by default, get all the significant OTUs
        otu_list = set(bayes_res.get_data(compare).sig_otu.index) | set(deseq_res.get_data(compare).sig_otu.index)

        def sort_on(otu, method):
            if otu in method.get_data(compare).sig_otu.index:
                return method.get_data(compare).get_data(otu).loc[otu]['mid']
            else:
                return -1e10

        if otu_sort_by == 'BGLMM':
            otu_list = sorted(list(otu_list), key=lambda row: (sort_on(row, bayes_res), sort_on(row, deseq_res)))
        elif otu_sort_by == 'DESeq2':
            otu_list = sorted(list(otu_list), key=lambda row: (sort_on(row, deseq_res), sort_on(row, bayes_res)))
        elif otu_sort_by.lower() in ['abundance', 'abun', 'rel abun']:
            rel_abun = otu_table.otu_stats()['rel_abun_avg_on_all'].to_dict()
            otu_list = sorted(list(otu_list), key=lambda row: rel_abun[row], reverse=False)

    if otu_scope is not None:
        if otu_scope in ['abun', '>0.1%']:
            otu_scope = otu_table.otu_stats()[otu_table.otu_stats()['rel_abun_avg_on_all'] >= 0.001].index
        elif otu_scope in ['not abun', '<0.1%']:
            otu_scope = otu_table.otu_stats()[otu_table.otu_stats()['rel_abun_avg_on_all'] <= 0.001].index
        otu_list = [otu for otu in otu_list if otu in otu_scope]

    if add_otu is not None:
        for otu in add_otu:
            if otu not in otu_list:
                otu_list.append(otu)

    if otu_group_bg_color is None:
        otu_group_bg_color = {}
    if otu_group is not None:
        otu_groups = [{'name': group_key,
                       'otu_list': [otu for otu in otu_list if otu in group_member],
                       'bg_color': otu_group_bg_color.get(group_key, None)}
                      for group_key, group_member in otu_group.items()]
        classified = []
        for group in otu_groups:
            classified += group['otu_list']

        otu_groups.append({'name': 'Other',
                           'otu_list': [otu for otu in otu_list if otu not in classified],
                           'bg_color': otu_group_bg_color.get('Other', None)})
        otu_groups = otu_groups[::-1]
    else:
        otu_groups = [{'name': None, 'otu_list': otu_list, 'bg_color': None}]

    otu_list = []
    for group in otu_groups:
        otu_list += list(group['otu_list'])

    if dry_run:
        return otu_list

    if label_mapper in ['genus name', 'genus_name']:
        from ..data import table
        label_mapper = table.get_otu_id_to_taxo_name().to_dict()

    # Plot
    if 'figsize' not in kwargs.keys():
        kwargs['figsize'] = (12, len(otu_list) / 4)
    if 'gridspec_kw' not in kwargs.keys():
        kwargs['gridspec_kw'] = {'width_ratios': [1, 1, 1]}

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=kwargs['figsize'], gridspec_kw=kwargs['gridspec_kw'])
    fig.subplots_adjust(wspace=0.02)
    ax_map = {data: ix for ix,data in enumerate(panel_order)}

    # Plot heatmap
    if ax_map['diff_abun'] == 0:
        otu_label_off = False
    else:
        otu_label_off = True
        axes[ax_map['diff_abun']].set_yticks([])
    ax = axes[ax_map['diff_abun']]
    if subject_to_show is None:
        subject_to_show = sorted(list(set([int(sample[:-1]) for sample in bayes_res.sample_modeled])))
    diff_abun_res.heatmap(compare=compare, ax=ax, otu_list=otu_list, subject_list=subject_to_show, adjust_patient_id=True,
                          y_label_off=otu_label_off, cbar_pos='right top' if ax_map['diff_abun'] == 2 else 'top right',
                          label_mapper=label_mapper, midpoint=0, cbar_label=cbar_label, z_log=z_log, **kwargs)

    # Plot beta on axes[1]
    if ax_map['bayes'] == 0:
        otu_label_off = False
    else:
        otu_label_off = True
        axes[ax_map['bayes']].set_yticks([])
    ax = axes[ax_map['bayes']]

    if bglmm_label is None:
        if compare in ['pre_vs_skin']:
            bglmm_label = r'BGLMM ($\beta_{j1}$)'
        elif compare in ['pre_vs_post']:
            bglmm_label = r'BGLMM ($\beta_{j1} - \beta_{j2}$)'
        else:
            bglmm_label = r'BGLMM ($\beta$)'

    bayes_res.scatter_plot_with_bar(data_to_plot=compare, ax=ax, otu_list=otu_list, marker_colors=otu_marker_colors,
                                    otu_label_off=otu_label_off, label_mapper=label_mapper, value_label=bglmm_label)
    if len(otu_groups) > 1:
        _group_bg_colorer(ax, otu_groups, show_label= otu_label_off)

    # Plot beta on axes[2]
    if ax_map['deseq'] == 0:
        otu_label_off = False
    else:
        otu_label_off = True
        axes[ax_map['deseq']].set_yticks([])
    ax = axes[ax_map['deseq']]
    deseq_res.scatter_plot_with_bar(data_to_plot=compare, ax=ax, otu_list=otu_list, marker_colors=otu_marker_colors,
                                    otu_label_off=otu_label_off, label_mapper=label_mapper)
    if len(otu_groups) > 1:
        _group_bg_colorer(ax, otu_groups, show_label= otu_label_off, vertical=False)

    from ..visualizers.PlotTools import make_legends

    if compare == 'pre_vs_skin':
        handle_name = ['Enriched in wound', 'Enriched in skin']
    elif compare == 'pre_vs_post':
        handle_name = ['Enriched in pre-debridement', 'Enriched in post-debridement']
    else:
        raise ValueError("Compare has to be 'pre_vs_skin' or 'pre_vs_post'")

    handles, labels = make_legends(legend_configs={handle_name[0]: {'marker': '', 'color': '#B2112A'},
                                                   handle_name[1]: {'marker': '', 'color': '#2C73B4'},
                                                   'Not significant': {'marker': '', 'color': '#AEAEAE'}})
    ax.legend(handles=handles, labels=labels,
              bbox_to_anchor=(0.2, 0.875, 0.6, 0.5), loc='lower center',
              ncol=min(len(handles), 4), frameon=False, fontsize=12, bbox_transform=fig.transFigure)
    fig.align_xlabels(axes)
    plt.show()
    return otu_list


def _group_bg_colorer(ax, otu_groups, show_label=True, vertical=True):
    import numpy as np
    group_sizes = [len(group['otu_list']) for group in otu_groups]
    group_bounds = [(np.sum(group_sizes[:group_ix]) - 0.5, np.sum(group_sizes[:group_ix + 1]) - 0.5)
                    for group_ix, group_size in enumerate(group_sizes)]
    xlim = ax.get_xlim()
    for ix, group in enumerate(otu_groups):
        ax.plot(xlim, [group_bounds[ix][1], group_bounds[ix][1]], '#FFFFFF', lw=1.5)
        if group['bg_color'] is not None:
            ax.axhspan(group_bounds[ix][0], group_bounds[ix][1], facecolor=group['bg_color'], alpha=0.3, zorder=0)
        if show_label:
            if vertical:
                ax.text(s=group['name'], x=xlim[1] + 0.5, y=np.mean(group_bounds[ix]),
                        fontsize=10, ha='left', va='center', rotation=90, alpha=0.6)
            else:
                ax.text(s=group['name'] + ' ', x=xlim[1], y= group_bounds[ix][1] - 0.5,
                        fontsize=8, ha='right', va='center', alpha=0.6)
    ax.set_xlim(xlim)


def results_robustness(model='BGLMM', compare='pre_vs_skin', show_only_abun=False, use_phylo_name=True,
                       figsize=None):
    from ..data.table import OtuTable

    if show_only_abun:
        otu_table = OtuTable.load_default('filtered')
        otu_scope = otu_table.otu_stats()[otu_table.otu_stats()['rel_abun_avg_on_all'] >= 0.001].index
    else:
        otu_scope = None

    if use_phylo_name:
        from ..data.table import get_otu_id_to_taxo_name
        label_mapper = get_otu_id_to_taxo_name().to_dict()
    else:
        label_mapper = None

    if model == 'BGLMM':
        from ..data.bayesian import BayesianRes
        model_handle = BayesianRes
        if compare in ['pre_vs_skin']:
            value_label = r'BGLMM ($\beta_{j1}$)'
        elif compare in ['pre_vs_post']:
            value_label = r'BGLMM ($\beta_{j1} - \beta_{j2}$)'
        else:
            value_label = r'BGLMM ($\beta$)'
    elif model == 'DESeq2':
        from ..data.deseq2 import DeseqResult
        model_handle = DeseqResult
        value_label = r'DESeq2 ($\log_2$ fold change)'
    else:
        raise ValueError("model should be 'BGLMM' or 'DESeq2'")

    # Plot
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(5, 10) if figsize is None else figsize)
    ax = fig.add_subplot(111)

    from ..visualizers.ScatterErrorBar import multi_series
    series_to_plot = {'With patient 16': model_handle.load_default('filtered').get_data(compare),
                      'Without patient 16': model_handle.load_default('filtered_no_p16').get_data(compare)}

    return multi_series(series_to_plot=series_to_plot, ax=ax,
                        plot_config={'With patient 16': {'marker': 'o'},
                                     'Without patient 16': {'marker': '^'}},
                        entry_list=None, plot_scope=otu_scope, pos_mapper=None,
                        vertical=True, label_mapper=label_mapper,
                        value_label=value_label, entry_label_off=False,
                        legend_config={'With patient 16': {'marker': 'o', 'color':'k'},
                                       'Without patient 16': {'marker': '^', 'color': 'k'}})


def sig_otu_beta_deseq2_heal_unheal(dataset='filtered_no_p16', heal_vs_unheal=True,
                                    otu_list=None, otu_scope=None, label_mapper=None,
                                    otu_marker_colors=None, metabolic_groups=None, group_bg_colors=None,
                                    panel_order=('bayes', 'deseq'), **kwargs):
    from ..data import BayesianRes, DeseqResult

    bayes_res = BayesianRes.load_default(dataset)
    deseq_res = DeseqResult.load_default(dataset)

    # Plot
    if 'figsize' not in kwargs.keys():
        kwargs['figsize'] = (10, 10)
    if 'gridspec_kw' not in kwargs.keys():
        kwargs['gridspec_kw'] = {'width_ratios': [1, 1]}

    if heal_vs_unheal:
        compares = ['healed_vs_unhealed_pre', 'healed_vs_unhealed_post']
        handle_name = ['Healed', 'Unhealed', 'Pre-debridement', 'Post-debridement']
    else:
        compares = ['pre_vs_post_healed', 'pre_vs_post_unhealed']
        handle_name = ['Enriched in Pre-debridement', 'Enriched in Post-debridement', 'Healed', 'Unhealed']

    if otu_list is None:
        otu_list = set()
        for compare in compares:
            otu_list.update(set(bayes_res.get_data(compare).sig_otu.index) | set(deseq_res.get_data(compare).sig_otu.index))

    if otu_scope is not None:
        if otu_scope in ['abun']:
            from ..data import OtuTable

            otu_table = OtuTable.load_default(dataset)
            otu_scope = otu_table.otu_stats()[otu_table.otu_stats()['rel_abun_avg_on_all'] >= 0.001].index
        otu_list = [otu for otu in otu_list if otu in otu_scope]

    if label_mapper in ['genus_name']:
        from ..data.table import get_otu_id_to_taxo_name
        label_mapper = get_otu_id_to_taxo_name().to_dict()

    if metabolic_groups is not None:
        if label_mapper is None:
            raise ValueError('label_mapper should be given if use metabolic_groups')

        if group_bg_colors is None:
            group_bg_colors = {}
        otu_groups = [{'name': group_key,
                       'otu_list': [otu for otu in otu_list if label_mapper[otu] in group_items],
                       'bg_color': group_bg_colors.get(group_key, None)}
                      for group_key, group_items in metabolic_groups.items()]
        classified = []
        for group in otu_groups:
            classified += group['otu_list']
        otu_groups.append({'name': 'Other',
                           'otu_list':[otu for otu in otu_list if otu not in classified],
                           'bg_color': group_bg_colors.get('Other', None)})
        otu_groups = otu_groups[::-1]
    else:
        otu_groups = [{'name': None, 'otu_list': otu_list, 'bg_color': None}]

    otu_list = []
    for group in otu_groups:

        def sort_fn(otu):
            import numpy as np
            return np.sum([np.sum([otu in res.get_data(compare).sig_otu.index
                                   for res in [bayes_res, deseq_res]])
                           for compare in compares])

        otu_list += sorted(list(group['otu_list']), key=sort_fn)

    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, **kwargs)
    fig.subplots_adjust(wspace=0.02)

    ax_map = {data: ix for ix,data in enumerate(panel_order)}
    if otu_marker_colors is None:
        otu_marker_colors = ['#B2112A', '#AEAEAE', '#2C73B4']

    # Plot beta on axes[1]
    if ax_map['bayes'] == 0:
        otu_label_off = False
    else:
        otu_label_off = True
        axes[ax_map['bayes']].set_yticks([])
        axes[ax_map['bayes']].set_yticks([])
    ax = axes[ax_map['bayes']]

    bayes_res.scatter_plot_with_bar(
        data_to_plot=compares,
        plot_config={
            compares[0]:{'marker': 'o'},
            compares[1]:{'marker': '^'}
        },
        ax=ax, otu_list=otu_list, marker_colors=otu_marker_colors, otu_label_off=otu_label_off, label_mapper=label_mapper)
    ax.set_xlabel(r'BGLMM ($\beta$)', fontsize=14)

    # Plot beta on axes[2]
    if ax_map['deseq'] == 0:
        otu_label_off = False
    else:
        otu_label_off = True
        axes[ax_map['deseq']].set_yticks([])
        axes[ax_map['deseq']].set_yticks([])

    ax = axes[ax_map['deseq']]
    deseq_res.scatter_plot_with_bar(
        data_to_plot=compares,
        plot_config={
            compares[0]: {'marker': 'o'},
            compares[1]: {'marker': '^'}
        },
        ax=ax, otu_list=otu_list, marker_colors=otu_marker_colors, otu_label_off=otu_label_off, label_mapper=label_mapper)

    ax.set_xlabel(r'DESeq2 (Log2 fold change)', fontsize=14)

    if len(otu_groups) > 1:
        _group_bg_colorer(axes[0], otu_groups, show_label=False)
        _group_bg_colorer(axes[1], otu_groups, show_label=True)

    from ..visualizers.PlotTools import make_legends

    handles, labels = make_legends(legend_configs={
        handle_name[0]: {'marker': '', 'color': otu_marker_colors[0]},
        handle_name[2]: {'marker': 'o', 'color': '#1C1C1C', 'type': 'scatter'},
        handle_name[1]: {'marker': '', 'color': otu_marker_colors[2]},
        handle_name[3]: {'marker': '^', 'color': '#1C1C1C', 'type': 'scatter'},
        'Not significant': {'marker': '', 'color': otu_marker_colors[1]}
    })
    ax.legend(handles=handles, labels=labels,
              bbox_to_anchor=(0.2, 0.875, 0.6, 0.5), loc='lower center',
              ncol=min(len(handles), 3), frameon=False, fontsize=12, bbox_transform=fig.transFigure)
    fig.align_xlabels(axes)
    plt.show()
    return otu_list


def phylo_tree_model_res(compare_sets=None, use_phylo_names=True, legend_label=None, color=None,
                         path_to_tree=None, otu_list=None, abun_otu_only=True, figsize=None):
    from ..data import BayesianRes, DeseqResult

    if compare_sets is None:
        compare_sets = [{'dataset': 'filtered', 'model': BayesianRes,
                         'compare': 'pre_vs_skin', 'label': 'With Patient 16, BGLMM'},
                        {'dataset': 'filtered_no_p16', 'model': BayesianRes,
                         'compare': 'pre_vs_skin', 'label': 'Without Patient 16, BGLMM'},
                        {'dataset': 'filtered', 'model': DeseqResult,
                         'compare': 'pre_vs_skin', 'label': 'With Patient 16, DESeq2'},
                        {'dataset': 'filtered_no_p16', 'model': DeseqResult,
                         'compare': 'pre_vs_skin', 'label': 'Without Patient 16, DESeq2'}]
    if use_phylo_names:
        from skin_mb.data.table import get_otu_id_to_taxo_name
        otu_name_mapper = get_otu_id_to_taxo_name().to_dict()
    else:
        otu_name_mapper = None
    if legend_label is None:
        legend_label = ['Enriched in wound', 'Enriched in skin']
    if color is None:
        color = ['#B2112A', '#2C73B4']
    if path_to_tree is None:
        from os import environ
        if 'DATA_PATH' not in environ:
            raise EnvironmentError('Environmental variable DATA_PATH does not exist, please indicate path_to_tree')
        else:
            path_to_tree = environ['DATA_PATH'] + '/otu_tables/rep_set.tre'

    if otu_list is None:
        otu_list = set()
        for compare_set in compare_sets:
            otu_list.update(list(
                compare_set['model'].load_default(compare_set['dataset']).get_data(compare_set['compare']).sig_otu.index
            ))
        otu_list = list(otu_list)

    if abun_otu_only:
        from ..data import OtuTable
        otu_table = OtuTable.load_default('filtered')
        otu_scope = otu_table.otu_stats()[otu_table.otu_stats()['rel_abun_avg_on_all'] >= 0.001].index
        otu_list = [otu for otu in otu_list if otu in otu_scope]

    from ..visualizers.DendroPlot import dendro_plot, sig_otu_chart
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.style.use('default')

    if figsize is None:
        figsize = [18, 6]
    fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=False)
    fig.subplots_adjust(hspace=0, wspace=0)

    otu_reordered, chart_lim, pos = dendro_plot(ax=axes[0], path_to_tree=path_to_tree, otu_list=otu_list)
    sig_otu_chart(ax=axes[1], pos=pos, otu_list=otu_reordered, chart_lim=chart_lim, compare_sets=compare_sets,
                  label_mapper=otu_name_mapper, color=color, legend_label=legend_label)

    plt.show()


def sig_otu_abun_bar_plot(dataset='filtered_no_p16', sample_to_show=None):

    import pandas as pd
    from ..data import OtuTable, DeseqResult, BayesianRes

    otu_table = OtuTable.load_default(dataset=dataset)
    if sample_to_show is None:
        sample_to_show = otu_table.sample_list

    deseq2 = DeseqResult.load_default(dataset=dataset)
    bglmm = BayesianRes.load_default(dataset=dataset)

    sig_otu_commom = set(deseq2.pre_vs_skin.sig_otu.index.values).intersection(set(bglmm.pre_vs_skin.sig_otu.index.values))
    sig_otu_deseq2 = set(deseq2.pre_vs_skin.sig_otu.index.values).difference(sig_otu_commom)
    sig_otu_bglmm = set(bglmm.pre_vs_skin.sig_otu.index.values).difference(sig_otu_commom)

    data_to_plot = pd.DataFrame(index=sample_to_show)
    data_to_plot['all_otu'] = otu_table.count_table[sample_to_show].sum(axis=0)
    data_to_plot['DESeq2 only'] = otu_table.count_table.loc[sig_otu_deseq2][sample_to_show].sum(axis=0)/ \
                                            data_to_plot['all_otu']
    data_to_plot['BGLMM only'] = otu_table.count_table.loc[sig_otu_bglmm][sample_to_show].sum(axis=0) / \
                                            data_to_plot['all_otu']
    data_to_plot['Both'] = otu_table.count_table.loc[sig_otu_commom][sample_to_show].sum(axis=0) / \
                                            data_to_plot['all_otu']
    data_to_plot['Non-significant'] = 1 - data_to_plot['DESeq2 only'] - data_to_plot['BGLMM only'] - data_to_plot['Both']
    data_to_plot = data_to_plot[['DESeq2 only', 'Both', 'BGLMM only', 'Non-significant']] * 100

    def patient_id_corrector(pid):
        return str(int(pid[:-1]) - 20) + pid[-1]

    data_to_plot.rename(index=patient_id_corrector, inplace=True)

    import matplotlib.pyplot as plt
    from ..visualizers.ThreePanelBarplot import three_panel_barplot
    three_panel_barplot(data_to_plot=data_to_plot, percent=True,
                        stacked=True, width=0.8, color=['#2C73B4', '#F39730', '#B2112A',  '#AEAEAE'], ncol=4)
    plt.show()