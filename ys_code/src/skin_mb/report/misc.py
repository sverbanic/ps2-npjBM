
def _read_pipeline_history(count_files_root=None):
    """Read tab files of reads in each steps"""

    # import raw reads for each sample (forward)
    def read_tab_raw(file_path):
        """Read raw read count info from count tab files, only forward reads are recorded
        """

        import pandas as pd
        sample_reads = pd.read_csv(file_path, sep='\t', index_col=0)
        ix_fwd = [ix for ix in sample_reads.index if 'R1' in ix]
        name_mapper = {ix: ix.split('_')[0] for ix in ix_fwd}
        return sample_reads.loc[ix_fwd].rename(index=name_mapper)

    # import trimmed reads for each sample
    def read_tab_trimmed(file_path):
        """Read trimmed count info from count tab file, only forward reads are recorded
        """

        import pandas as pd
        sample_reads = pd.read_csv(file_path, sep='\t', index_col=0)
        ix_fwd = [ix for ix in sample_reads.index if 'R1' in ix]
        name_mapper = {ix: ix.split('_')[0] for ix in ix_fwd}
        return sample_reads.loc[ix_fwd].rename(index=name_mapper)

    # import joined reads for each sample
    def read_tab_joined(file_path):
        """Read joined (fastq-join) count info from count tab file
        """

        import pandas as pd
        sample_reads = pd.read_csv(file_path, sep='\t', index_col=0)
        name_mapper = {ix: ix.split('_')[0] for ix in sample_reads.index}
        return sample_reads.rename(index=name_mapper)

    if count_files_root is None:
        from os import environ
        if 'DATA_PATH' not in environ:
            raise EnvironmentError('Environmental variable DATA_PATH not found, please indicate count_files_roots')
        else:
            count_files_root = environ['DATA_PATH'] + '/seq_prep/'

    # import OTU picked reads for each sample
    from ..data import OtuTable

    read_history = read_tab_raw(count_files_root + '/rawread.counts.tab').rename(columns={'Reads': 'raw'})
    read_history.loc[:, 'trimmed'] = read_tab_trimmed(count_files_root + '/trimmedread.counts.tab')['Reads']
    read_history.loc[:, 'joined'] = read_tab_joined(count_files_root + '/joinedread.counts.tab')['Reads']
    read_history = read_history.rename(index={'MC7-5': 'MC7.5', 'Mcplus': 'MCplus'})
    otu_table = OtuTable.load_default('full_table')
    read_history.loc[:, 'otu_picking'] = otu_table.count_table.sum(axis=0).rename({'Mcplus': 'MCplus'})

    return read_history


def passing_rate_history_plot(count_files_root=None):
    import matplotlib.pyplot as plt

    read_history = _read_pipeline_history(count_files_root=count_files_root)

    fig = plt.figure(figsize=[10, 5])
    ax = fig.add_subplot(111)

    def get_color(row):
        if row.name[-1] == 'A':
            return '#2C73B4'
        elif row.name[-1] == 'B':
            return '#F39730'
        elif row.name[-1] == 'C':
            return '#1C7725'
        else:
            return '#AEAEAE'

    def plot_single(row):
        x = [0, 1, 2, 3]
        y = (row / row['raw'] * 100).values
        ax.plot(x, y, get_color(row), marker='o', markersize=2, lw=1, alpha=0.5)

        if y[2] < 20:
            ax.text(s=row.name, x=2, y=y[2] * 1.02, ha='right', va='bottom', fontsize=12)

    read_history.apply(plot_single, axis=1)

    lgd = [ax.plot([], [], '#2C73B4', marker='o', markersize=2, lw=1, label='Pre-debridement'),
           ax.plot([], [], '#F39730', marker='o', markersize=2, lw=1, label='Post-debridement'),
           ax.plot([], [], '#1C7725', marker='o', markersize=2, lw=1, label='Skin'),
           ax.plot([], [], '#AEAEAE', marker='o', markersize=2, lw=1, label='Controls')]

    labels = ['Pre-debridement', 'Post-deberidement', 'Skin', 'Controls']
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(labels=['Raw reads', 'Trimmed reads', 'Joined reads', 'Reads passed\nOTU picking'], fontsize=14,
                       rotation=0, ha='center')
    ax.set_ylabel('Percent of reads (%)', fontsize=14)
    ax.legend(frameon=False, fontsize=14)
    plt.show()


def otu_distribution_scatter_plot(dataset='full_table'):
    from ..data import OtuTable
    otu_table = OtuTable.load_default(dataset=dataset)

    otu_table.filter(axis='sample', ids_to_keep=otu_table.sample_list[:60], remove_empty=True, inplace=True)

    otu_abun = otu_table.otu_stats()['rel_abun_avg_on_all']
    otu_prev = otu_table.otu_stats()['smpl_detected']

    import matplotlib.pyplot as plt
    import numpy as np

    fig, axes = plt.subplots(2, 2, figsize=[8, 6],
                         gridspec_kw={'width_ratios': (5, 1),
                                      'height_ratios': (1, 3)})
    fig.subplots_adjust(hspace=0, wspace=0)
    axes[0,1].set_axis_off()

    axes[1,0].scatter(otu_abun, otu_prev, color='#2C73B4', alpha=0.2, s=5)
    axes[1,0].set_xscale('log')
    axes[1,0].set_xlim([np.min(otu_abun), np.max(otu_abun)])
    axes[1,0].set_ylim([0, 61])

    bins = np.logspace(np.log10(np.min(otu_abun)), np.log10(np.max(otu_abun)), 35)
    axes[0,0].hist(otu_abun, bins=bins, color='#2C73B4', edgecolor='white')
    axes[0,0].set_xlim([np.min(otu_abun), np.max(otu_abun)])
    axes[0,0].set_xscale('log')
    axes[0,0].set_xticks([])
    axes[0,0].set_yticks([])
    axes[0,0].spines['right'].set_visible(False)
    axes[0,0].spines['top'].set_visible(False)
    axes[0,0].spines['left'].set_visible(False)
    axes[0,0].set_facecolor('#FFFFFF')

    bins = np.linspace(0, 61, 62)
    axes[1,1].hist(otu_prev, bins=bins, color='#2C73B4', edgecolor='white', align='left', orientation='horizontal')
    axes[1,1].set_ylim([0, 61])
    axes[1,1].set_xticks([])
    axes[1,1].set_yticks([])
    axes[1,1].spines['right'].set_visible(False)
    axes[1,1].spines['top'].set_visible(False)
    axes[1,1].spines['bottom'].set_visible(False)
    axes[1,1].set_facecolor('#FFFFFF')

    axes[1,0].set_xlabel('Average relative abundance', fontsize=14)
    axes[1,0].set_ylabel('Number of samples detected', fontsize=14)

    plt.show()


def per_sample_top_otu_curve():
    def survey_unique_vs_abun(sample_id, otu_table, norm_otu_number=False):
        import numpy as np

        otu_counts = sorted(otu_table.count_table[sample_id], reverse=True)
        unique_abun_list = []
        if norm_otu_number:
            otu_num = np.sum(np.array(otu_counts) > 0)
            total_counts = np.sum(otu_counts)
            for otu_ix in range(otu_num):
                unique_abun_list.append([otu_ix / otu_num, np.sum(otu_counts[:otu_ix]) / total_counts])
        else:
            total_counts = np.sum(otu_counts)
            for otu_ix in range(len(otu_counts)):
                if otu_counts[otu_ix] > 0:
                    unique_abun_list.append([otu_ix, np.sum(otu_counts[:otu_ix] / total_counts)])
                else:
                    unique_abun_list.append([otu_ix, 1])
        return np.array(unique_abun_list).T

    from ..data import OtuTable
    import matplotlib.pyplot as plt

    otu_table = OtuTable.load_default(dataset='full_table')
    otu_table.filter(axis='sample', ids_to_keep=otu_table.sample_list[:60], remove_empty=True, inplace=True)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    for sample in otu_table.sample_list[:60]:
        unique_abun = survey_unique_vs_abun(sample, otu_table, norm_otu_number=False)
        if 'A' in sample or 'B' in sample:
            ax.plot(unique_abun[0], unique_abun[1], ls='-', color='#FC820D', lw=0.5, alpha=0.5)
        else:
            ax.plot(unique_abun[0], unique_abun[1], ls='-', color='#2C73B4', lw=0.5, alpha=0.5)

    ax.plot([], [], ls='-', color='#FC820D', lw=0.5, label='Wound samples')
    ax.plot([], [], ls='-', color='#2C73B4', lw=0.5, label='Skin samples')

    ax.set_xlim([-50, 1000])
    ax.set_xlabel('Number of most abundant unique OTUs', fontsize=14)
    ax.set_ylabel('Accumulated relative abundance', fontsize=14)
    ax.text(s='Patient 16', x=800, y=0.45, fontsize=12)
    ax.legend(frameon=False, fontsize=12)
    plt.show()


def survey_blast_results_all_sample(root=None):

    def plot_stack_bars(result_1, result_2, y, ax):
        label_formatter = lambda label, bar: label + ' ({:,})'.format(bar)
        get_length = lambda otu_set: len(otu_set)

        set_1 = set([query[0] for query in result_1.items() if query[1]['hit_num'] > 0])
        set_2 = set([query[0] for query in result_2.items() if query[1]['hit_num'] > 0])
        set_12 = set_1 & set_2
        bar_0 = get_length(set(result_1.keys()) - (set_1 | set_2))
        bar_1 = get_length(set_1 - set_12)
        bar_2 = get_length(set_12)
        bar_3 = get_length(set_2 - set_12)
        ax.barh(y=[y], width=[bar_0], left=[0],
                align='center', height=0.4, color='#AEAEAE',
                label=label_formatter('No hits', bar_0))
        ax.barh(y=[y], width=[bar_1], left=[bar_0],
                align='center', height=0.4, color='#2C73B4',
                label=label_formatter('Hits only in 16s rRNA', bar_1))
        ax.barh(y=[y], width=[bar_2], left=[bar_0 + bar_1],
                align='center', height=0.4, color='#1C7725',
                label=label_formatter('Hits in both', bar_2))
        ax.barh(y=[y], width=[bar_3], left=[bar_0 + bar_1 + bar_2],
                align='center', height=0.4, color='#FC820D',
                label=label_formatter('Hits only in human genome', bar_3))
        ax.set_xlim([0, len(result_1)])

    from os import environ

    if root is None:
        if 'DATA_PATH' not in environ:
            raise EnvironmentError('Environmental variable DATA_PATH is not found, please indicate the root to blastn results')
        else:
            root = environ['DATA_PATH'] + '/blastn/'
    blast_res = ['fj_silva_aligned_16s_trimmed.json', 'fj_silva_aligned_human_trimmed.json',
                 'fj_silva_failed_16s_trimmed.json', 'fj_silva_failed_human_trimmed.json']

    from ..data.blastn import survey_blast_result

    results = [survey_blast_result(root + res) for res in blast_res]
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, figsize=[12, 5])
    plot_stack_bars(result_1=results[0], result_2=results[1], y=1, ax=axes[0])
    axes[0].text(s='OTUs aligned', x=-0.004 * (axes[0].get_xlim()[1] - axes[0].get_xlim()[0]), y=1, ha='right',
                 va='center', rotation=90, fontsize=14)
    axes[0].axis('off')
    axes[0].legend(loc=(-0.01, -0.2), ncol=4, frameon=False, fontsize=11)
    plot_stack_bars(result_1=results[2], result_2=results[3], y=1, ax=axes[1])
    axes[1].text(s='OTUs failed', x=-0.004 * (axes[1].get_xlim()[1] - axes[1].get_xlim()[0]), y=1, ha='right',
                 va='center', rotation=90, fontsize=14)
    axes[1].axis('off')
    axes[1].legend(loc=(-0.01, -0.2), ncol=4, frameon=False, fontsize=11)
    plt.show()


def filter_table_barplot(figsize=None):
    from ..data import OtuTable

    if figsize is None:
        figsize= [12, 8]

    full_table = OtuTable.load_default(dataset='full_table').count_table.sum(axis=0)
    pre_filtered_table = OtuTable.load_default(dataset='filtered').count_table.sum(axis=0)
    full_table = full_table[pre_filtered_table.index]

    import pandas as pd
    data_to_plot = pd.DataFrame(index=full_table.index)
    data_to_plot.loc[:, 'In filtered table'] = pre_filtered_table/full_table * 100
    data_to_plot.loc[:, 'Not in filtered table'] = 100 - pre_filtered_table/full_table * 100

    def patient_id_corrector(pid):
        return str(int(pid[:-1]) - 20) + pid[-1]
    data_to_plot.rename(index=patient_id_corrector, inplace=True)

    from ..visualizers import ThreePanelBarplot
    ThreePanelBarplot.three_panel_barplot(data_to_plot, width=0.8, figsize=figsize,
                                  colors=['#2C73B4', '#AEAEAE'], stacked=True,
                                  bbox_to_anchor=(0, -0.3, 1, 0.2))
    import matplotlib.pyplot as plt
    plt.show()


def add_break(ax, x, y, length=3.0, rotation=45, color='#151515', **kwargs):
    """length and width is currently on data coordinate"""
    import numpy as np
    if y == 0:
        rotation = 90 - 45
    ax.plot([x - np.sin(rotation) * length / 2, x + np.sin(rotation / 180 * np.pi) * length / 2],
            [y - np.cos(rotation) * length / 2, y + np.cos(rotation / 180 * np.pi) * length / 2], color=color,
            **kwargs)


def add_break_axis(ax, x=None, y=None, break_pos=0.4, break_space=0.1, color='#151515', lw=1):
    if x is not None:
        ax.plot([x[0], break_pos - break_space / 2], [0, 0], color=color, lw=lw, zorder=1)
        ax.plot([break_pos + break_space / 2, x[1]], [0, 0], color=color, lw=lw, zorder=1)
        add_break(ax, x=break_pos - break_space / 2, y=0, lw=lw, length=0.1, zorder=1, color=color)
        add_break(ax, x=break_pos + break_space / 2, y=0, lw=lw, length=0.1, zorder=1, color=color)

    if y is not None:
        ax.plot([0, 0], [y[0], break_pos - break_space / 2], color=color, lw=lw, zorder=1)
        ax.plot([0, 0], [break_pos + break_space / 2, y[1]], color=color, lw=lw, zorder=1)
        add_break(ax, x=0, y=break_pos - break_space / 2, lw=lw, length=0.1, zorder=1, color=color)
        add_break(ax, x=0, y=break_pos + break_space / 2, lw=lw, length=0.1, zorder=1, color=color)


def bglmm_prediction(group='pre', data_root=None, fold_upper=10, exclude_zeros=False):

    import numpy as np
    import pandas as pd
    from ..data import OtuTable
    import matplotlib.pyplot as plt

    if data_root is None:
        from os import environ
        if 'DATA_PATH' not in environ:
            raise EnvironmentError(
                'Environmental variable DATA_PATH is not found, please indicate the root to blastn results')
        else:
            data_root = environ['DATA_PATH']

    dataset = 'filtered_no_p16'
    otu_table = OtuTable.load_default(dataset)
    prediction = pd.read_csv(data_root + '/bglmm/bglmm_pred/{}-median.txt'.format(dataset), sep='\t',header=None).values
    prediction = pd.DataFrame(prediction, columns=otu_table.otu_list, index=otu_table.sample_list).T

    label, label_name = {'pre': ('A', 'Pre-debridement'),
                         'post': ('B', 'Post-debridement'),
                         'skin': ('C', 'Skin')}[group]

    col_to_select = [col for col in prediction.columns if label in col]
    pred = prediction[col_to_select].values.flatten()
    truth = otu_table.count_table[col_to_select].values.flatten()
    flags = np.repeat(True, len(pred))
    if exclude_zeros:
        flags = (pred > 0) & (truth > 0) & flags
        print("Excluding zero passing rate: {}".format(np.sum(flags) / len(flags)))
    if fold_upper is not None:
        flags = ((pred <= fold_upper * truth) | (truth == 0) | (pred == 0)) & flags
        print("Fold upper bound passing rate: {}".format(np.sum(flags) / len(flags)))

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(truth[flags], pred[flags], alpha=0.5, s=5, color='#2C73B4', zorder=2)
    ax.scatter(truth[~flags], pred[~flags], alpha=0.5, s=5, color='#AEAEAE', zorder=2)
    ax.plot([0, 1e6], [0, 1e6], color='#595959', lw=1, ls='--', zorder=0)

    add_break_axis(ax=ax, x=[0, 1e6], y=[0, 1e6], break_pos=0.4, break_space=0.1, color='#595959')
    ax.set_xlim([-0.1, 1e6])
    ax.set_ylim([-0.1, 1e6])
    ax.set_xscale('symlog', linthreshx=1)
    ax.set_yscale('symlog', linthreshy=1)
    ax.set_xlabel('True counts (log scale)', fontsize=12)
    ax.set_ylabel('Predicted counts (log scale)', fontsize=12)

    from scipy.stats import pearsonr, spearmanr
    ax.text(
        s='Pearson = %.5f (%.5f, pred < %i$\\times$truth)\nSpearman = %.5f (%.5f, pred < %i$\\times$truth)' % (
            pearsonr(pred, truth)[0],
            pearsonr(pred[flags], truth[flags])[0],
            fold_upper,
            spearmanr(pred, truth)[0],
            spearmanr(pred[flags], truth[flags])[0],
            fold_upper
        ),
        x=0.3, y=3e5, ha='left', va='top',
        color='#565656', fontsize=10)

    ax.text(s=label_name, x=0.3, y=6e5, ha='left', va='center', color='#151515', fontsize=12)

    plt.show()