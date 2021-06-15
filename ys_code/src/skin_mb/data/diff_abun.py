from .result import Result
import numpy as np
import pandas as pd


class DiffAbunRes(Result):

    def __init__(self, otu_table, transform_pipe=None, percent=False, **kwargs):
        super().__init__()

        self.pre_vs_skin = diff_rel_abun(otu_table, compare='pre_vs_skin', transform_pipe=transform_pipe,
                                         percent=percent, **kwargs)
        self.post_vs_skin = diff_rel_abun(otu_table, compare='post_vs_skin', transform_pipe=transform_pipe,
                                         percent=percent, **kwargs)
        self.pre_vs_post = diff_rel_abun(otu_table, compare='pre_vs_post', transform_pipe=transform_pipe,
                                         percent=percent, **kwargs)

    @classmethod
    def load_default(cls, dataset='filtered', transform_pipe=None, percent=False, data_root=None, **kwargs):
        from os import environ
        if data_root is None:
            if 'DATA_PATH' not in environ:
                raise EnvironmentError('Please indicate the root to data folder')
            else:
                root = environ['DATA_PATH']
        else:
            root = data_root

        from skin_mb.data import OtuTable
        return cls(otu_table=OtuTable.from_biom(root + 'otu_tables/' + dataset + '.biom'),
                  transform_pipe=transform_pipe, percent=percent, **kwargs)

    def heatmap(self, compare, ax=None, otu_list=None, subject_list=None,
                adjust_patient_id=True, label_mapper=None, z_log=False, **kwargs):
        dataset = self.get_data(compare)
        if otu_list is None:
            otu_list = dataset.index
        else:
            dataset = dataset.loc[otu_list]

        if subject_list is None:
            subject_list = dataset.columns

        dataset = dataset[subject_list]

        if ax is None:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(len(subject_list)/5, len(otu_list)/5))
            ax = fig.add_subplot(111)

        if adjust_patient_id:
            id_mapper = {col: col-20 for col in dataset.columns}
            dataset.rename(columns=id_mapper, inplace=True)

        if label_mapper:
            dataset.rename(index=label_mapper, inplace=True)

        if 'figsize' in kwargs.keys():
            _ = kwargs.pop('figsize')
        if 'gridspec_kw' in kwargs.keys():
            _ = kwargs.pop('gridspec_kw')

        from ..visualizers.Heatmap import heatmap
        from ..visualizers.PlotTools import CMapsDi
        cmap = CMapsDi.BluWhtRed(reverse=True)
        im = heatmap(values=dataset, ax=ax, origin='lower', z_log=z_log, zorder=5, nan_color='#AEAEAE', cmap=cmap, **kwargs)
        ax.set_xlabel('Patient', fontsize=14)
        return im, dataset


def rel_abun(table):
    return table.apply(lambda row: row/row.sum(), axis=0)


def clr(table, pseudo_count=0.1):
    def clr_row(row):
        from scipy.stats import gmean
        return np.log((row + pseudo_count)/gmean(row + pseudo_count))
    return table.apply(clr_row, axis=0)


def diff(table, subject_list, postfix1, postfix2, percent=False):
    diff_mtx = pd.DataFrame(index=table.index, columns=subject_list)
    for subject in subject_list:
        if percent:
            diff_mtx[subject] = (table[str(subject) + postfix1] - table[str(subject) + postfix2])/table[str(subject) + postfix1]
        else:
            diff_mtx[subject] = table[str(subject) + postfix1] - table[str(subject) + postfix2]
    return diff_mtx


def diff_rel_abun(otu_table, compare='wound_vs_skin', transform_pipe=None, pseudo_count=0.1,
                  percent=False, otu_list=None, subject_list=None):

    if subject_list is None:
        subject_list = set([int(sample[:2]) for sample in otu_table.sample_list])
    if otu_list is None:
        otu_list = otu_table.otu_list

    if compare.lower() in ['pre_vs_skin']:
        postfix1 = 'A'
        postfix2 = 'C'
    elif compare.lower() in ['pre_vs_post']:
        postfix1 = 'A'
        postfix2 = 'B'
    elif compare.lower() in ['post_vs_skin']:
        postfix1 = 'B'
        postfix2 = 'C'
    else:
        raise ValueError("Compare should be 'pre_vs_skin', 'pre_vs_post', or 'post_vs_skin'")

    if transform_pipe is None:
        transform_pipe = ['rel abun', 'diff']

    from functools import partial

    transformations = {
        'rel abun': rel_abun,
        'clr': partial(clr, pseudo_count=pseudo_count),
        'diff': partial(diff, subject_list=subject_list,
                        postfix1=postfix1, postfix2=postfix2, percent=percent)
    }

    table = otu_table.count_table.loc[otu_list]
    for transform in transform_pipe:
        table = transformations[transform](table)
    return table

