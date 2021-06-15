from .result import Result


class DeseqResultSingle:

    def __str__(self):
        import numpy as np
        return 'DESeq2 single compare: {} total OTUs, {} positive, {} negative, {} not significant'.format(
            self.full_res.shape[0],
            np.sum(self.sig_otu['log2FoldChange'] > 0),
            np.sum(self.sig_otu['log2FoldChange'] < 0),
            self.full_res.shape[0] - self.sig_otu.shape[0]
        )

    def __init__(self, csv_path, header=0, index=0, p_cutoff=0.05, p_col='padj', **kwargs):
        import pandas as pd
        self.full_res = pd.read_csv(csv_path, header=header, index_col=index)
        self.sig_otu = self.full_res[self.full_res[p_col] <= p_cutoff]

    def get_data(self, otu_list=None, get_all=False):
        if otu_list is None:
            data = self.sig_otu
        else:
            if isinstance(otu_list, str):
                otu_list = [otu_list]
            data = self.full_res.loc[otu_list]
        if get_all:
            return data
        else:
            data = data[['log2FoldChange', 'lfcSE']]
            data.loc[:, 'mid'] = data['log2FoldChange']
            data.loc[:, 'lower'] = data['log2FoldChange'] - 1.96 * data['lfcSE']
            data.loc[:, 'upper'] = data['log2FoldChange'] + 1.96 * data['lfcSE']

            def sig_mapper(row):
                if row.name in self.sig_otu.index:
                    return 1 if row['mid'] > 0 else -1
                else:
                    return 0

            data.loc[:, 'sig'] = data.apply(sig_mapper, axis=1)
            return data[['mid', 'lower', 'upper', 'sig']]

    def cutoff_pval(self, p_cutoff, p_col='padj'):
        self.sig_otu = self.full_res[self.full_res[p_col] <= p_cutoff]


class DeseqResult(Result):

    def __init__(self, csv_paths=None, name=None, **kwargs):
        super().__init__()
        if csv_paths is None:
            self.load_default(**kwargs)
        elif isinstance(csv_paths, dict):
            self.__dict__.update({
                data_id: DeseqResultSingle(data_path, **kwargs)
                for data_id, data_path in csv_paths.items()
            })
        else:
            raise TypeError('`csv_paths` should be dictionary')
        self.name = name

    def get_data(self, name):
        return self.__dict__[name]

    def get_significance(self, compare, otu_list):
        dataset = self.get_data(name=compare)

        def significance_mapper(row):
            if row['sig_otu']:
                if row['mid'] < 0:
                    return -1
                else:
                    return 1
            else:
                return 0
        return dataset.get_data(otu_list=otu_list).apply(significance_mapper, axis=1)

    @classmethod
    def from_dataset(cls, dataset_path, content_map=None, name=None, **kwargs):

        if content_map is None:
            content_map = {
                'healed_vs_unhealed': 'healed.vs.unhealed.overall.csv',
                'healed_vs_unhealed_pre': 'healed.vs.unhealed.pre.csv',
                'healed_vs_unhealed_post': 'healed.vs.unhealed.post.csv',
                'pre_vs_post': 'pre.vs.post.overall.csv',
                'pre_vs_post_healed': 'pre.vs.post.healed.csv',
                'pre_vs_post_unhealed': 'pre.vs.post.unhealed.csv',
                'pre_vs_skin': 'pre.vs.skin.csv',
                'wound_vs_skin': 'wound.vs.skin.overall.csv',
            }

        return cls(csv_paths={data_id: dataset_path + '/results/' + file_name for data_id, file_name in content_map.items()},
                   name=name, **kwargs)

    def scatter_plot_with_bar(self, data_to_plot, plot_config=None, ax=None, otu_list=None, plot_scope=None,
                              pos_mapper=None, vertical=True, value_label=r'DESeq2 ($\log_2$ fold change)', label_mapper=None,
                              otu_label_off=False, legend_config=None, figsize=None, **kwargs):

        if isinstance(data_to_plot, str):
            data_to_plot = [self.get_data(data_to_plot)]
        elif isinstance(data_to_plot, list):
            data_to_plot = {data: self.get_data(data) for data in data_to_plot}

        from ..visualizers.ScatterErrorBar import multi_series
        return multi_series(series_to_plot=data_to_plot, ax=ax, plot_config=plot_config, entry_list=otu_list,
                            plot_scope=plot_scope, pos_mapper=pos_mapper, figsize=figsize, vertical=vertical,
                            label_mapper=label_mapper, value_label=value_label, entry_label_off=otu_label_off,
                            legend_config=legend_config, **kwargs)

    @classmethod
    def load_default(cls, dataset=None, data_root=None):
        """

        Args:
            dataset: select from 'filtered', 'filtered_no_p16'

        Returns:

        """
        from os import environ
        if data_root is None:
            if 'DATA_PATH' not in environ:
                raise EnvironmentError('Please indicate the root to data folder')
            else:
                root = environ['DATA_PATH']
        else:
            root = data_root

        if dataset is None:
            dataset = 'filtered'

        return cls.from_dataset(dataset_path=root + 'deseq2/' + dataset, name=dataset)
