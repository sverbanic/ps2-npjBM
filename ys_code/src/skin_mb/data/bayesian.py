"""
Functions to analyze Bayesian modeling results
"""
from .result import Result


class BayesResSingle:
    """Class storing results for single comparison (\beta, e.g. pre_vs_skin)"""

    def __str__(self):
        import numpy as np
        return 'BGLMM single compare: {} total OTUs, {} positive, {} negative, {} not significant'.format(
            self.full_res.shape[0],
            np.sum(self.sig_otu['mid'] > 0),
            np.sum(self.sig_otu['mid'] < 0),
            self.full_res.shape[0] - self.sig_otu.shape[0]
        )

    def __init__(self, beta_table, prefix, flip=False):
        self.full_res = beta_table[[prefix + '_Lower', prefix + '_Mid', prefix + '_Upper']].rename(
            columns={
                prefix + '_Lower': 'lower',
                prefix + '_Mid': 'mid',
                prefix + '_Upper': 'upper'
            }
        )
        if flip:
            self.full_res = -self.full_res
            self.full_res.rename(columns={'lower': 'upper', 'upper': 'lower'}, inplace=True)
        self.sig_otu = self.full_res[(self.full_res['lower'] > 0) | (self.full_res['upper'] < 0)]

    def get_data(self, otu_list=None):
        if otu_list is None:
            data = self.sig_otu
        else:
            if isinstance(otu_list, str):
                otu_list = [otu_list]
            data = self.full_res.loc[otu_list]

        def significance_mapper(row):
            if row.name in self.sig_otu.index:
                if row['mid'] < 0:
                    return -1
                else:
                    return 1
            else:
                return 0

        data.loc[:, 'sig'] = data.apply(significance_mapper, axis=1)
        return data[['mid', 'lower', 'upper', 'sig']]


class BayesianRes(Result):

    def __init__(self, path_to_betas, path_to_otu_modeled=None, path_to_alpha_hat=None, path_to_r_hat=None,
                 betas_tab='beta-all', otu_modeled_tab='list of OTUs', healed_vs_unhealed=None, name=None, note=None,
                 use_legacy=False):
        super().__init__()
        if path_to_betas is not None:
            beta_table, _ = self._load_betas(xlsx_path=path_to_betas, tab=betas_tab, use_legacy=use_legacy)
            _, self.sample_modeled = self._load_betas(xlsx_path=path_to_betas, tab='beta', index_col='OUT_name',
                                                      use_legacy=use_legacy)
            self.pre_vs_skin = BayesResSingle(beta_table=beta_table, prefix='Pre')
            self.post_vs_skin = BayesResSingle(beta_table=beta_table, prefix='Post')
            self.pre_vs_post = BayesResSingle(beta_table=beta_table, prefix='Diff')
            if healed_vs_unhealed is not None:
                beta_table, _ = self._load_betas(xlsx_path=healed_vs_unhealed,
                                                 tab=betas_tab,
                                                 healed_vs_unhealed=True,
                                                 use_legacy=use_legacy)
                self.pre_vs_post_unhealed = BayesResSingle(beta_table=beta_table, prefix='diff1', flip=True)
                self.pre_vs_post_healed = BayesResSingle(beta_table=beta_table, prefix='diff2', flip=True)
                self.healed_vs_unhealed_pre = BayesResSingle(beta_table=beta_table, prefix='diff3')
                self.healed_vs_unhealed_post = BayesResSingle(beta_table=beta_table, prefix='diff4')
            else:
                self.pre_vs_post_unhealed = None
                self.pre_vs_post_healed = None
                self.unhealed_vs_healed_pre = None
                self.unhealed_vs_healed_post = None

        if path_to_otu_modeled is not None:
            self.otu_modeled = self._load_otu_modeled(xlsx_path=path_to_otu_modeled, tab=otu_modeled_tab,
                                                      use_legacy=use_legacy)

        if path_to_alpha_hat is not None:
            self.alpha_hat = self._load_alpha_hat(file_path=path_to_alpha_hat, use_legacy=use_legacy)

        if path_to_r_hat is not None:
            self.r_hat = self._load_r_hat(file_path=path_to_r_hat, use_legacy=use_legacy)

        self.name = name
        self.note = note

    def scatter_plot_with_bar(self, data_to_plot, plot_config=None, ax=None, otu_list=None, plot_scope=None,
                              pos_mapper=None, vertical=True, value_label=r'BGLMM ($\beta$)', label_mapper=None,
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

    @staticmethod
    def _otu_name_parser(compound):
        return compound.split(',')[0]

    def _load_betas(self, xlsx_path, tab='beta-all', index_col='OTU name', clean_index=True,
                    healed_vs_unhealed=False, include_cols=None, use_legacy=False):
        """parse beta values from Bayesian model results in excel tables

        Args:
            xlsx_path (`str`): directory to `.xlsx` file of bayesian model results
            tab (`str`): name of the tab in the file containing the results, default 'beta'
            index_col (`str`): Column of OTU name, default 'OUT_name'
            include_cols (list of str): list of column names to include, include all if None

        Returns: a `pd.DataFrame` containing the model results

        """
        import pandas as pd
        if use_legacy:
            tab = 'beta'
            index_col = 'OUT_name'

        if include_cols is None:
            res_df = pd.read_excel(xlsx_path, sheet_name=tab, index_col=index_col, header=0)
        else:
            res_df = pd.read_excel(xlsx_path, sheet_name=tab, index_col=index_col, usecols=include_cols, header=0)
        if clean_index:
            res_df.index = res_df.index.map(self._otu_name_parser)

        if healed_vs_unhealed:
            beta_cols = ['diff1_Mid', 'diff1_Lower', 'diff1_Upper',
                         'diff2_Mid', 'diff2_Lower', 'diff2_Upper',
                         'diff3_Mid', 'diff3_Lower', 'diff3_Upper',
                         'diff4_Mid', 'diff4_Lower', 'diff4_Upper']
            sample_names = pd.Series(data=[col for col in res_df.columns if (col != 'No') and (col not in beta_cols)],
                                     name='sample_modeled')
        else:
            beta_cols = ['Pre_Mid', 'Pre_Lower', 'Pre_Upper',
                         'Post_Mid', 'Post_Lower', 'Post_Upper',
                         'Diff_Mid', 'Diff_Lower', 'Diff_Upper']
            sample_names = pd.Series(data=[col for col in res_df.columns if (col != 'No') and (col not in beta_cols)],
                                     name='sample_modeled')
        return res_df[beta_cols], sample_names

    def _load_otu_modeled(self, xlsx_path, tab='list of OTUs', use_legacy=False):
        import pandas as pd

        if use_legacy:
            tab = 'list of OTUs'
        res_df = pd.read_excel(xlsx_path, sheet_name=tab, header=None)

        return res_df[0].map(self._otu_name_parser)

    def _load_alpha_hat(self, file_path, tab='alpha', index_col='OTU name', use_legacy=False):
        import pandas as pd

        if use_legacy:
            return pd.read_csv(file_path, sep='\t')[['x']]
        else:
            df = pd.read_excel(file_path, tab, index_col=index_col)
            df.index = df.index.map(self._otu_name_parser)
            return df

    def _load_r_hat(self, file_path, tab='r', index_col='sample', use_legacy=False):
        import pandas as pd

        if use_legacy:
            return pd.read_csv(file_path, sep='\t')[['x']]
        else:
            df = pd.read_excel(file_path, tab, index_col=index_col)
            df.index = df.index.map(self._otu_name_parser)
            return pd

    @property
    def default_datasets(self):
        default_datasets = ['pre_filtered',
                            'pre_filtered_no_p36',
                            'pre_filtered_no_p36_abun_otu']
        return default_datasets

    @classmethod
    def load_default(cls, dataset=None, data_root=None):
        from os import environ
        if data_root is None:
            if 'DATA_PATH' not in environ:
                raise EnvironmentError('Please indicate the root to data folder')
            else:
                root = environ['DATA_PATH'] + 'bglmm/'
        else:
            root = data_root + 'bglmm/'

        dataset_mapper = {
            'filtered':{
                'path_to_betas': root + 'filtered/run-ALL-F-inf-corrected.xlsx',
                'path_to_otu_modeled': root + 'filtered/run-ALL-F-inf-corrected.xlsx',
                'path_to_alpha_hat': root + 'filtered/run-ALL-F-inf-corrected.xlsx',
                'path_to_r_hat': root + 'filtered/run-ALL-F-inf-corrected.xlsx',
                'name': 'filtered',
                'use_legacy': False
            },
            'filtered_no_p16': {
                'path_to_betas': root + 'filtered_no_p16/run-filtered-NEWNEW.xlsx',
                'path_to_otu_modeled': root + 'filtered_no_p16/run-filtered-NEWNEW.xlsx',
                'path_to_alpha_hat': root + 'filtered_no_p16/run-filtered-NEWNEW.xlsx',
                'path_to_r_hat': root + 'filtered_no_p16/run-filtered-NEWNEW.xlsx',
                'healed_vs_unhealed': root + 'filtered_no_p16/run-HUN-filtered-NEWNEW.xlsx',
                'name': 'filtered_no_p16',
            },
            'filtered_no_p16_abun_otu': {
                'path_to_betas': root + 'filtered_no_p16_abun_otu/run-NEWNEW.xlsx',
                'path_to_otu_modeled': root + 'filtered_no_p16_abun_otu/run-NEWNEW.xlsx',
                'path_to_alpha_hat': root + 'filtered_no_p16_abun_otu/run-NEWNEW.xlsx',
                'path_to_r_hat': root + 'filtered_no_p16_abun_otu/run-NEWNEW.xlsx',
                'healed_vs_unhealed': root + 'filtered_no_p16_abun_otu/run-HUN-NEWNEW.xlsx',
                'name': 'filtered_no_p16_abun_otu'
            }
        }

        if dataset not in dataset_mapper.keys():
            raise ValueError('dataset not found, please select from {}'.format(dataset_mapper.keys()))
        return cls(**dataset_mapper[dataset])
