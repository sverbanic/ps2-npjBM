#!/usr/bin/env python3

import numpy as np
import pandas as pd

class OtuTable(object):
    """A wrapper over ``biom.Table`` with extra functions to work with otu table"""

    def __init__(self, table, meta_otu=None, meta_sample=None, file_path=None, sort_sample_by='name', name=None, note=None):
        import pandas as pd

        self._table = table
        self._sample_list = pd.Series(data=self._table.ids(axis='sample'), dtype=str)
        self._otu_list = pd.Series(data=self._table.ids(axis='observation'), dtype=str)
        if sort_sample_by is None or sort_sample_by is False:
            pass
        else:
            self.sort_sample_id(sort_by=sort_sample_by)
        if file_path:
            self.source = file_path
        self.name = name
        self.note = note
        self.vis = OtuTableVis(self)

    def __repr__(self):
        return 'OTU Table (<class skin_mb.data.OtuTable>): current sample: {}; otu: {}'.format(len(self._sample_list),
                                                                                               len(self._otu_list))

    @property
    def otu_list(self):
        return self._otu_list

    @otu_list.setter
    def otu_list(self, value):
        self._otu_list = value

    @property
    def sample_list(self):
        return self._sample_list

    @sample_list.setter
    def sample_list(self, value):
        self._sample_list = value

    @property
    def shape(self):
        return len(self._otu_list), len(self._sample_list)

    @property
    def count_table(self):
        table = self._table.filter(ids_to_keep=self._sample_list, axis='sample').filter(ids_to_keep=self._otu_list,
                                                                                        axis='observation')
        table = table.to_dataframe(dense=True).loc[self._otu_list][self._sample_list]
        table.index.name = 'OTU'
        return table

    def metadata_to_dataframe(self, ids=None, axis='otu'):
        if axis.lower() in ['otu', 'observation', 'feature']:
            if ids is None:
                ids = self._otu_list
            metadata = {id_: dict(self._table.metadata(id=id_, axis='observation')) for id_ in ids}
        elif axis.lower() in ['sample']:
            if ids is None:
                ids = self._sample_list
            metadata = {id_: dict(self._table.metadata(id=id_, axis='sample')) for id_ in ids}
        else:
            raise ValueError('`axis` should be either otu or sample')

        if metadata != {}:
            import pandas as pd
            import numpy as np

            def flatten(d):
                if not d:
                    return {'dummy': np.nan}
                d_flatten = {}
                for key in d.keys():
                    if isinstance(d[key], list) or isinstance(d[key], tuple) or isinstance(d[key], np.ndarray):
                        for ix,value in enumerate(d[key]):
                            d_flatten[str(key) + '_' + str(ix)] = value
                    else:
                        d_flatten[str(key)] = d[key]
                return d_flatten

            metadata = {id_:flatten(metadata[id_]) for id_ in metadata.keys()}
            df = pd.DataFrame.from_dict(data=metadata, orient='index')
            df.index.name = 'OTU'
            columns = sorted([col for col in df.columns if col != 'dummy'])
            return df[columns]
        else:
            raise Warning('no metadata for {}'.format(axis))

    def reset(self):
        self._sample_list = pd.Series(data=self._table.ids(axis='sample'), dtype=str)
        self._otu_list = pd.Series(data=self._table.ids(axis='observation'), dtype=str)

    def purge(self):
        self._table = self._table.filter(ids_to_keep=self._sample_list, axis='sample').filter(ids_to_keep=self._otu_list,
                                                                                              axis='observation')

    def get_rel_abun(self, on_original=False):
        if on_original:
            rel_abun = self._table.to_dataframe(dense=True)
        else:
            rel_abun = self.count_table
        return rel_abun / rel_abun.sum(axis=0)

    def sort_sample_id(self, sort_by=None, reverse=False):
        if sort_by is None:
            raise ValueError('`sort_by` can not be None')
        elif sort_by == 'name':
            self._sample_list = sorted(self._sample_list, reverse=reverse)
        elif callable(sort_by):
            self._sample_list = sorted(self._sample_list, key=sort_by, reverse=reverse)
        else:
            raise ValueError('`sort_by` should either be name or a callable')

    def get_data(self, ids, axis='otu', rel_abun=False, on_original=False):
        if isinstance(ids, str):
            ids = [ids]
        if axis.lower() in ['otu', 'observation', 'feature']:
            if rel_abun:
                return self.get_rel_abun(on_original=on_original).loc[ids]
            else:
                return self.count_table.loc[ids]
        elif axis.lower() in ['sample']:
            if rel_abun:
                return self.get_rel_abun(on_original=on_original)[ids]
            else:
                return self.count_table[ids]
        else:
            raise ValueError('`axis` should be either otu or sample')

    def filter(self, ids_to_keep=None, ids_to_remove=None, axis='otu', inplace=False, remove_empty=True):
        import pandas as pd
        if inplace:
            handle = self
        else:
            import copy
            handle = copy.deepcopy(self)
        if ids_to_keep is None:
            if ids_to_remove is None:
                raise ValueError('`ids_to_keep` and `ids_to_remove` are both None')
            else:
                if axis.lower() in ['otu', 'observation', 'feature']:
                    ids_to_keep = [otu for otu in handle._otu_list if otu not in ids_to_remove]
                elif axis.lower() in ['sample']:
                    ids_to_keep = [sample for sample in handle._sample_list if sample not in ids_to_remove]

        if axis.lower() in ['otu', 'observation', 'feature']:
            for otu in ids_to_keep:
                if otu not in handle._otu_list:
                    raise ValueError('{} is not in current table'.format(otu))
            handle._otu_list = pd.Series(data=ids_to_keep)
        elif axis.lower() in ['sample']:
            for sample in ids_to_keep:
                if sample not in handle._sample_list:
                    print(sample)
                    print(handle._sample_list)
                    raise ValueError('{} is not in current table'.format(sample))
            handle._sample_list = pd.Series(data=ids_to_keep)
        else:
            raise ValueError('`axis` should be either otu or sample')
        if remove_empty:
            handle.remove_empty(axis='both', inplace=True)
        if not inplace:
            handle.purge()
            return handle

    def remove_empty(self, axis='both', inplace=True):

        if inplace:
            handle = self
        else:
            import copy
            handle = copy.deepcopy(self)
        if axis.lower() in ['otu', 'observation', 'feature', 'both']:
            handle._otu_list = handle.count_table.index[handle.count_table.sum(axis=1) > 0]
        if axis.lower() in ['sample', 'both']:
            handle._sample_list = handle.count_table.columns[handle.count_table.sum(axis=0) > 0]
        if not inplace:
            handle.purge()
            return handle

    def filter_objects(self, obj_to_keep=None, obj_to_remove=None, inplace=True, remove_empty=True):
        if obj_to_keep is None:
            if obj_to_remove is None:
                raise ValueError('Please indicate either `obj_to_keep` or `obj_to_remove`')
            else:
                samples = []
                for obj in obj_to_remove:
                    samples += [str(obj) + 'A', str(obj) + 'B', str(obj) + 'C']
                remove = True
        else:
            samples = []
            for obj in obj_to_keep:
                samples += [str(obj) + 'A', str(obj) + 'B', str(obj) + 'C']
            remove = False
        if inplace:
            if remove:
                self.filter(ids_to_remove=samples, axis='sample', inplace=inplace, remove_empty=remove_empty)
            else:
                self.filter(ids_to_keep=samples, axis='sample', inplace=inplace, remove_empty=remove_empty)
        else:
            if remove:
                return self.filter(ids_to_remove=samples, axis='sample', inplace=inplace, remove_empty=remove_empty)
            else:
                return self.filter(ids_to_keep=samples, axis='sample', inplace=inplace, remove_empty=remove_empty)

    def otu_stats(self, on_original=False):
        import pandas as pd

        if on_original:
            count_table = self._table.to_dataframe(dense=True)
        else:
            count_table = self.count_table
        rel_abun = self.get_rel_abun()
        df = pd.DataFrame(index=count_table.index)
        df['total_counts'] = count_table.sum(axis=1)
        df['smpl_detected'] = (count_table > 0).sum(axis=1)
        df['counts_avg_on_all'] = df['total_counts']/count_table.shape[1]
        df['counts_avg_on_detected'] = df['total_counts']/df['smpl_detected']
        df['rel_abun_avg_on_all'] = rel_abun.sum(axis=1)/rel_abun.shape[1]
        df['avg_abun_on_detected'] = rel_abun.sum(axis=1)/df['smpl_detected']
        return df

    def get_count_pct(self, otu_list):
        counts = self.otu_stats()['total_counts']
        return counts[otu_list].sum()/counts.sum()

    @classmethod
    def from_biom(cls, file_path, note=None):
        import biom
        table = biom.load_table(file_path)
        return cls(table=table, file_path=file_path, note=note)

    def to_biom(self, file_path, with_meta=True, save_original=False, note=None):
        """Save as `.biom` file that could be read by `Biom`"""
        if save_original:
            table = self._table
        else:
            table = self._table.filter(ids_to_keep=self._sample_list, axis='sample').filter(ids_to_keep=self._otu_list,
                                                                                            axis='observation')
        if note is None:
            if self.note is not None:
                note = self.note
            else:
                note = 'biom table processed by skim_mb package'
        if not with_meta:
            table.del_metadata()
        import h5py
        with h5py.File(file_path, 'w') as handle:
            table.to_hdf5(handle, generated_by=note)
        print('OTU table (sample: {}, OTU: {}) saved to {}'.format(table.shape[1], table.shape[0], file_path))

    @classmethod
    def from_csv(cls, file_path, sep=',', index=0, note=None, **kwargs):
        from biom.table import Table
        import pandas as pd
        data = pd.read_csv(file_path, sep=sep, index_col=index)
        return cls(table=Table(data=data.values, observation_ids=data.index, sample_ids=data.columns),
                   file_path=file_path, note=note)

    def to_csv(self, file_path, include_meta=True, no_header=False, save_original=False, sep=',', **kwargs):
        """Save as csv file(s)"""
        if save_original:
            table = self._table.to_dataframe(dense=True)
        else:
            table = self.count_table
        from pathlib import Path

        file_path = Path(file_path)
        if no_header:
            table.to_csv(str(file_path.parent/file_path.stem) + '_counts.csv', header=False, index=False, sep=sep)
            table.index.to_csv(str(file_path.parent/file_path.stem) + '_otus.csv', index=False)
            table.column.to_csv(str(file_path.parent/file_path.stem) + '_samples.csv', index=False, header='sample')
        else:
            table.to_csv(str(file_path.parent/file_path.stem) + '.csv', sep=sep)

        print('OTU table (sample: {}, OTU: {}) saved to {}'.format(table.shape[1], table.shape[0], file_path))
        if include_meta:
            if self._table.metadata(axis='observation') is not None:
                self.metadata_to_dataframe(axis='otu').to_csv(str(file_path.parent/file_path.stem) + '_otu_meta.csv')
            if self._table.metadata(axis='sample') is not None:
                self.metadata_to_dataframe(axis='sample').to_csv(str(file_path.parent/file_path.stem) + '_sample_meta.csv')

    @classmethod
    def read_table(cls, file_path, note=None):
        from pathlib import Path
        if Path(file_path).suffix in ['.biom', '.bm']:
            return cls.from_biom(file_path=file_path, note=note)
        elif Path(file_path).suffix in ['.tsv', '.txt', '.csv']:
            return cls.from_csv(file_path=file_path, note=note)

    @classmethod
    def load_default(cls, dataset='full_table', data_root=None):
        from os import environ
        if data_root is None:
            if 'DATA_PATH' not in environ:
                raise EnvironmentError('Please indicate the root to data folder')
            else:
                root = environ['DATA_PATH']
        else:
            root = data_root
        table_path = root + '/otu_tables/' + dataset + '.biom'
        return cls.from_biom(table_path, note='BIOM table loaded from {}'.format(table_path))


def taxonomy_corrector(otu_meta):
    col_mapper = {
        'taxonomy_0': 'Domain',
        'taxonomy_1': 'Phylum',
        'taxonomy_2': 'Class',
        'taxonomy_3': 'Order',
        'taxonomy_4': 'Family',
        'taxonomy_5': 'Genus',
        'taxonomy_6': 'Species'
    }

    def taxo_name_corrector(name):
        import numpy as np
        if isinstance(name, float):
            return np.nan
        elif isinstance(name, str):
            if 'unknown' in name.lower() or 'ambiguous' in name.lower() or 'uncultured' in name.lower() or 'unassign' in name.lower():
                return np.nan
            else:
                return name.split('__')[-1]
        else:
            return np.nan

    return otu_meta.applymap(func=taxo_name_corrector).rename(columns=col_mapper)


def get_otu_id_to_taxo_name(otu_table=None, use_builtin='lowest_taxa', formatter=None):
    from os import environ
    if 'DATA_PATH' not in environ:
        root = '/Users/randal/research/projects/skin_wound_microbiome/data/'
    else:
        root = environ['DATA_PATH']

    if otu_table is None:
        otu_table = 'full_table'

    if isinstance(otu_table, str):
        otu_table = OtuTable.from_biom(root + 'otu_tables/' + otu_table + '.biom')
    taxonomy = taxonomy_corrector(otu_table.metadata_to_dataframe(axis='otu'))

    if formatter is None:
        if use_builtin == 'genus_species':
            def genus_species(row):
                taxonomy_rank = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
                for domain in taxonomy_rank[::-1]:
                    if domain in row.index:
                        if isinstance(row[domain], str):
                            if domain == 'Species':
                                return '{} {}'.format(row['Genus'], row['Species'].lower())
                            elif domain == 'Genus':
                                return row[domain].title()
                            else:
                                return row[domain].title() + '({})'.format(domain[0])
                return row.name
            formatter = genus_species
        elif use_builtin=='lowest_taxa':
            def lowest_taxa(row):
                taxonomy_rank = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
                for domain in taxonomy_rank[::-1]:
                    if domain in row.index:
                        if isinstance(row[domain], str):
                            if domain == 'Genus':
                                return row[domain]
                            else:
                                return row[domain] + '({})'.format(domain[0])

                return row.name
            formatter = lowest_taxa
    return taxonomy.apply(formatter, axis=1)


class OtuTableVis:
    
    def __init__(self, otu_table):
        self.otu_table = otu_table
    
    def alpha_diversity(self, sample_list=None, method='shannon', plot_off=False, ax=None, pos_map=None, **kwargs):
        from skbio.diversity import alpha
        import numpy as np

        def simpson_inv(counts):
            return 1/(1-alpha.simpson(counts))

        method_mapper = {
            'shannon': (alpha.shannon, 'bit'),
            'chao1': (alpha.chao1, None),
            # 'equitability': (alpha.equitability, None),
            # 'entropy_efficienty': (alpha.equitability, None),
            'simpson': (alpha.simpson, None),
            'simpson_inv': (simpson_inv, None)
        }

        if sample_list is None:
            sample_list = self.otu_table.sample_list

        alphas = pd.Series(data=self.otu_table.count_table[sample_list].apply(method_mapper[method][0], axis=0),
                           index=sample_list)
        if not plot_off:
            if ax is None:
                import matplotlib.pyplot as plt
                fig = plt.figure(figsize=(len(sample_list)/4, 6))
                ax = fig.add_subplot(111)
            if pos_map is None:
                pos = np.linspace(0, len(sample_list) - 1, len(sample_list))
            else:
                pos = [pos_map[sample] for sample in sample_list]
            ax.scatter(pos, alphas, **kwargs)
            return ax, alphas
        else:
            return alphas


class OtuTableSet:

    def __init__(self, table_list):
        if isinstance(table_list, list):
            if isinstance(table_list[0], OtuTable):
                table_list = {ix:table for ix,table in enumerate(table_list)}
            else:
                raise TypeError('Input should be a list/dict of `OtuTable` objects')
        elif isinstance(table_list, dict):
            if isinstance(list(table_list.values())[0], OtuTable):
                pass
            else:
                raise TypeError('Input should be a list/dict of `OtuTable` objects')
        self.__dict__.update(table_list)

    @classmethod
    def from_records(cls, paths_to_files, alias=None, all_file_with_pattern=None):
        from pathlib import Path
        if all_file_with_pattern is not None:
            if isinstance(paths_to_files, str):
                paths_to_files = [paths_to_files]
            file_paths = []
            for path in paths_to_files:
                if Path(path).is_dir():
                    root = Path(path)
                else:
                    root = Path(path).parent
                file_paths += [file.absolute() for file in root.glob(all_file_with_pattern)]

            file_paths = list(set(file_paths))
            return cls(table_list={Path(file).stem: OtuTable.read_table(file) for file in file_paths})
        else:
            if isinstance(paths_to_files, str):
                paths_to_files = [paths_to_files]
            if alias is not None:
                if isinstance(alias, str):
                    alias = [alias]
                if len(alias) != len(paths_to_files):
                    raise ValueError('Lengths of `paths_to_files` and `alias` do not match')
                else:
                    return cls(table_list={alia:OtuTable.read_table(path) for path,alia in zip(paths_to_files, alias)})
            else:
                return cls(table_list={Path(path).stem: OtuTable.read_table((path)) for path in paths_to_files})

    def shared_otu(self, datasets=None, otu_list=None, ax=None):
        if datasets is None:
            datasets = self.__dict__.keys()
        values = {dataset: self.__dict__[dataset] for dataset in datasets}
        common_otu = 'init'
        for dataset in datasets:
            if common_otu == 'init':
                common_otu = set(values[dataset].otu_list)
            else:
                common_otu = common_otu & set(values[dataset].otu_list)

        if otu_list:
            common_otu = common_otu & set(otu_list)
        common_otu = list(common_otu)

        pct = [values[dataset].get_count_pct(common_otu) for dataset in datasets]
        if ax is None:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(len(datasets), 4))
            ax = fig.add_subplot(111)
        pos = np.linspace(0, len(datasets) - 1, len(datasets))
        ax.bar(pos, height=pct, width=0.8, color='#2C73B4')
        ax.set_ylabel('Percent of reads in common OTUs', fontsize=12)
        ax.set_xticks(pos)
        ax.set_xticklabels(datasets, fontsize=12, rotation=60, ha='right')


    def compare_alpha(self, datasets=None, samples=None, method='shannon', ax=None, pos_map=None, label_map=None, **kwargs):
        if datasets is None:
            datasets = self.__dict__.keys()

        values = {dataset:self.__dict__[dataset] for dataset in datasets}

        if samples is None:
            samples = []
            for dataset in datasets:
                samples += list(values[dataset].sample_list)
            samples = sorted(list(set(samples)))

        if pos_map is None:
            pos_map = {sample:pos for pos,sample in enumerate(samples)}

        if label_map is None:
            label_map = {dataset:dataset for dataset in datasets}

        if ax is None:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(len(samples)/4, 6))
            ax = fig.add_subplot(111)

        for dataset in datasets:
            _ = values[dataset].vis.alpha_diversity(ax=ax, method=method, pos_map=pos_map, label=label_map[dataset], **kwargs)
        ax.set_ylabel(method, fontsize=14)
        ax.set_xticks([pos for pos in pos_map.values()])
        ax.set_xticklabels([label for label in pos_map.keys()], fontsize=12, rotation=90)
        ax.legend(frameon=False)


def _otu_info_to_csv(otu_table, export_dir=None, include_meta=True, include_raw_counts=False):
    otu_df = OtuTable.survey_otu_properties(otu_table)
    if include_meta:
        otu_df = pd.concat([otu_df, otu_table.metadata_otu()], axis=1, join='inner')
    if include_raw_counts:
        otu_df = pd.concat([otu_df, otu_table.table.to_dataframe(dense=True)], axis=1, join='inner')
    if export_dir:
        otu_df.to_csv(export_dir)
    else:
        return otu_df


def _otu_table_filter(otu_table, filters, sample_select_first=False, inplace=False):

    import copy

    otu_table_filtered = copy.deepcopy(otu_table)
    if sample_select_first:
        if filters.samples_to_keep is not None:
            print('Filter samples: table shape {} --->'.format(otu_table_filtered.shape), end='')
            otu_table_filtered.filter_samples(samples_to_keep=filters.samples_to_keep,
                                              exclude_zeros=True, inplace=True)
            print('{}'.format(otu_table_filtered.shape))
        if filters.objects_to_keep is not None:
            print('Filter objects: table shape {} --->'.format(otu_table_filtered.shape), end='')
            otu_table_filtered.filter_objects(objects_to_keep=filters.objects_to_keep,
                                              exclude_zeros=True, inplace=True)
            print('{}'.format(otu_table_filtered.shape))

    if filters.exclude_controls == True:
        print('Exclude control samples: table shape {} --->'.format(otu_table_filtered.shape), end='')
        samples_to_keep = [
            sample for sample in otu_table_filtered.sample_list
            if ('CL' not in sample.upper())and('MC' not in sample.upper())and('WC' not in sample.upper())
        ]
        otu_table_filtered.filter_samples(samples_to_keep=samples_to_keep, exclude_zeros=True, inplace=True)
        print('{}'.format(otu_table_filtered.shape))

    if filters.count_cutoff is not None:
        print('Apply count cutoff: table shape {} --->'.format(otu_table_filtered.shape), end='')
        otu_table_filtered.cutoff_count(filters.count_cutoff, inplace=True, exclude_zeros=True)
        print('{}'.format(otu_table_filtered.shape))
    if filters.abun_cutoff is not None:
        print('Apply relative abun cutoff: table shape {} --->'.format(otu_table_filtered.shape), end='')
        otu_table_filtered.cutoff_abun(filters.abun_cutoff, inplace=True, exclude_zeros=True)
        print('{}'.format(otu_table_filtered.shape))

    otu_df = otu_table_filtered.otu_stats()

    print('Apply OTU filters... starting shape: {}'.format(otu_table_filtered.shape))
    if filters.min_sample_num is not None:
        print('\tApply minimal sample num: OTU {} --->'.format(len(otu_df)), end='')
        otu_df = otu_df[otu_df['num_detected'] >= filters.min_sample_num]
        print('{}'.format(len(otu_df)))

    if filters.avg_count_on_all is not None:
        print('\tApply average counts on all: OTU {} --->'.format(len(otu_df)), end='')
        otu_df = otu_df[otu_df['avg_counts_on_all'] >= filters.avg_count_on_all]
        print('{}'.format(len(otu_df)))

    if filters.avg_count_on_detected is not None:
        print('\tApply average counts on detected: OTU {} --->'.format(len(otu_df)), end='')
        otu_df = otu_df[otu_df['avg_counts_on_detected'] >= filters.avg_count_on_detected]
        print('{}'.format(len(otu_df)))

    if filters.avg_abun_on_all is not None:
        print('\tApply average relative abun counts on all: OTU {} --->'.format(len(otu_df)), end='')
        otu_df = otu_df[otu_df['avg_abun_on_all'] >= filters.avg_abun_on_all]
        print('{}'.format(len(otu_df)))

    if filters.avg_abun_on_detected is not None:
        print('\tApply average relative abun on detected: OTU {} --->'.format(len(otu_df)), end='')
        otu_df = otu_df[otu_df['avg_abun_on_detected'] >= filters.avg_abun_on_detected]
        print('{}'.format(len(otu_df)))

    otu_table_filtered.filter_otus(otus_to_keep=otu_df.index, inplace=True)
    print('end shape: {}'.format(otu_table_filtered.shape))

    if not sample_select_first:
        if filters.samples_to_keep is not None:
            print('Filter samples: table shape {} --->'.format(otu_table_filtered.shape), end='')
            otu_table_filtered.filter_samples(samples_to_keep=filters.samples_to_keep,
                                              exclude_zeros=True, inplace=True)
            print('{}'.format(otu_table_filtered.shape))
        if filters.objects_to_keep is not None:
            print('Filter objects: table shape {} --->'.format(otu_table_filtered.shape), end='')
            otu_table_filtered.filter_objects(objects_to_keep=filters.objects_to_keep,
                                              exclude_zeros=True, inplace=True)
            print('{}'.format(otu_table_filtered.shape))

    if inplace:
        otu_table.table = otu_table_filtered.table
    else:
        return otu_table_filtered


class OtuTableFilter(object):

    def __init__(self):
        self.exclude_controls = None
        self.count_cutoff = None
        self.abun_cutoff = None
        self.min_sample_num = None
        self.avg_count_on_all = None
        self.avg_count_on_detected = None
        self.avg_abun_on_all = None
        self.avg_abun_on_detected = None

        self.samples_to_keep = None
        self.objects_to_keep = None

    def print_filters(self):
        print('Filters added:')
        for filter_name, value in self.__dict__.items():
            if value is not None:
                print('\t{}: {}'.format(filter_name, value))
