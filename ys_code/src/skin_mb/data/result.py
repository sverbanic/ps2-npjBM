class Result:

    def __init__(self):
        self.pre_vs_post = None
        self.pre_vs_post_healed = None
        self.pre_vs_post_unhealed = None
        self.pre_vs_skin = None
        self.post_vs_skin = None
        self.healed_vs_unhealed_pre = None
        self.healed_vs_unhealed_post = None
        self.name = None

    def get_data(self, compare):
        return self.__dict__[compare]

    def summary(self, compare, dataset=None):
        from ..data.table import OtuTable
        import numpy as np

        if dataset is None:
            dataset = self.name
            otu_table = OtuTable.load_default(dataset=dataset)
        else:
            from pathlib import Path
            if Path(dataset).is_file():
                otu_table = OtuTable.from_biom(dataset)
            else:
                otu_table = OtuTable.load_default(dataset=dataset)
        print('Dataset {}:\n\t{}:'.format(dataset, compare))
        print('\tTotal number of OTUs: {}, number of samples: {}'.format(otu_table.shape[0], otu_table.shape[1]))
        otu_scope = otu_table.otu_stats()[otu_table.otu_stats()['rel_abun_avg_on_all'] >= 0.001].index
        print('\tNumber of OTUs with avg. rel. abun. > 0.1%: {}'.format(len(otu_scope)))
        sig_otu = self.get_data(compare).get_data()
        print('\tNumber of significant otus: {}'.format(sig_otu.shape[0]))
        print('\tNumber of positive otus:{}'.format(np.sum(sig_otu['sig'] == 1)))
        print('\tNumber of positive otus (> 0.1%):{}'.format(len([otu for otu in sig_otu[sig_otu['sig'] == 1].index
                                                                if otu in otu_scope])))
        print('\tNumber of negative otus:{}'.format(np.sum(sig_otu['sig'] == -1)))
        print('\tNumber of negative otus (> 0.1%):{}'.format(len([otu for otu in sig_otu[sig_otu['sig'] == -1].index
                                                                if otu in otu_scope])))

