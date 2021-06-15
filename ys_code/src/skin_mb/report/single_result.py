
def sig_otu_bar_plot(otu_table, model_res, sample_to_show=None):

    import pandas as pd

    sig_otus = set(model_res.pre_vs_skin.sig_otu.index.values)
    print('function reloaded')

    if sample_to_show is None:
        sample_to_show = otu_table.sample_list

    data_to_plot = pd.DataFrame(index=sample_to_show)
    data_to_plot['all_otu'] = otu_table.count_table[sample_to_show].sum(axis=0)
    data_to_plot['Significant OTUs'] = otu_table.count_table.loc[sig_otus][sample_to_show].sum(axis=0)/data_to_plot['all_otu']
    data_to_plot['Non-significant OTUs'] = 1 - data_to_plot['Significant OTUs']
    data_to_plot = data_to_plot[['Significant OTUs', 'Non-significant OTUs']]
    data_to_plot = data_to_plot.divide(data_to_plot.sum(axis=1), axis=0) * 100

    def patient_id_corrector(pid):
        return str(int(pid[:-1]) - 20) + pid[-1]

    data_to_plot.rename(index=patient_id_corrector, inplace=True)

    import matplotlib.pyplot as plt
    from ..visualizers.ThreePanelBarplot import three_panel_barplot
    three_panel_barplot(data_to_plot=data_to_plot, percent=True,
                        stacked=True, width=0.8, color=['#2C73B4', '#F39730'])
    plt.show()