def tree_to_z_matrix(path_to_tree, otu_list):

    import dendropy
    import numpy as np

    def get_z_matrix(root_node, base_num):

        def return_height(submtx):
            if isinstance(submtx, np.ndarray):
                return submtx[-1, -3]
            else:
                return 0

        if root_node.is_leaf():
            return taxon_num_map[root_node.taxon.label]

        if len(root_node.child_nodes()) != 2:
            child_a = root_node.child_nodes()[0]
            child_b = dendropy.datamodel.treemodel.Node(label='tmp_node', edge_length=0)
            for node in root_node.child_nodes()[1:]:
                child_b.add_child(node=node)
        else:
            child_a, child_b = root_node.child_nodes()

        leave_num = len(root_node.leaf_nodes())
        submtx_a = get_z_matrix(child_a, base_num)
        submtx_b = get_z_matrix(child_b, base_num + len(child_a.leaf_nodes()) - 1)
        root_height = max(return_height(submtx_a) + child_a.edge_length,
                          return_height(submtx_b) + child_b.edge_length)
        if isinstance(submtx_a, np.ndarray):
            if isinstance(submtx_b, np.ndarray):
                return np.concatenate((submtx_a,
                                       submtx_b,
                                       np.array([[submtx_a[-1, -1], submtx_b[-1, -1], root_height, leave_num,
                                                  base_num + leave_num - 2]])))
            else:
                return np.concatenate((submtx_a, np.array(
                    [[submtx_a[-1, -1], int(submtx_b), root_height, leave_num, base_num + leave_num - 2]])))
        else:
            if isinstance(submtx_b, np.ndarray):
                return np.concatenate((submtx_b, np.array(
                    [[submtx_a, submtx_b[-1, -1], root_height, leave_num, base_num + leave_num - 2]])))
            else:
                return np.array([[submtx_a, submtx_b, root_height, leave_num, base_num + leave_num - 2]])

    if isinstance(path_to_tree, str):
        tree = dendropy.Tree.get(path=path_to_tree, schema='newick')
        if otu_list is None:
            pass
        taxon_num_map = {otu: int(ix) for ix, otu in enumerate(otu_list)}
        tree.retain_taxa_with_labels(labels=otu_list)
        return get_z_matrix(tree.seed_node, base_num=len(otu_list)).astype(np.double)[:, :-1], taxon_num_map


def dendro_plot(ax, path_to_tree, otu_list):

    from scipy.cluster.hierarchy import dendrogram, set_link_color_palette

    Z, otu_to_num = tree_to_z_matrix(path_to_tree=path_to_tree, otu_list=otu_list)
    num_to_otu = {num: otu for otu, num in otu_to_num.items()}

    set_link_color_palette(['#AEAEAE'])
    dg = dendrogram(Z=Z, leaf_rotation=0, orientation='top', ax=ax, no_labels=False, color_threshold=0,
                    above_threshold_color='k')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_yaxis().set_ticks([])
    xlim = ax.get_xlim()
    xticks = ax.get_xticks()
    ax.set_xlim([xlim[0], xlim[1] + 5])
    ax.get_xaxis().set_ticks([])

    otu_reordered = [num_to_otu[num] for num in dg['leaves']]

    return otu_reordered, xlim, xticks


def sig_otu_chart(ax, pos, otu_list, chart_lim, compare_sets, label_mapper=None,
                  color=['#B2112A', '#2C73B4'], legend_label=['Positive', 'Negative']):
    import numpy as np

    for ix, compare_set in enumerate(compare_sets):
        data = compare_set['model'].load_default(compare_set['dataset'])
        sig = data.get_data(compare_set['compare']).get_data(otu_list)['sig']

        color_map = {
            -1: color[1],
            0: '#FFFFFF',
            1: color[0]
        }
        ax.scatter(x=pos, y=np.repeat(ix, len(pos)), color=sig.map(color_map).values, zorder=5)
        ax.plot(chart_lim, [ix + 0.5, ix + 0.5], color='#AEAEAE', lw=0.5)
        for p in pos:
            ax.plot([p - 5, p - 5], [-0.5, len(compare_sets) - 0.5], color='#AEAEAE', lw=0.5)
        ax.plot([p + 5, p + 5], [-0.5, len(compare_sets) - 0.5], color='#AEAEAE', lw=0.5)

    ax.plot(chart_lim, [-0.5, -0.5], color='#AEAEAE', lw=1)
    ax.set_xlim([chart_lim[0], chart_lim[1] + 5])
    ax.set_xticks(pos)
    ax.set_ylim(-0.5, len(compare_sets) - 0.49)
    ax.set_yticks(np.linspace(0, len(compare_sets) - 1, len(compare_sets)))
    ax.set_yticklabels([compare_set['label'] for compare_set in compare_sets], fontsize=12)
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis=u'both', which=u'both', length=0)

    if label_mapper is None:
        ax.set_xticklabels(otu_list, rotation=270, ha='center')
    else:
        ax.set_xticklabels([label_mapper[otu] for otu in otu_list], rotation=270, ha='center')

    lgd = [ax.scatter([], [], color=color[0], label=legend_label[0]),
           ax.scatter([], [], color=color[1], label=legend_label[1])]
    ax.legend(handles=lgd, labels=legend_label, bbox_to_anchor=(1, 0.5), loc='left', frameon=False)