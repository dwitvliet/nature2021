import os

import numpy as np

from src.data import data_manager
from src.data.dataset_info import all_datasets, timepoint

from src.plotting import plotter


categories = {
    'increase': 'Developmentally dynamic',
    'decrease': 'Developmentally dynamic',
    'stable': 'Stable',
    'noise': 'Variable',
    'remainder': 'Variable'
}


class Figure(object):

    def __init__(self, output_path, page_size=7.20472):
        self.path = output_path
        self.plt = plotter.Plotter(output_path=output_path, page_size=page_size)

    def list_connections_between_neurons(self, fname, edge_classifications,
                                         dataset='Dataset8'):

        G = data_manager.get_connections()['count'].copy()
        G = data_manager.remove_postemb(G)
        G_pair = data_manager.to_npair(G)

        ns = 'AWA, AWC, AWB, ASH, ADL, AIZ, AIA, AVA, AVD, AIB, AIY, AFD, RIF, AVH, AVB'
        ns = ns.replace(' ', '').split(',')
        with open(os.path.join(self.path, f'{fname}_edges_{dataset}.txt'), 'w') as f:

            for (pre, post), c in sorted(edge_classifications.iteritems()):
                e = (pre, post)

                if pre not in ns or post not in ns:
                    continue

                s = np.mean(G_pair.loc[e][dataset])
                if s == 0:
                    continue

                s_normalized = min(s ** (1 / 1.08), 2 * 8)

                f.write(' '.join([pre, post, str(s_normalized / 8), c,
                                  str(G_pair.loc[e].to_numpy())]) + '\n')

    def edge_classification_count(self, f, edge_classifications,
                                  edge_pair_classifications, to_sum='connections',
                                  datasets=all_datasets, paired=False, small=True):

        datasets = list(datasets)

        colors = {
            'Variable': '#eeeeee', 'Stable': '#aaaaaa',
            'Developmentally dynamic': '#6DB6FF',
            'Variable_edge': '#aaaaaa', 'Stable_edge': '#4C4D4C',
            'Developmentally dynamic_edge': '#006DDB'
        }

        G = data_manager.get_connections()['size' if to_sum == 'size' else 'count'][datasets].copy()
        G = data_manager.remove_postemb(G)
        classifications = edge_classifications

        if paired:
            G = data_manager.to_npair(G)
            classifications = edge_pair_classifications

        classifications = classifications.map(categories)
        G_classifications = G.merge(classifications, left_index=True,
                                    right_index=True)

        if to_sum == 'connections':
            count = G_classifications.groupby('edge_classifications').agg(
                lambda s: s.astype(bool).sum())

        else:
            count = G_classifications.groupby('edge_classifications').sum()

        x = [timepoint[d] for d in datasets]
        y = tuple(sorted(count.iterrows()))

        y_label = 'Connections (#)'
        name = 'connections'
        ylim = 2000
        yticks = range(0, 2500, 500)
        if paired:
            y_label = 'Pair connections (#)'
            name = 'pair_connections'
            ylim = 1200
            yticks = range(0, 1500, 400)
        if to_sum == 'synapses':
            y_label = 'Synapses (#)'
            name = 'synapses'
            ylim = 8000
            yticks = range(0, 10000, 2000)

        size = 0.41
        markersize = 7
        if small:
            size = 0.275
            markersize = 5

        self.plt.plot(
            'xy_graph', (x, y), y_label=y_label, x_label='Developmental age',
            cumulative=True, rev_legend=True, ylim=(0, ylim),
            save=f + '_edge_classification_counts_' + name,
            yticks=yticks, clipon=True,
            size=size,
            margin={'left': 0.07, 'right': 0.01, 'top': 0.03, 'bottom': 0.05},
            colors=colors, legend_shift_right=0.13, legend_shift_top=-0.02,
            legendcol=3, markersize=markersize
        )

        # print(count)
        # print(count / count.sum())
