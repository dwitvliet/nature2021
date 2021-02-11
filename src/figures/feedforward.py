from collections import defaultdict

from src.data import data_manager
from src.data.neuron_info import ntype
from src.data.dataset_info import all_datasets, datasets_with_adj, timepoint

from src.plotting import plotter


class Figure(object):

    def __init__(self, output_path, page_size=7.20472):
        self.plt = plotter.Plotter(output_path=output_path, page_size=page_size)


    def _feed_type(self, edge, only_sensory=False, with_muscle=False):
        edge_type = (ntype(edge[0]), ntype(edge[1]))
        if only_sensory:
            if edge_type in (('sensory', 'inter'), ('sensory', 'motor'), ('sensory', 'modulatory')):
                return 'Feed-forward'
            if edge_type in (('inter', 'sensory'), ('motor', 'sensory'), ('modulatory', 'sensory')):
                return 'Feed-back'
            if edge_type in (('sensory', 'sensory'), ):
                return 'Recurrent'
            return None
        if edge_type in (('sensory', 'inter'), ('inter', 'motor'), ('sensory', 'motor'), ('modulatory', 'inter'), ('sensory', 'modulatory'), ('modulatory', 'motor')):
            return 'Feed-forward'
        if edge_type in (('inter', 'sensory'), ('motor', 'inter'), ('motor', 'sensory'), ('inter', 'modulatory'), ('modulatory', 'sensory')):
            return 'Feed-back'
        if edge_type in (('sensory', 'sensory'), ('inter', 'inter'), ('motor', 'motor'), ('modulatory', 'modulatory')):
            return 'Recurrent'
        if with_muscle and edge_type[1] == 'muscle':
            return 'Feed-forward'
        return None

    def feedforward_stable_increase(self, f, edge_classifications, use_size=False):

        G = data_manager.get_connections()['size' if use_size else 'count'].copy()
        G = G[G.sum(axis=1) > 0]  # remove edges without size if need be
        G = data_manager.remove_postemb(G)

        edge_classifications = edge_classifications.copy()

        edges = [e for e in G.index if edge_classifications[e] == 'stable']

        feed_types = ['Feed-forward', 'Recurrent', 'Feed-back']
        feed_colors = {'Feed-forward': '#C7EAE4', 'Feed-back': '#EAC9C1', 'Recurrent': 'white'}
        syn_increases = {ft: [] for ft in feed_types}
        
        for edge in edges:

            feed_type = self._feed_type(edge)
            if not feed_type:
                continue

            syns = G.loc[edge]
            if syns[0] == 0 and syns[1] == 0:
                continue

            syn_increase_relative = syns[['Dataset7', 'Dataset8']].mean() / syns[['Dataset1', 'Dataset2']].mean()
            syn_increases[feed_type].append(syn_increase_relative)

        data, c, l = [], [], []
        for ft in feed_types:
            data.append(syn_increases[ft])
            l.append(ft)
            c.append(feed_colors[ft])
            
        if use_size:
            ylim = (0, 15)
            yticks = (0, 5, 10, 15)
            y_label = 'Relative synapse volume increase'
            size = (0.15, 0.15)
        else:
            ylim = (0, 12)
            yticks = (0, 4, 8, 12)
            y_label = 'Relative synapse addition'
            size = (0.15, 0.15)


        self.plt.plot(
            'box_plot', data, size=size,
            margin={'left': 0.04, 'right': 0.01, 'top': 0.05, 'bottom': 0.04},
            colors=c, xticklabels=l, xtickpad=3, xpad=5, ylim=ylim, yticks=yticks,
            y_label=y_label, x_label='Stable connection directionality',
            show_outliers=False, stats=((2, 3), (1, 2), (1, 3)),
            save=f+'_feedforward_stable_increase' + ('_size' if use_size else '')
        )


    def feedforward_edge_proportion(self, f, edge_classifications, use_size=False):

        G = data_manager.get_connections()['size' if use_size else 'count'].copy()
        G = G[G.sum(axis=1) > 0] # remove edges without size if need be
        G = data_manager.remove_postemb(G)

        edge_classifications = edge_classifications.copy()

        feed_types = ['Feed-back', 'Recurrent', 'Feed-forward']
        edge_types = ('stable', 'increase', 'decrease')
        feed_colors = {'Feed-forward': '#C7EAE4', 'Feed-back': '#EAC9C1', 'Recurrent': 'white'}

        connections = defaultdict(lambda: {ft: 0 for ft in feed_types})

        edges_per_type = {
            'stable': [e for e in G.index if edge_classifications[e] == 'stable'],
            'increase': [e for e in G.index if edge_classifications[e] == 'increase'],
            'decrease': [e for e in G.index if edge_classifications[e] == 'decrease'],
#            'Variable': [e for e in G if edge_classifications[(npair(e[0]), npair(e[1]))] in ('remainder', 'noise')]
        }

        xlabels = (
            'Stable',
            'Strengthened',
            'Weakened',
        )

        for edge_type, edges in edges_per_type.items():
            for edge in edges:
                feed_type = self._feed_type(edge)
                if not feed_type:
                    continue
                connections[edge_type][feed_type] += 1



        data = tuple((ft, [connections[et][ft] for et in edge_types]) for ft in feed_types)
        print(data)

        self.plt.plot(
            'stacked_bar_graph', data, stats=((1, 2), (1, 3)), size=(0.15, 0.15),
            margin={'left': 0.04, 'right': 0.08, 'top': 0.05, 'bottom': 0.04},
            y_label='Proportion of connections',
            colors=feed_colors, xlabels=xlabels, x_label='Connection classification',
            xtickpad=3, xpad=5,
            legendpos='right', legendcol=1, legendreverse=True, width=0.5,
            save=f+'_feedforward_edge_proportion'
        )


    def feedforward_global_shift(self, f, only_sensory=False, use_size=False):
        
        datasets = list(datasets_with_adj if use_size else all_datasets)

        G = data_manager.get_connections()['size' if use_size else 'count'].copy()
        G = G[G.sum(axis=1) > 0] # remove edges without size if need be
        G = G[datasets]
        G = data_manager.remove_postemb(G)

        y_label = 'Proportion of synapses'
        ylim = (0, 0.6)
        
        if use_size:
            y_label = 'Proportion of synapse volume'
            

        if only_sensory:
            y_label += ' to\nor from sensory neurons'
            ylim = (0, 0.8)

        feed_types = ['Feed-back', 'Feed-forward', 'Recurrent']

        G['feed_types'] = G.index.map(self._feed_type)
        feed_type_counts = G.groupby('feed_types').sum()
        feed_type_counts = feed_type_counts / feed_type_counts.sum()

        xs = [timepoint[d] for d in datasets]
        colors = {
            'Feed-forward': '#C7EAE4', 'Feed-back': '#EAC9C1', 'Recurrent': 'white',
            'Feed-forward_edge': '#7ccfc1', 'Feed-back_edge': '#d18876', 'Recurrent_edge': '#999999'
        }
        
        data = (xs, [(ft, feed_type_counts.loc[ft]) for ft in feed_types])

        self.plt.plot(
            'xy_graph', data, size=(0.10, 0.19),
            margin={'left': 0.04, 'right': 0.10, 'top': 0.01, 'bottom': 0.04},
            y_label=y_label, ylim=ylim, stats='spearmanr', colors=colors,
            x_label='Developmental age',
            legendpos='right', rev_legend=True, legend_shift_top=0.03, legend_shift_right=0.05,
            save=f+'_feedforward_global_shift', linkpoints=False,
#            hlines=(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
        )

