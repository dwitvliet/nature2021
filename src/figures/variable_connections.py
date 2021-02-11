from collections import defaultdict

import numpy as np
import pandas as pd

from src.data import data_manager
from src.data.dataset_info import datasets_with_adj
from src.data.neuron_info import ntype, dark_colors

from src.plotting import plotter

categories = {
    'increase': 'Developmental change',
    'decrease': 'Developmental change',
    'stable': 'Stable',
    'noise': 'Variable',
    'remainder': 'Variable'
}


class Figure(object):

    def __init__(self, output_path, page_size=7.20472):
        self.plt = plotter.Plotter(output_path=output_path, page_size=page_size)

    def variable_edge_prevalence_quantified(self, f, edge_classifications,
                                            inputs=False, variable=None,
                                            nonvariable=None):

        G = data_manager.get_connections()['count'].copy()
        G = data_manager.remove_postemb(G)
        
        classifications = edge_classifications.copy()

        variable_defined_by_adult = bool(variable or nonvariable)

        if not variable_defined_by_adult:
            nonvariable = defaultdict(float)
            variable = defaultdict(float)

            for (pre, post), syns in G.iterrows():
                if np.max(syns[-2:]) <= 0:
                    continue

                s = np.count_nonzero(syns[-2:])
                if s == 0:
                    continue

                classification = classifications[(pre, post)]
                if classification in ('stable', 'increase', 'decrease'):
                    nonvariable[post if inputs else pre] += s

                if classification in ('remainder', 'noise'):
                    variable[post if inputs else pre] += s



        y_label='Proportion of variable connection'

        if inputs:
            y_label += ' inputs'
            x_label='Cell type'
            size = (0.2*1.25, 0.12)
            margin_top = 0.1
            ylim = (0, 1.01)
        else:
            y_label += ' outputs'
            x_label='Neuron type'
            size = (0.2, 0.12)
            ylim = (0, 1)
            margin_top = 0.05

        ns = list(set(list(nonvariable.keys()) + list(variable.keys())))


        data, l, c = [], [], []
        types = ['sensory', 'modulatory', 'inter', 'motor',]
        if inputs:
            types += ['muscle']
        for typ in types:
            data.append([variable[n]/(nonvariable[n]+variable[n]) for n in ns if ntype(n) == typ])
            c.append(dark_colors[typ])
            l.append(typ.capitalize())

        # for d in data:
        #     print(np.mean(d), len(d))

        stats = ((1,2), (1,3), (4,3), (4, 5), (3,5), (2,5), (5, 1), ) if inputs else ((3, 2), (2, 1), (3, 4), (2, 4), (1, 4))
        self.plt.plot(
            'box_plot', data,
            margin={'left': 0.04, 'right': 0.01, 'top': margin_top, 'bottom': 0.02},
            colors=c, xticklabels=l, y_label=y_label, x_label=x_label, ylim=ylim,
            xtickpad=2, xpad=3,
            show_outliers=False, stats=stats, size=size,
            save=f+'_variable_edge_prevalence_quantified'+('_in' if inputs else ''), darkmedian=True
        )

    def variable_contacts(self, f, adj_pair_classifications, edge_pair_classifications):

        adj_classifications = adj_pair_classifications.copy().rename('classification')
        classifications = edge_pair_classifications.copy().rename('classification')
        
        G = data_manager.get_connections(datasets_with_adj)['count'].copy()
        G = data_manager.to_npair(G)
        G = data_manager.remove_postemb(G)
        G = G[G.sum(axis=1) > 0]
        
        G_adj = data_manager.get_adjacency(datasets_with_adj)
        G_adj = data_manager.to_npair(G_adj)
        G_adj = data_manager.remove_postemb(G_adj)
        G_adj = G_adj[G_adj.sum(axis=1) > 0]

        G_adj = G_adj.join(adj_classifications, how='outer')
        G = G.join(classifications, how='outer')
        
        p_nonvariable_connection_at_nonvariable_contact = []
        p_variable_connection_at_nonvariable_contact = []
        
        for dataset in datasets_with_adj:
            G_dataset = G[G[dataset] > 0]
        
            variable_connections = G_dataset.index[G_dataset['classification'].isin(['noise', 'remainder'])]
            stable_connections = G_dataset.index[G_dataset['classification'].isin(['stable', 'increase', 'decrease'])]
            
            variable_contacts = G_adj.index[G_adj['classification'].isin(['noise', 'remainder'])]
            
            p_nonvariable_connection_at_nonvariable_contact.append(1 - stable_connections.isin(variable_contacts).sum() / stable_connections.size)
            p_variable_connection_at_nonvariable_contact.append(1 - variable_connections.isin(variable_contacts).sum() / variable_connections.size)
            
        xy = (((''), [
            [1, 2],
            [
                np.mean(p_nonvariable_connection_at_nonvariable_contact),
                np.mean(p_variable_connection_at_nonvariable_contact), 
            ]
        ]),)
        colors = ['#666666', '#dddddd']
        self.plt.plot(
            'bar_graph',
            xy,
            colors=colors,
            xticks=[1, 2],
            xticklabels=['Non-variable', 'Variable'],
            x_label='Connection classification',
            y_label='Proportion of connections\nat non-variable physical contacts',
            ylim=(0, 1),
            yticks=(0, 0.5, 1),
            legend=False,
            size=(0.09, 0.135),
            save=f+'_at_variable_contacts',
        )

    def sections_per_synapse(self, f, edge_classifications):

        edge_classifications = edge_classifications.copy()
        synapses = data_manager.get_synapses_one_to_one().copy()
        
        variable_categories = {
            'increase': 'Non-variable',
            'decrease': 'Non-variable',
            'stable': 'Non-variable',
            'noise': 'Variable',
            'remainder': 'Variable'
        }

        df = synapses.merge(
            edge_classifications.map(variable_categories),
            left_on=['pre', 'post'],
            right_index=True
        )

        nonvariable_section_counts = df[(df['edge_classifications'] != 'Variable')]['sections'].value_counts(normalize=True).sort_index()
        variable_section_counts = df[(df['edge_classifications'] == 'Variable')]['sections'].value_counts(normalize=True).sort_index()
        
        xy = [
            (
                str(s.start or s.stop) + ('+' if s.start else ''), 
                [
                    nonvariable_section_counts[s].sum(),
                    variable_section_counts[s].sum(),
                ]
            )
            for s in [slice(5, None), slice(4), slice(3), slice(2), slice(1)]
        ]
        colors = {
            '1': ['#888888', '#ffffff'],
            '2': ['#777777', '#eeeeee'],
            '3': ['#666666', '#dddddd'],
            '4': ['#555555', '#cccccc'],
            '5+': ['#444444', '#bbbbbb'],
            '1_edge': ['#000000', '#000000'],
            '2_edge': ['#000000', '#000000'],
            '3_edge': ['#000000', '#000000'],
            '4_edge': ['#000000', '#000000'],
            '5+_edge': ['#000000', '#000000'],
        }
        self.plt.plot(
            'stacked_bar_graph',
            xy,
            ylim=(0, 1),
            x_label='Connection classification',
            y_label='Proportion of synapses',
            yticks=(0, 0.5, 1),
            colors=colors,
            legend_title='Sections', legendpos='right', legendcol=1, legendreverse=True,
            size=(0.09, 0.135),
            xlabels=['Non-variable', 'Variable'],
            save=f+'_number_of_sections',
        )    
        print('Non-variable one-section synapses:', nonvariable_section_counts[1])
        print('Variable one-section synapses:', variable_section_counts[1])

    def synapse_post_size(self, f, edge_classifications):

        edge_classifications = edge_classifications.copy()
        synapses = data_manager.get_synapses_one_to_one().copy()

        variable_categories = {
            'increase': 'Non-variable',
            'decrease': 'Non-variable',
            'stable': 'Non-variable',
            'noise': 'Variable',
            'remainder': 'Variable'
        }

        df = synapses.merge(
            edge_classifications.map(variable_categories),
            left_on=['pre', 'post'],
            right_index=True
        )

        # Do variable connections take up less of the synapses?
        df_for_plot = df[df['dataset'].isin(datasets_with_adj)]

        data_to_plot = (
            ('Non-variable', df_for_plot.loc[df_for_plot['edge_classifications'] == 'Non-variable', 'post_weight']),
            ('Variable', df_for_plot.loc[df_for_plot['edge_classifications'] == 'Variable', 'post_weight']),
        )
        
        self.plt.plot(
            'kde', data_to_plot, colors=['#000000', '#aaaaaa'],
            size=(0.15, 0.135),
            x_label='Relative postsynaptic\ncontact area at synapse',
            yticks=[], 
            y_label='Frequency', ypad=5,
            linewidth=1.5,
            save=f+'_variable_post_synaptic_size',
            xlim=(0, 1),
            fill=False, clear_x=False,
            margin={'left': 0.05, 'right': 0.01, 'top': 0.05, 'bottom': 0.03},
        )

    def variable_proportion(self, f, edge_pair_classifications, ratio=True):

        G = data_manager.get_connections()['count'].copy()
        G = data_manager.remove_postemb(G)
        G = data_manager.to_npair(G)
        
        cl = edge_pair_classifications.copy()

        irange = range(1, 14)

        types = ['sensory', 'modulatory', 'inter', 'motor']

        x, ys = [], []
        for i in irange:
            x.append(i)
        for t in ['all cells'] + types:

            variable = [G.loc[e][-2:] for e, c in cl.items() if c in ('noise', 'remainder') and (ntype(e[0]) == t or t == 'all cells')]
            nonvariable = [G.loc[e][-2:] for e, c in cl.items() if c in ('stable', 'increase', 'decrease') and (ntype(e[0]) == t or t == 'all cells')]

            y = []
            y2 = []
            for i in irange:

                variable_c = np.sum(np.array(variable) >= i)*0.5
                nonvariable_c = np.sum(np.array(nonvariable) >= i)*0.5
                if ratio:
                    y.append(variable_c/(variable_c+nonvariable_c+0.0))
                else:
                    y.append(nonvariable_c)
                    y2.append(variable_c)

            if ratio:
                ys.append((t.capitalize(), y))
            if not ratio and t == 'all cells':
                ys.append(('Non-variable', y))
                ys.append(('Variable', y2))

        if ratio:
            ylim = (0, 0.7)
            colors = ['k'] + [dark_colors[t] for t in types]
            y_label = 'Proportion of variable\nconnections between cell pairs'
        else:
            ylim = (0, 600)
            colors = ['k', 'gray']
            y_label = 'Connections between cell pairs'

        self.plt.plot(
            'xy_graph', (x, ys), larval_stages=False,
            colors=colors,
            ylim=ylim, xlim=(min(irange), max(irange)), markersize=4,
            xticks=irange, xticklabels=[('' if not i % 2 else i) for i in irange],
            y_label=y_label,
            x_label='Synapse number threshold',
            save=f'{f}_variable_proportion_{int(ratio)}',
            size=0.18, margin={'left': 0.04, 'right': 0.01, 'top': 0.1, 'bottom': 0.05},
        )

    def number_of_postsynaptic_cells(self, f, edge_classifications):

        edge_classifications = edge_classifications.copy()
        synapses = data_manager.get_synapses_one_to_one().copy()

        variable_categories = {
            'increase': 'Non-variable',
            'decrease': 'Non-variable',
            'stable': 'Non-variable',
            'noise': 'Variable',
            'remainder': 'Variable'
        }

        df = synapses.merge(
            edge_classifications.map(variable_categories),
            left_on=['pre', 'post'],
            right_index=True
        )

        # Are variable synapses more likely to have more posts?
        freq_table = pd.crosstab(df['num_post'], df['edge_classifications'])
        freq_table['Variable'] /= freq_table.index  # otherwise dyadic synapses are counted twice, etc
        freq_table['Non-variable'] /= freq_table.index  # otherwise dyadic synapses are counted twice, etc
        freq_table /= freq_table.sum()
        freq_table = freq_table.transpose()
        x = freq_table.columns
        y = tuple(sorted(freq_table.iterrows()))
        colors = {'Non-variable': '#000000', 'Variable': '#aaaaaa'}
        self.plt.plot(
            'xy_graph', (x, y), y_label='Proportion of\nsynapses', x_label='Postsynaptic cells at synapse',
            colors=colors,
            size=0.18,
            linewidth=1.5,
            xticks=range(1, 7), yticks=[0, 0.3, 0.6], xlim=(1, 6.5),
            margin={'left': 0.07, 'right': 0.01, 'top': 0.03, 'bottom': 0.05},
            markersize=5, larval_stages=False, dots=False,
            save=f+'_number_of_postsynaptic_cells',
        )

    def polyadic_partners(self, f, edge_classifications):

        edge_classifications = edge_classifications.copy()
        synapses = data_manager.get_synapses_one_to_one().copy()

        variable_categories = {
            'increase': 'Non-variable',
            'decrease': 'Non-variable',
            'stable': 'Non-variable',
            'noise': 'Variable',
            'remainder': 'Variable'
        }

        df = synapses.merge(
            edge_classifications.map(variable_categories),
            left_on=['pre', 'post'],
            right_index=True
        )

        # Are variable more likely to be part of synapses with stable post?
        df['variable'] = (df['edge_classifications'] == 'Variable').astype(int)
        df_syns = df.groupby(['dataset', 'syn_id'])['variable'].agg(['count', 'sum'])
        df_syns = df_syns.rename({'count': 'post', 'sum': 'variable'}, axis=1)
        df_syns['stable'] = df_syns['post'] - df_syns['variable']
        df_syns = df_syns[df_syns['post'] > 1]
        
        total = df_syns['post'].sum()
        
        stable_total = df_syns['stable'].sum()
        stable_and_stable_neighbor = df_syns[df_syns['stable'] > 1]['stable'].sum()
        variable_total = df_syns['variable'].sum()
        variable_and_stable_neighbor = df_syns[df_syns['stable'] > 0]['variable'].sum()
        
        p_stable = (stable_and_stable_neighbor + variable_and_stable_neighbor) / total
        p_stable_neighbor_given_stable = stable_and_stable_neighbor / stable_total
        p_stable_neighbor_given_variable = variable_and_stable_neighbor / variable_total
        
        xy = (((''), [
            [1, 2],
            [
                p_stable_neighbor_given_stable,
                p_stable_neighbor_given_variable, 
            ]
        ]),)
        colors = ['#666666', '#dddddd']
        self.plt.plot(
            'bar_graph',
            xy,
            colors=colors,
            xticks=[1, 2],
            xticklabels=['Non-variable\nsynapses', 'Variable\nsynapses'],
            y_label='Proportion of synapses with\na non-variable polyadic partner',
            ylim=(0, 1),
            yticks=(0, 0.5, 1),
            legend=False,
            size=(0.09, 0.135),
            save=f+'_polyadic_partners',
        )

    def filling_fraction_per_type(self, f, inputs=False):

        G = data_manager.get_connections()['count'].copy()
        G_adj = data_manager.get_adjacency().copy()

        synapses = defaultdict(float)
        contacts = defaultdict(float)
        for (pre, post), syns in G.iterrows():
            synapses[post if inputs else pre] += 1
        for (pre, post), adj in G_adj.iterrows():
            contacts[pre] += 1

        filling_fractions = defaultdict(list)
        for n in synapses:
            typ = ntype(n)
            if typ == 'other':
                continue
            filling_fractions[typ].append(synapses[n]/contacts[n])

        data, l, c = [], [], []
        types = ['sensory', 'modulatory', 'inter', 'motor']
        if inputs:
            types += ['muscle']
        for typ in types:
            data.append(filling_fractions[typ])
            c.append(dark_colors[typ])
            l.append(typ.capitalize())

        y_label = 'Proportion of physical contacts with synapses'
        x_label = ' cell type'

        if inputs:
            x_label = 'Postsynaptic' + x_label
            size = (0.2*1.25, 0.12)
            ylim = (0, 0.6)
        else:
            x_label = 'Presynaptic' + x_label
            size = (0.2, 0.12)
            ylim = (0, 0.6)

        inputs_or_outputs = 'inputs' if inputs else 'outputs'
        self.plt.plot(
            'box_plot', data, colors=c, darkmedian=True, xticklabels=l, y_label=y_label, ylim=ylim,
            stats=None, size=size, show_outliers=False,
            margin={'left': 0.04, 'right': 0.01, 'top': 0.1, 'bottom': 0.05},
            x_label=x_label, xtickpad=2,
            save=f'{f}_filling_fraction_{inputs_or_outputs}'
        )



    def _edge_types_per_neuron(self, edge_classifications, outputs=True, inputs=False):

        classifications = edge_classifications.copy()
        G = data_manager.get_connections()['count'].copy()
        G = data_manager.remove_postemb(G)

        stable = defaultdict(float)
        variable = defaultdict(float)

        for (pre, post), syns in G.iterrows():

            classification = classifications[(pre, post)]
            if classification in ('stable', 'increase', 'decrease'):
                if inputs:
                    stable[post] += np.count_nonzero(syns[-2:])
                if outputs:
                    stable[pre] += np.count_nonzero(syns[-2:])

            if classification in ('remainder', 'noise'):
                if inputs:
                    variable[post] += np.count_nonzero(syns[-2:])
                if outputs:
                    variable[pre] += np.count_nonzero(syns[-2:])

        for n in stable:
            stable[n] *= 0.5
        for n in variable:
            variable[n] *= 0.5

        return stable, variable

    def variable_correlation(self, f, edge_classifications, versus='synapse_number', box=False):
        
        classifications = edge_classifications.copy()
        G = data_manager.get_connections()['count'].copy()
        G = data_manager.remove_postemb(G)

        y_label = 'Variable connections'

        ylim = (0, 20)
        fname = '_variable_edge_prevalence_vs_' + versus

        if versus == 'synapse_number':
            xlim = (0, 300)
            xticks = range(0, 250, 50)
            x_label = 'Synapses in stable connections'
        if versus == 'stable_growth':
            xlim = (0, 20)
            xticks = range(0, 25, 5)
            x_label = 'Relative number of synapses added to stable connections'

        stable, variable = self._edge_types_per_neuron(classifications)

        stable_synapses_start = defaultdict(float)
        stable_synapses_end = defaultdict(float)
        stable_synapse_increase = {}

        for (pre, post), syns in G.iterrows():
            if classifications[(pre, post)] != 'stable':
                continue
            n = pre
            stable_synapses_start[n] += syns[0]
            stable_synapses_end[n] += np.mean(syns[-2:])

        for n in stable_synapses_start:
            if stable_synapses_start[n] == 0:
                continue
            if versus == 'synapse_number':
                stable_synapse_increase[n] = stable_synapses_end[n]
            if versus == 'stable_growth':
                stable_synapse_increase[n] = stable_synapses_end[n]/stable_synapses_start[n]

        x, y, c, l = [], [], [], []
        types = ['sensory', 'modulatory', 'inter', 'motor',]

        for typ in types:
            x.append([stable_synapse_increase[n] for n in stable_synapse_increase if ntype(n) == typ])
            y.append([variable[n] for n in stable_synapse_increase if ntype(n) == typ])
            c.append(dark_colors[typ])
            l.append(typ.capitalize())

        if box:
            self.plt.plot(
                'box_plot', x, colors=c, yticklabels=l, x_label=x_label, y_label='Presynaptic cell type',
                show_outliers=False, size=(0.2*1.25, 0.07), vert=False, #stats=stats,
                darkmedian=True, xtickpad=2, xlim=xlim,
                xticks=xticks,
                margin={'left': 0.04, 'right': 0.01, 'top': 0.1, 'bottom': 0.05},
                save=f+fname+'_box',
            )
        else:
            self.plt.plot(
                'scatter', [x, y], x_label=x_label, y_label=y_label, colors=c, legend=l,
                size=(0.2*1.25, 0.17), xlim=xlim, ylim=ylim, crop=True, stats='spearman_combined',
                xticks=xticks, yticks=range(0, 25, 5),
                margin={'left': 0.04, 'right': 0.01, 'top': 0.1, 'bottom': 0.05},
                save=f+fname, legendpos='floatright'
            )
