from collections import defaultdict

import numpy as np
import pandas as pd
import scipy as sc
from scipy.stats import linregress

from src.data import data_manager
from src.data.neuron_info import contralateral
from src.data.dataset_info import all_datasets, datasets_with_adj, timepoint

from src.plotting import plotter


class Figure(object):

    def __init__(self, output_path, page_size=7.20472):
        self.plt = plotter.Plotter(output_path=output_path, page_size=page_size)

    def number_of_connections(self, f):

        G = data_manager.get_connections()['count'].copy()
        G = data_manager.remove_postemb(G)

        con_count = G.astype(bool).sum()

        x = [timepoint[d] for d in all_datasets]
        y = [
            ('All cell', con_count)
        ]

        y_label = 'Connections between\ncells existing from birth'
        self.plt.plot(
            'xy_graph', (x, y),
            y_label=y_label, save=f+'_number_of_connections', ylim=(0, 2000),
            linkpoints=False, stats='regression_only', clipon=True,
            size=0.19, x_label='Developmental age',
            margin={'left': 0.045, 'right': 0.015, 'top': 0.01, 'bottom': 0.03},
        )

        connection_at_birth = 0
        connection_after_birth = 0
        for syns in G.values:
            if syns[0] or syns[1] or syns[2] or syns[3]:
                connection_at_birth += np.mean(syns[-2:]) - syns[0]
            else:
                connection_after_birth += np.mean(syns[-2:])

        print('Synapses enforcing existing connections:', connection_at_birth, connection_at_birth/(connection_at_birth+connection_after_birth))
        print('Synapses making new connections:', connection_after_birth, connection_after_birth/(connection_at_birth+connection_after_birth))

        print('{:.2f} -> {:.2f}'.format(con_count[0], con_count[-1]))
        print('Factor: {:.2f}'.format(con_count[-1]/con_count[0]))

    def edge_changes(self, f):

        G = data_manager.get_connections()['count'].copy()
        G = data_manager.remove_postemb(G)
        G_from_birth =  G[G['Dataset1'] > 0]


        syn_count = G_from_birth.apply(lambda s: s[s > 0].mean())
        syn_err = G_from_birth.apply(lambda s: sc.stats.sem(s[s > 0]))

        x = [timepoint[d] for d in all_datasets]
        y = [
            ('All cell', syn_count)
        ]

        y_label = 'Synapses per connection\nexisting from birth'
        self.plt.plot('xy_graph', (x, y),# err=[syn_err],
            y_label=y_label, save=f+'_synapses_per_connection', ylim=(0, 8),
            linkpoints=False, stats='regression_only', clipon=True,
            size=0.19, x_label='Developmental age',
            margin={'left': 0.045, 'right': 0.015, 'top': 0.01, 'bottom': 0.03},
        )



    def new_edge_contacts(self, f):
        
        G = data_manager.get_connections()['count'].copy()
        G = data_manager.to_npair(G)
        G = data_manager.remove_postemb(G)
        
        G_adj = data_manager.get_adjacency().copy()
        G_adj = data_manager.to_npair(G_adj)
        G_adj = data_manager.remove_postemb(G_adj)
        
        G_merge = G.merge(G_adj, 'outer', left_index=True, right_index=True, suffixes=('_syn', '_adj'))
        
        # Filter to available contacts at birth and in adults.
        G_merge = G_merge[(G_merge['Dataset1_adj'] > 0) & (G_merge['Dataset8_adj'] > 0)]
        
        # .. without a connection already at birth.
        G_merge = G_merge.fillna(0)
        G_merge = G_merge[(G_merge['Dataset1_syn'] == 0) & (G_merge['Dataset2_syn'] == 0)]
        
        # Determine if a new connection is built.
        G_merge['new_connection'] = (G_merge['Dataset7_syn'] > 0) & (G_merge['Dataset8_syn'] > 0)
        
        # Bin to 0.1 um^3
        bins = list(np.arange(0, 1.05, 0.1)) + [np.inf]
        G_merge['Dataset1_adj_binned'] = pd.cut(
            G_merge['Dataset1_adj']/10**6, 
            bins
        )
        
        # Calculate probability of new connection.
        bin_midpoints = [b+0.05 for b in bins[:-1]]
        prob_of_new_con = G_merge.groupby('Dataset1_adj_binned')['new_connection'].mean()
        
        

        xticks = np.arange(0, 1.11, 0.1)
        xticklabels = ['' if i%2 else f'{x:.1f}' for i, x in enumerate(xticks)]

        self.plt.plot(
            'scatter', ([bin_midpoints], [prob_of_new_con*100]), colors=['k'],
            x_label=u'Contact area at birth (μm$^2$)',
            y_label='Probability of new connection (%)',
            size=0.19,
            ylim=(0, 30), xlim=(0, 1.05),
            xgrid=True, stats='spearman_combined',
            xticks=xticks, xticklabels=xticklabels,
            yticks=range(0, 31, 10),
            marker='.', markersize=12, markeredgewidth=1, alpha=1,
            margin={'left': 0.045, 'right': 0.015, 'top': 0.01, 'bottom': 0.03},
            save=f+'_connection_probability'
        )



    def correlate_centrality_with_growth(self, f, x_axis='degree', increase=None, input_output=None):

        assert x_axis in ('degree', 'degree_in', 'degree_out')
        assert increase in ('degree', 'synapse', 'synapse_size')
        assert input_output in ('output', 'input')

        G = data_manager.get_connections()['size' if increase == 'synapse_size' else 'count'].copy()

        degree_out = G.astype(bool).groupby('pre').sum()
        degree_in = G.astype(bool).groupby('post').sum()
        degree = degree_out.add(degree_in, fill_value=0)
        
        synapse_out = G[G['Dataset1'] > 0].groupby('pre').sum()
        synapse_in = G[G['Dataset1'] > 0].groupby('post').sum()

        degree_out = degree_out[degree_out['Dataset1'] > 0]
        degree_in = degree_in[degree_in['Dataset1'] > 0]
        synapse_out = synapse_out[synapse_out['Dataset1'] > 0]
        synapse_in = synapse_in[synapse_in['Dataset1'] > 0]
        
        if increase == 'degree':
            if input_output == 'input':
                df = degree_in
            else:
                df = degree_out
            degree_increase = df[['Dataset7', 'Dataset8']].mean(axis=1) - df['Dataset1']
            y = degree_increase
        
        if increase in ('synapse', 'synapse_size'):
            if input_output == 'input':
                df = synapse_in
            else:
                df = synapse_out
            synapse_increase = df[['Dataset7', 'Dataset8']].mean(axis=1) - df['Dataset1']
            y = synapse_increase / degree['Dataset1'].loc[synapse_increase.index]
            
        
        

        x = degree['Dataset1']
        if x_axis == 'degree_in':
            x = degree_in['Dataset1']
        if x_axis == 'degree_out':
            x = degree_out['Dataset1']

        neurons = list(set(x.index).intersection(y.index))
        x = x[neurons]
        y = y[neurons]
        
        x_label = 'Connections at birth (#)'
        if x_axis == 'degree_in':
            x_label = 'In-connections at birth (#)'
        if x_axis == 'degree_out':
            x_label = 'Out-connections at birth (#)'
        y_label = ''
        xlim = (0, 35)
        ylim = None
        xticks = (0, 15, 30)
        yticks = None
        colors = [(1, 0.5, 0.5) if input_output == 'input' else '#86A2A4']

        if increase == 'degree':
            y_label = 'New connections (#)'
            if input_output == 'input':
                ylim = (-2, 20)
                yticks = (0, 10, 20)
            else:
                ylim = (-3, 30)
                yticks = (0, 15, 30)

        if increase == 'synapse':
            y_label = 'Synapses added per existing connection (#)'
            if input_output == 'input':
                ylim = (-2, 12)
                yticks = (0, 6, 12)
            else:
                ylim = (-1, 12)
                yticks = (0, 6, 12)
                
        if increase == 'synapse_size':
            y /= 1e+6
            y_label = 'Synapse volume added per existing connection (10$^6$ nm$^3$)'
            if input_output == 'input':
                ylim = (-0.2, 1.4)
                yticks = (0, 0.7, 1.4)
            else:
                ylim = (-0.2, 1.4)
                yticks = (0, 0.7, 1.4)
        
        size = (0.09, 0.09)
        if x_axis in ('degree_in', 'degree_out'):
            xticks = (0, 10, 20)
            xlim = (0, 20)
            size = (0.07, 0.09)
        


        self.plt.plot(
            'scatter', (x, y), stats='spearman_combined',
            size=size, alpha=0.7,
            colors=colors,
            margin={'left': 0.05, 'right': 0.01, 'top': 0.01, 'bottom': 0.03},
            legend=None,
            marker='.', markersize=3, markeredgewidth=0.7,

            x_label=x_label, y_label=y_label, xlim=xlim, ylim=ylim, crop=True, xticks=xticks, yticks=yticks,
            save=f+'_correlate_centrality_' + input_output + '_' + increase,
        )





    def _stable_edge_pairs(self):

        G = data_manager.get_connections()['count'].copy()
        
        # Restrict to edges present across development.
        G = G[G.astype(bool).sum(axis=1) == G.shape[1]]
        
        edge_pairs = {}
        for edge1 in G.index:
            for edge2 in G.index:
                if edge1 == edge2 or edge1[0] == edge1[1] or edge2[0] == edge2[1]:
                    continue
                if (edge2, edge1) in edge_pairs:
                    continue
                if edge1[0] == edge2[0]:
                    typ = 'shared_pre'
                elif edge1[1] == edge2[1]:
                    typ = 'shared_post'
                else:
                    typ = 'other'

                edge_pairs[(edge1, edge2)] = typ

        return edge_pairs

    def hebbian_plasticity(self, f, quantify_slope=False, relative=False, use_size=False):
        
        datasets = list(datasets_with_adj if use_size else all_datasets)
        
        G = data_manager.get_connections()['size' if use_size else 'count'][datasets].copy()

        edge_pairs = self._stable_edge_pairs()

        times = [timepoint[d] for d in datasets]

        result = defaultdict(list)
        for i, d in enumerate(datasets + ['slope']):
            values = {}
            for edge in set([e for ep in edge_pairs for e in ep]):
                if d == 'slope':
                    if G.loc[edge][0]:
                        values[edge] = np.mean(G.loc[edge][-2:])/G.loc[edge][0]
                else:
                    values[edge] = float(G.loc[edge][i])


            cvs = defaultdict(list)
            for (edge1, edge2), typ in edge_pairs.items():
                s1 = values[edge1]
                s2 = values[edge2]
                if np.isnan(s1) or np.isinf(s1):
                    continue
                if np.isnan(s2) or np.isinf(s2):
                    continue
                if s1 + s2 == 0:
                    continue
                cv = abs(s1 - s2) / (s1 + s2)
                cvs[typ].append(cv)

            if d == 'slope':
                result['slope'].append(cvs['shared_pre'])
                result['slope'].append(cvs['shared_post'])
                result['slope'].append(cvs['other'])
            else:
                result['shared_pre'].append(cvs['shared_pre'])
                result['shared_post'].append(cvs['shared_post'])
                result['other'].append(cvs['other'])


#        print('Number of edges:')
#        print(len(result['shared_pre'][0]), len(result['shared_post'][0]), len(result['other'][0]))

        a1 = 'Outputs from\nthe same cell'
        a2 = 'Inputs to\nthe same cell'
        a3 = 'Connections to and\nfrom different cells'

        color = {a1: (73/255.0, 187/255.0, 168/255.0),
                 a2: (0, 0, 0),
                 a3: (0.5, 0.5, 0.5)}

        for c in list(color.keys()):
            color[c+'_edge'] = color[c]

        if not quantify_slope:

            means = [1]*len(datasets)
            if use_size:
                y_label = 'Coefficient of variation\n(CV) in synapse volume'
                ylim = (0.3, 0.5)
                yticks = (0.3, 0.4, 0.5)
            elif relative:
                means = np.array([np.mean(v) for v in result['other']])
                y_label = 'Relative variation of\nsynapse number'
                ylim = (0.7, 1.1)
                yticks = (0.7, 0.8, 0.9, 1.0, 1.1)
            else:
                y_label = 'Coefficient of variation\n(CV) in synapse number'
                ylim = (0.2, 0.4)
                yticks = (0.2, 0.3, 0.4)

            mean = ((a1, [(np.mean(vs)/means[i]) for i, vs in enumerate(result['shared_pre'])]),
                    (a2, [(np.mean(vs)/means[i]) for i, vs in enumerate(result['shared_post'])]),
                    (a3, [(np.mean(vs)/means[i]) for i, vs in enumerate(result['other'])]))
            err = ([(sc.stats.sem(vs)/means[i]) for i, vs in enumerate(result['shared_pre'])],
                   [(sc.stats.sem(vs)/means[i]) for i, vs in enumerate(result['shared_post'])],
                   [(sc.stats.sem(vs)/means[i]) for i, vs in enumerate(result['other'])])
    
            self.plt.plot(
                'xy_graph', (times, mean), colors=color, show=True, err=err,
                stats='regression_only', linkpoints=False, ylim=ylim, yticks=yticks,
                y_label=y_label,
                size=(0.19, 0.11), clipon=True,
                margin={'left': 0.045, 'right': 0.015, 'top': 0.01, 'bottom': 0.03},
                no_legend=True, #xticks='remove',
                x_label='Developmental age',
                save=f+'_hebbian_plasticity',
            )


            bars = [
                ('', [
                    times,
                    [(mean[1][1][i]-mean[0][1][i])/mean[2][1][i] for i in range(len(datasets))]
                ]),
            ]

            reg = linregress(*bars[0][1])
            a = reg.slope
            b = reg.intercept

            cor = sc.stats.spearmanr(*bars[0][1])
            print('spearmanr', cor.pvalue, cor.correlation, len(bars[0][1][0]))

            x_fit = np.arange(0.01, 55, 0.01)
            y_fit = [b + a * x for x in x_fit]
            customstats = (x_fit, y_fit, 0, 0, '')

            for series in bars:
                series[1][0] = series[1][0][:7]
                series[1][1] = series[1][1][:6] + [np.mean(series[1][1][6:])]

            self.plt.plot(
                'bar_graph', bars,
                size=(0.19, 0.08), width=1.5,
                margin={'left': 0.045, 'right': 0.015, 'top': 0.01, 'bottom': 0.03},
                ylim=(-0.05, 0.2),
                yticks=(0, 0.1, 0.2),
                y_label='(post-pre)/not',
                xintersect=0,
                colors=[[np.mean(cs) for cs in zip(color[a1], color[a2])]],
                xticks='remove',
                customstats=customstats,
                legend=False, larval_stages=True,
                save=f+'_hebbian_plasticity_inset',
            )

        else:

            y_label = 'Coefficent of variation (CV) in\nstrengthening of connections'
            x_label = 'Connection pairs'
            types_pretty = (a1, a2, a3)

            boxes = result['slope']

            c = [color[t] for t in types_pretty]

            self.plt.plot(
                'box_plot', boxes, colors=c, xticklabels=[t for t in types_pretty], ylim=(-0.05, 1),
                show_outliers=False, stats=((2, 3), (1, 2), (1, 3)), y_label=y_label, x_label=x_label,
                xtickpad=2, xpad=3,
                size=(0.25, 0.16),
                margin={'left': 0.05, 'right': 0.01, 'top': 0.05, 'bottom': 0.03},
                save=f+'_hebbian_plasticity_slope', darkmedian=False,
            )


    def connectivity_gaps(self, f):

        G = data_manager.get_connections()['count'].copy()

        y_label = 'Missing connections'

        xs = [timepoint[d] for d in all_datasets]
        ys = []
        for d in all_datasets:

            gaps = 0
            for edge in G.index:

                edge_c = (contralateral(edge[0]), contralateral(edge[1]))

                s = G.loc[edge]
                if edge_c in G.index:
                    s_c = G.loc[edge_c]
                else:
                    s_c = pd.Series([0]*len(all_datasets), index=all_datasets)

                if s[d] == 0 and s_c[d] > 0 and np.count_nonzero(s) == 7:
                    gaps += 1

            ys.append(gaps)

        self.plt.plot(
            'xy_graph', (xs, (('', ys),)),
            y_label=y_label, ylim=(0, 50), yticks=(0, 25, 50),
            smooth=True, save=f+'_connectivity_gaps',
            size=0.18, x_label='Developmental age',
            margin={'left': 0.05, 'right': 0.01, 'top': 0.03, 'bottom': 0.03},
        )



    def asymmetry(self, f):

        G = data_manager.get_connections()['count'].copy()
        G = data_manager.remove_noise(G)

        xs = [timepoint[d] for d in all_datasets]
        ys = []
        err = []
        for d in all_datasets:

            done = set()
            cvs = []
            for edge in G.index:

                if edge in done:
                    continue

                edge_c = (contralateral(edge[0]), contralateral(edge[1]))
                if edge == edge_c:
                    continue
                s1 = G.loc[edge, d]
                s2 = 0
                if edge_c in G.index:
                    s2 = G.loc[edge_c, d]

                if s1 == 0 and s2 == 0:
                    continue

                if G.loc[edge].astype(bool).sum() < 7:
                    continue

                cvs.append(abs(s1-s2)/(s1+s2))

                done.add(edge)
                done.add(edge_c)

            ys.append(np.mean(cvs))
            err.append(sc.stats.sem(cvs))

        self.plt.plot(
            'xy_graph', (xs, (('', ys),)),
            ylim=(0, 0.5), err=[err], colors=['k'],
            smooth=True,y_label='Left-right asymmetry (CV)',
            save=f+'_asymmetry',
            size=0.18, x_label='Developmental age',
            margin={'left': 0.05, 'right': 0.01, 'top': 0.03, 'bottom': 0.03},
        )

    def degree_distribution(self, f):

        G = data_manager.get_connections()['count'].copy()

        degree_l1 = defaultdict(int)
        degree_adult = defaultdict(int)

        for (pre, post), syn in G.iterrows():

            if syn[0]:
                degree_l1[post] += 1
                degree_l1[pre] += 1

            if syn[-1]:
                degree_adult[post] += np.count_nonzero(syn[-2:])
                degree_adult[pre] += np.count_nonzero(syn[-2:])

        data = (
            ('', list(degree_l1.values())),
#            ('ad', np.array(degree_adult.values())/2.0),
        )

        self.plt.plot(
            'hist', data, colors=[(0.3, 0.3, 0.3, 0.3)],
            size=0.18, xlim=(0, 35), x_label='Connections per neuron at birth',
            y_label='Neurons',
            fitbins=True,
            save=f+'_degree_distribution',
            margin={'left': 0.05, 'right': 0.01, 'top': 0.05, 'bottom': 0.03},
        )
