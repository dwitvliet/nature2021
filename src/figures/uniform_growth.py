from collections import defaultdict

import numpy as np

from src.data import data_manager
from src.data.dataset_info import all_datasets, datasets_with_adj, timepoint
from src.data.neuron_info import neuron_list, class_members, in_brain, ntype, is_postemb

from src.plotting import plotter


class Figure(object):
    
    def __init__(self, output_path, page_size=7.20472):
        self.plt = plotter.Plotter(output_path=output_path, page_size=page_size)

    def cable_length(self, f):
        
        skeletons = data_manager.get_skeletons()
        
        length = [sum([s['length']/1000000.0 for s in skeletons[d].values()]) for d in all_datasets]

        x = [timepoint[d] for d in all_datasets]
        y = [
            ('All cells', length),
        ]

        y_label = 'Neurite length (mm)'
        self.plt.plot('xy_graph', 
            (x, y), y_label=y_label, ylim=(0, 12), yticks=(0, 4, 8, 12), 
            save=f+'_cable_length', 
            size=0.15, x_label='Developmental age',
            margin={'left': 0.045, 'right': 0.005, 'top': 0.02, 'bottom': 0.03}
        )
        
        print('{:.2f} -> {:.2f}'.format(length[0], length[-1]))
        increase = [skeletons['Dataset8'][n]['length']/skeletons['Dataset1'][n]['length'] for n in skeletons['Dataset8'] if n in skeletons['Dataset1'] and skeletons['Dataset1'][n]['length'] > 0]
        print('Mean increase: {:.2f} +- {:.2f}'.format(np.mean(increase), np.std(increase)))

    def proportion_of_stable_contacts(self, f, uniquecontact=False):
        
        G_adj = data_manager.get_adjacency().copy()
        G_adj = data_manager.to_npair(G_adj)
        G_adj = data_manager.remove_postemb(G_adj)
        
        adj_stable = np.zeros(7)
        adj_nonstable = np.zeros(7)
        
        for adj in G_adj.values:
            a = adj[[0, 1, 2, 3, 4, 5, 7]]
            if uniquecontact:
                a = a.astype(bool)
            if np.count_nonzero(a) == 7:
                adj_stable += a
            else:
                adj_nonstable += a
                
        # print(adj_stable)
        # print(adj_nonstable)
        # print(adj_stable/(adj_stable+adj_nonstable))

        x = [timepoint[d] for d in all_datasets][:7]
        y = [
            ('Present throughout life', adj_stable/(adj_stable+adj_nonstable)*100),
        ]
        
        self.plt.plot('xy_graph', 
            (x, y),
            size=0.15, x_label='Developmental age',
            ylim=(0, 100),
            y_label='Persistent physical contact (%)',
            save=f+'_stable_contacts', legendcol=1, rev_legend=True,
            margin={'left': 0.045, 'right': 0.005, 'top': 0.02, 'bottom': 0.03}
            
        )

    def number_of_synapses(self, f):
        
        G = data_manager.get_connections()['count'].copy()
        syn_count = G.sum()
        
        x = [timepoint[d] for d in all_datasets]
        y = [
            ('All cell', syn_count.to_numpy()),
        ]
        
        y_label = 'Synapses (#)'
        self.plt.plot('xy_graph', 
            (x, y), y_label=y_label, ylim=(0, 9000), yticks=(0, 3000, 6000, 9000), 
            save=f+'_number_of_synapses', 
            size=0.15, x_label='Developmental age',
            margin={'left': 0.045, 'right': 0.005, 'top': 0.02, 'bottom': 0.03}
        )
        
        print('{:.2f} -> {:.2f}'.format(syn_count[0], syn_count[-1]))
        G_nonzero_at_birth = G[G['Dataset1'] > 0]
        increase = G_nonzero_at_birth['Dataset8'] / G_nonzero_at_birth['Dataset1']
        print('Mean factor: {:.2f} +- {:.2f}'.format(np.mean(increase), np.std(increase)))

    def synapse_density(self, f, synapse_type='count'):
        
        G = data_manager.get_connections()[synapse_type].copy()
        G = data_manager.remove_postemb(G)
        
        skeletons = data_manager.get_skeletons()
        
        syn_count = G.sum()
        length = np.array([sum([s['length']/1000 for s in skeletons[d].values() if not is_postemb(s['name'])]) for d in all_datasets])
        
        xs = [timepoint[d] for d in all_datasets]
        ys = [
            ('No post-embryonic', (syn_count/length).to_numpy())
        ]
        y_label = u'Synapse density (#/μm)'
        self.plt.plot('xy_graph', 
            (xs, ys), y_label=y_label, ylim=(0, 1.2), yticks=(0, 0.4, 0.8, 1.2), 
            save=f+'_synapse_density', 
            size=0.15, x_label='Developmental age',
            margin={'left': 0.045, 'right': 0.005, 'top': 0.02, 'bottom': 0.03}
        )

    def branches_correlation(self, f):
        skeletons = data_manager.get_skeletons()
        cells_with_skeletons = [n for c in neuron_list for n in class_members(c) if (in_brain(c) and ntype(c) not in ('other', 'muscle') and not is_postemb(c))]
        branches_l1 = defaultdict(list)
        branches_adult = defaultdict(list)
        for d in all_datasets:
            for n in cells_with_skeletons:
                if ntype(n) in ('other', 'muscle'):
                    continue

                for branch in skeletons[d][n]['branches']:
                    distances = [skeletons[d][n]['dist_to_root'][t] for t in branch]
                    l = max(distances) - min(distances)
                    if d == 'Dataset1':
                        branches_l1[n].append(l)
                    if d == 'Dataset8':
                        branches_adult[n].append(l)

        x, y = [], []
        for n in cells_with_skeletons:
            if ntype(n) in ('other', 'muscle'):
                continue

            length_l1 = sum(branches_l1[n])
            length_adult = sum(branches_adult[n])
            if length_l1 == 0 or length_adult == 0:
                continue
            for i in range(max(len(branches_l1[n]), len(branches_adult[n]))):
                if i >= len(branches_l1[n]):
                    length_rel_l1 = 0
                else:
                    length_rel_l1 = branches_l1[n][i]/length_l1
                if i >= len(branches_adult[n]):
                    length_rel_adult = 0
                else:
                    length_rel_adult = branches_adult[n][i]/length_adult
                x.append(length_rel_l1)
                y.append(length_rel_adult)
                
                if length_rel_l1 > 0.7 and length_rel_adult < 0.1:
                    print(n)

        c = ['#777777']

        self.plt.plot(
            'scatter', [[x], [y]], colors=c, stats='spearman_combined', crop=True,
            x_label='Relative neurite length at birth',
            y_label='Relative neurite length in adult',
            xticks=[0, 0.5, 1], yticks=[0, 0.5, 1],
            alpha=0.3, markersize=3, marker='.',
            xlim=(0, 1), ylim=(0, 1),
            save=f+'_cable_correlation',
            size=(0.15, 0.15),
            margin={'left': 0.04, 'right': 0.01, 'top': 0.01, 'bottom': 0.04},

        )

    def proportion_of_contacts_with_synapses(self, f):
        G = data_manager.get_connections(datasets_with_adj)['count'].copy()
        G_adj = data_manager.get_adjacency(datasets_with_adj).copy()
        
        G = data_manager.remove_postemb(G)
        G_adj = data_manager.remove_postemb(G_adj)
        
        filling_fraction = G.astype(bool)[G_adj.astype(bool)].sum()/G_adj.astype(bool).sum()
    
        print('Filling fractions:', ', '.join([f'{ff:.2f}' for ff in filling_fraction.to_numpy()]))

        x = [timepoint[d] for d in datasets_with_adj]
        y = [
            ('', filling_fraction.to_numpy()),
        ]
        
        self.plt.plot('xy_graph', 
            (x, y),
            size=(0.15, 0.15), x_label='Developmental age',
            ylim=(0, 0.2), 
            y_label='Proportion of contacts\nwith one or more synapses',
            save=f+'_proportion_of_contacts_with_synapses',
            margin={'left': 0.045, 'right': 0.005, 'top': 0.02, 'bottom': 0.03}
            
        )
