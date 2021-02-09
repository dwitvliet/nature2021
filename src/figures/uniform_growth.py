import numpy as np

from src.data import data_manager
from src.data.dataset_info import all_datasets, timepoint
from src.data.neuron_info import is_postemb

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
