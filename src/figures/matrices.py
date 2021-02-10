import numpy as np

from src.data import data_manager
from src.data.neuron_info import neuron_list, class_members, in_brain, is_postemb, ntype, npair
from src.data.dataset_info import datasets_with_adj

from src.plotting import plotter


class Figure(object):

    def __init__(self, output_path, page_size=7.20472):
        self.plt = plotter.Plotter(output_path=output_path, page_size=page_size)
    
    def adjacency_matrix(self, f, edge_pair_classifications, contact=False):

        if contact:
            G = data_manager.get_adjacency(datasets_with_adj).copy()
        else:
            G = data_manager.get_connections()['count'].copy()
        
        G_pair = data_manager.to_npair(G)
        
        classification = edge_pair_classifications

        neurons = set([npair(n) for c in neuron_list for n in class_members(c) if (in_brain(c) and not is_postemb(c) and c not in ('excgl', 'CAN'))])
        neuron_order = {'sensory': 0, 'modulatory': 1, 'inter': 2, 'motor': 3, 'muscle': 4, 'other': 5}
        
        if contact:
            post_neurons = sorted(neurons, reverse=True, key=lambda n: (neuron_order[ntype(n)], n))
            post_neurons = [n for n in post_neurons if ntype(n) not in ('other',)]
            pre_neurons = post_neurons[::-1]

        else:
            post_neurons = sorted(neurons, reverse=True, key=lambda n: (neuron_order[ntype(n)], n))
            pre_neurons = [n for n in post_neurons if ntype(n) not in ('other', 'muscle')][::-1]

        data = np.zeros((len(pre_neurons), len(post_neurons)), dtype=float)
        colors = np.ones((len(pre_neurons), len(post_neurons), 3), dtype=float)
        borders = np.ones((len(pre_neurons), len(post_neurons), 3), dtype=float)
        text = np.empty((len(pre_neurons), len(post_neurons)), dtype='str')

        # Set colors by classifications.
        for i, pre in enumerate(pre_neurons):
            for j, post in enumerate(post_neurons):
                if (pre, post) not in G_pair.index:
                    continue
                
                edge = G_pair.loc[(pre, post)]
                if edge.sum() == 0:
                    continue

                if contact:
                    data[i][j] = min(edge.max(), 1e7)
                    colors[i][j] = [1 - (edge > 0).sum() / edge.count()] * 3
                    
                else:
                    data[i][j] = min(edge.max(), 30)
 
                    colors[i][j] = (1, 1, 1)
                    colors[i][j] = {
                        'noise': (0.8, 0.8, 0.8),
                        'remainder': (0.8, 0.8, 0.8),
                        'stable': (0, 0, 0),
                        'increase': (0, 0.5, 1),
                        'decrease': (0, 0.5, 1),
                    }[classification[(pre, post)]]
                    
                    if classification[(pre, post)] == 'decrease':
                        text[i][j] = '.'

        # Insert gaps between neuron types.
        types = [ntype(n) for n in post_neurons]
        for t in set(types):
            i = types.index(t)
            if i == 0:
                continue
            types.insert(i, '')
            post_neurons.insert(i, '')
            data = np.insert(data, i, 0, axis=1)
            colors = np.insert(colors, i, 1, axis=1)
            text = np.insert(text, i, '', axis=1)

        post_spacers = sorted([i for i, t in enumerate(types) if not t])

        types = [ntype(n) for n in pre_neurons]
        for t in set(types):
            i = types.index(t)
            if i == 0:
                continue
            types.insert(i, '')
            pre_neurons.insert(i, '')
            data = np.insert(data, i, 0, axis=0)
            colors = np.insert(colors, i, 1, axis=0)
            text = np.insert(text, i, '', axis=0)

        pre_spacers = sorted([i for i, t in enumerate(types) if not t])

        if contact:
            data = data/np.max(data)*1.2
            x_label = ''
            y_label = ''
        else:
            data = data/np.max(data)*1.2
            x_label = 'Presynaptic neuron pair'
            y_label = 'Postsynaptic cell pair'


        self.plt.plot(
            'adjacency_matrix', data, colors=colors, borders=None, no_x=True, no_y=True,
            yticklabels=post_neurons, xticklabels=pre_neurons, text=text,
            yticks=range(len(post_neurons)), xticks=range(len(pre_neurons)),
            post_spacers=post_spacers, pre_spacers=pre_spacers,
            xtickpad=2, ytickpad=2, xpad=7, ypad=7,
            x_label=x_label, y_label=y_label,
            size=(0.92, 0.92/len(pre_neurons)*len(post_neurons)),
            margin={'left': 0.075, 'right': 0.005, 'top': 0.075, 'bottom': 0.005},
            save=f+'_adjacency_matrix'
        )