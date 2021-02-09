import os

import numpy as np
import scipy as sc
import networkx as nx

from src.data import data_manager
from src.data.neuron_info import ntype, is_neuron, is_postemb
from src.data.dataset_info import all_datasets


save_to = 'figures/exports'


def export_graphs_for_cytoscape(edge_classifications):

    categories = {
        'increase': 'Developmental change',
        'decrease': 'Developmental change',
        'stable': 'Stable',
        'noise': 'Variable',
        'remainder': 'Variable'
    }

    classifications = edge_classifications

    G = data_manager.get_connections()['count'].copy()
    G_normalized = G / G.sum() * G.sum().mean()

    graph = nx.DiGraph()

    for (pre, post), syns in G.iterrows():
        if pre == post:
            continue
        if is_postemb(pre) or is_postemb(post):
            classification = 'Postemb'
        else:
            classification = categories[classifications[(pre, post)]]

        edge_properties = {
            'weight': max(syns), 'classification': classification,
            'weight_normalized': G_normalized.loc[(pre, post)].max(),
            'transparency_stable': max(syns) if classification == 'Stable' else 0,
            'transparency_variable': max(syns) if classification == 'Variable' else 0,
            'transparency_changing': max(syns) if classification == 'Developmental change' else 0,
            'nomodule': int(is_postemb(pre) or is_postemb(post) or ntype(post) == 'other')
        }

        for i, d in enumerate(all_datasets):
            edge_properties[d] = syns[i]

        graph.add_edge(pre, post, **edge_properties)

    for n in ('CANL', 'CANR', 'excgl'):
        if graph.has_node(n):
            graph.remove_node(n)

    modules, module_names = data_manager.get_communities()

    for n in graph.nodes():
        typ = ntype(n)
        graph.nodes[n]['type'] = ntype(n)
        graph.nodes[n]['celltype'] = 'neuron' if is_neuron(n) else typ
        graph.nodes[n]['is_postemb'] = int(is_postemb(n))
        graph.nodes[n]['is_neuron'] = int(is_neuron(n))
        graph.nodes[n]['no_module'] = int(is_postemb(n) or typ == 'other')
        graph.nodes[n]['hidden_l2-l3'] = int(n in ('HSNL', 'HSNR', 'PVNL', 'PVNR', 'PLNL', 'PLNR'))
        graph.nodes[n]['hidden_latel1'] = int(n in ('HSNL', 'HSNR', 'PVNL', 'PVNR', 'PLNL', 'PLNR', 'AVFL', 'AVFR', 'AVM', 'RMFL' 'RMFR'))
        graph.nodes[n]['hidden_l1'] = int(is_postemb(n))
        for module, mname in zip(modules['Dataset8'], module_names):
            if n in module:
                graph.nodes[n]['module'] = mname

    graphs = {'combined': graph}
    for d in all_datasets:
        graphs[d] = graph.copy()
        for (pre, post) in list(graphs[d].edges()):
            weight = graphs[d][pre][post][d]
            if weight == 0:
                graphs[d].remove_edge(pre, post)
            else:
                graphs[d][pre][post]['weight'] = weight
        graphs[d].remove_nodes_from(list(nx.isolates(graphs[d])))

    for dataset, G in graphs.items():
        #number of nodes
        nodelist = list(G.nodes())
        A = np.array(nx.adjacency_matrix(G, nodelist=nodelist, weight='weight_normalized').todense()).astype(float)

        #symmetrize the adjacency matrix
        c = (A + np.transpose(A))/2.0

        #degree matrix
        d = np.diag(np.sum(c, axis=0))
        df = sc.linalg.fractional_matrix_power(d, -0.5)

        #Laplacian matrix
        l = d - c

        #compute the vertical coordinates
        b = np.sum(c * np.sign(A - np.transpose(A)), 1)
        z = np.matmul(np.linalg.pinv(l), b)

        #degree-normalized graph Laplacian
        q = np.matmul(np.matmul(df, l), df)

        #coordinates in plane are eigenvectors of degree-normalized graph Laplacian
        _, vx = np.linalg.eig(q)
        x = np.matmul(df, vx[:,1])
        y = np.matmul(df, vx[:,2])

        for n in graph.nodes():
            key = '' if dataset == 'combined' else dataset + '_'
            if n not in nodelist:
                graph.nodes[n][key+'x'] = 0.0
                graph.nodes[n][key+'y'] = 0.0
                graph.nodes[n][key+'z'] = 0.0
                graph.nodes[n][key+'zorder'] = 0.0
                graph.nodes[n][key+'visibility'] = 0
                continue

            i = nodelist.index(n)
            graph.nodes[n][key+'x'] = x[i] * 11500.0
            graph.nodes[n][key+'y'] = y[i] * 11500.0
            graph.nodes[n][key+'z'] = z[i] * -150.0
            graph.nodes[n][key+'zorder'] = -z[i] if ntype(n) != 'other' else -z[i]-100000  # put glia in back
            graph.nodes[n][key+'visibility'] = 1

    print('x range:', min(x), max(x))
    print('y range:', min(z), max(z))

    for pre, post in graph.edges():
        graph[pre][post]['bend'] = int(graph.nodes[pre]['z'] < graph.nodes[post]['z'])

    fpath = os.path.join(save_to, 'graphs_for_cytoscape.graphml')
    nx.write_graphml(graph, fpath)
    print(f'Saved to `{fpath}`')
