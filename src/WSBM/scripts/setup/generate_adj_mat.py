#!/usr/bin/env python
import os
import sys

import scipy.io
sys.path.append('/home/witvliet/worm-connectomics/')
#sys.path.append('/home/witvliet/worm-connectomics/')

from _data.dataset_info import all_datasets
from _data import data_manager
from _data.neuron_info import is_postemb, ntype, contralateral


DATA_DIR = '../../data/'





def generate_pair_info_file(node_list, filename):
    pairs = {}
    node_to_int = {n: i + 1 for i, n in enumerate(node_list)}
    for node in node_to_int:
        if node not in pairs and contralateral(node) not in pairs:
            if node != contralateral(node):
                pairs[node] = contralateral(node)
    with open(filename, "w") as f:
        for k, v in pairs.items():
            f.write("{},{}\n".format(node_to_int[k], node_to_int[v]))





def save_to_mat(arr, name):
    scipy.io.savemat(os.path.join(DATA_DIR, f'{name}.mat'), {'arr': arr})


if __name__ == "__main__":
    
    # Load data.
    connectivity = data_manager.get_connections()['count']
    
    for dataset_i, dataset in enumerate(connectivity.columns):
        
        matrix = data_manager.series_to_dense_matrix(connectivity[dataset]).sort_index().sort_index(axis=1)
        
        # Drop glia
        glia = [n for n in matrix.columns if ntype(n) == 'other']
        matrix = matrix.drop(glia, axis=0).drop(glia, axis=1)
        
        node_list = list(matrix.columns)
        
        # All + postemb.
        save_to_mat(matrix.to_numpy(), f'A-all-postemb_adj_mat_{dataset_i}')

        # All.
        postemb = [n for n in node_list if is_postemb(n)]
        matrix[postemb] = 0
        matrix.loc[postemb] = 0
        save_to_mat(matrix.to_numpy(), f'B-all-{dataset_i}')
        
        # No variable.
        
        
        print(matrix.sum().sum())
        
        break
    
    

#    nws = [loader.remove_nodes(nw, ['other']) for nw in nws]
#    node_list = generate_node_int_map(os.path.join(DATA_DIR, 'node_int.csv'))
#    
#    print([item for item, count in collections.Counter(node_list).items() if count > 1])
#
#
#    generate_pair_info_file(node_list, os.path.join(DATA_DIR, 'node_pairs.csv'))
#    
#    # All + postemb
#    write_edge_lists(nws, node_list, 'A-all-postemb')
#    
#    # All
#    nws = [loader.remove_postemb(nw) for nw in nws]
#    write_edge_lists(nws, node_list, 'B-all')
#    
#    # No variable.
#    nws = [loader.remove_edges(nw, type_of_edges_to_remove=['remainder']) for nw in nws]
#    write_edge_lists(nws, node_list, 'D-no-variable')
#    
#    # Only stable.
#    nws = [loader.remove_edges(nw, type_of_edges_to_remove=['increase', 'decrease']) for nw in nws]
#    write_edge_lists(nws, node_list, 'E-stable')
    
    