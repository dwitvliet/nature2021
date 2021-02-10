import os
import random
import itertools

import numpy as np
import pandas as pd
import scipy.io

from src.data.neuron_info import ntype, is_postemb


from src.data import data_manager


def export_to_matlab(path, edge_classifications):

    #DATA_DIR = '../WSBM/data/'
    G = data_manager.get_connections()['count'].copy()
    classifications = edge_classifications.copy()
    random.seed(0)

    get_matrix = lambda d: data_manager.series_to_dense_matrix(G[d]).sort_index().sort_index(axis=1)
    save_to_mat = lambda arr, name: scipy.io.savemat(os.path.join(path, f'{name}.mat'), {'arr': arr})

    l1_matrix = get_matrix('Dataset1').to_numpy().copy()
    adult_matrix = get_matrix('Dataset8').to_numpy().copy()

    l1_synapses = l1_matrix.sum()
    l1_edges = l1_matrix.astype(bool).sum()
    adult_synapses = adult_matrix.sum()
    adult_edges = adult_matrix.astype(bool).sum()


    for dataset_i, dataset in enumerate(G.columns):

        matrix = get_matrix(dataset)

        # Drop glia
        glia = [n for n in matrix.columns if ntype(n) == 'other']
        matrix = matrix.drop(glia, axis=0).drop(glia, axis=1)

        node_list = list(matrix.columns)

        # All + postemb.
        save_to_mat(matrix.to_numpy(), f'A-all-postemb_{dataset_i}')

        # All.
        postemb = [n for n in node_list if is_postemb(n)]
        matrix[postemb] = 0
        matrix.loc[postemb] = 0
        save_to_mat(matrix.to_numpy(), f'B-all_{dataset_i}')

        # Prune synapses until L1 numbers.
        matrix_pruned = matrix.to_numpy().copy()
        while matrix_pruned.astype(bool).sum() > l1_edges:
            edge_to_prune = random.choice(list(zip(*np.where(matrix_pruned > 0))))
            matrix_pruned[edge_to_prune] -= 1
        while matrix_pruned.sum() > l1_synapses:
            edge_to_prune = random.choice(list(zip(*np.where(matrix_pruned > 1))))
            matrix_pruned[edge_to_prune] -= 1
        save_to_mat(matrix_pruned, f'E-all-pruned-to-l1_{dataset_i}')

        # Randomly add synapses until adult numbers.
        matrix_grown_random = matrix.to_numpy().copy()
        all_edges = list(itertools.product(
            range(matrix_grown_random.shape[1]),
            range(matrix_grown_random.shape[0])
        ))
        while matrix_grown_random.astype(bool).sum() < adult_edges:
            edge_to_grow = random.choice(all_edges)
            matrix_grown_random[edge_to_grow] += 1
        nonzero_edges = list(zip(*np.where(matrix_grown_random > 0)))
        while matrix_grown_random.sum() < adult_synapses:
            edge_to_grow = random.choice(nonzero_edges)
            matrix_grown_random[edge_to_grow] += 1
        save_to_mat(matrix_grown_random, f'F-all-grown-to-adult-random_{dataset_i}')

        # Add synapses to existing cells until adult adult numbers.
        matrix_grown = matrix.to_numpy().copy()
        pres = range(matrix_grown.shape[1])
        posts = range(matrix_grown.shape[0])
        pre_synapses = matrix_grown.sum(axis=0)
        post_synapses = matrix_grown.sum(axis=1)
        while matrix_grown.astype(bool).sum() < adult_edges:
            pre_to_grow = random.choices(pres, weights=pre_synapses)[0]
            post_to_grow = random.choices(posts, weights=post_synapses)[0]
            matrix_grown[post_to_grow, pre_to_grow] += 1
        nonzero_edges = list(zip(*np.where(matrix_grown > 0)))
        nonzero_edge_weights = [pre_synapses[n1] + post_synapses[n2] for n1, n2 in nonzero_edges]
        while matrix_grown.sum() < adult_synapses:
            edge_to_grow = random.choices(nonzero_edges, nonzero_edge_weights)[0]
            matrix_grown[edge_to_grow] += 1
        save_to_mat(matrix_grown, f'G-all-grown-to-adult_{dataset_i}')

        # No variable.
        variable_edges = classifications[classifications.isin(['noise', 'remainder'])].index
        for pre, post in variable_edges:
            if pre not in node_list or post not in node_list:
                continue
            matrix.loc[post, pre] = 0
        save_to_mat(matrix.to_numpy(), f'C-no-variable_{dataset_i}')

        # Only stable.
        dynamic_edges = classifications[classifications.isin(['increase', 'decrease'])].index
        for pre, post in dynamic_edges:
            if pre not in node_list or post not in node_list:
                continue
            matrix.loc[post, pre] = 0
        save_to_mat(matrix.to_numpy(), f'D-only_stable_{dataset_i}')


    # Save cell list.
    pd.Series(range(1, len(node_list)+1), index=node_list).to_csv(os.path.join(path, 'node_int.csv'), header=False)

    # # Save list of left-right pairs (why?)
    # pairs = defaultdict(list)
    # for i, n in enumerate(node_list, 1):
    #     pair = npair(n)
    #     if pair == n:
    #         continue
    #     pairs[pair].append(i)
    # with open(os.path.join(path, 'node_pairs.csv'), 'w') as f:
    #     for pair in pairs.values():
    #         f.write(f'{pair[0]},{pair[1]}\n')


