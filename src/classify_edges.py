from itertools import compress
from functools import lru_cache

import numpy as np
import pandas as pd
import scipy as sc
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from scipy.stats import mannwhitneyu

from src.data import data_manager
from src.data.dataset_info import all_datasets, timepoint, stage
from src.data.neuron_info import npair


@lru_cache(maxsize=None)
def classify_edges(datasets=all_datasets, edge_type='count'):
    """

    Classifies connections as stable ('stable'), developmentally dynamic
    ('increase' and 'decrease'), or variable ('noise' and 'remainder').

    Args:
        datasets (tuple): datasets include for running the analysis. If data is
            not present for the selected edge_type, a dataset is automatically
            excluded. E.g. Dataset7 is excluded when edge_type is 'contact' or
            'size', as it does not have volumetric data for those measurements.
        edge_type:
            Which parameter to use for measuring changes across development.
            'count' is number of synapses, 'size' is the summed size of
            synapses, 'density' is the number of synapses over the length of the
            presynaptic neuron, 'contact' is the physical contact area.

    Returns:
        edge_classifications (pd.Series), edge_pair_classifications (pd.Series)

    """

    assert edge_type in ('count', 'size', 'density', 'contact')

    # Gather data.
    if edge_type == 'contact':
        G_postemb = data_manager.get_adjacency(datasets).copy()
    elif edge_type == 'density':
        G_postemb = data_manager.get_connections(datasets)['count'].copy()
    else:
        G_postemb = data_manager.get_connections(datasets)[edge_type].copy()

    # Remove post-embryonic cells and datasets where data is missing.
    G_postemb = G_postemb[(G_postemb.sum(axis=1) > 0)]
    G = data_manager.remove_postemb(G_postemb)

    # Normalize data so that the total sum of each dataset is comparable.
    G_normalized = G / G.sum() * G.sum().mean()

    # Filter out very weak connections.
    G_filtered = data_manager.remove_noise(G)

    # Convert synapse counts to density if requested.
    if edge_type == 'density':
        skeletons = data_manager.get_skeletons()
        G = G.astype(float)
        for dataset in G:
            for (pre, post) in G[dataset].index:
                G[dataset][(pre, post)] /= (skeletons[dataset][pre]['length'] / 1e6)
        G_normalized = G / G.sum() * G.sum().mean()
        G_filtered = G.loc[G_filtered.index]

    # Pool left-right cell pairs to get the pair connections.
    G_pair = data_manager.to_npair(G)
    G_pair_normalized = data_manager.to_npair(G_normalized)
    G_pair_filtered = data_manager.to_npair(G_filtered)
    pair_edges_filtered = G_pair_filtered[G_pair_filtered.sum(axis=1) > 0].index
    pair_edges_all = G_pair[G_pair.sum(axis=1) > 0].index

    # Calculate p-values for each pair connection.
    timepoints = [timepoint[d] for d in datasets]
    p_values = []
    for edge in pair_edges_filtered:
        synapses = G_normalized[G_normalized.index.map(lambda e: (npair(e[0]), npair(e[1])) == edge)]
        x = timepoints * synapses.shape[0]
        p_values.append(sc.stats.spearmanr(x, synapses.values.flatten()).pvalue)

    # Correct p-values for multiple comparisons by adjusting FDR.
    fdr = fdrcorrection0(p_values)
    significant_edges = list(compress(pair_edges_filtered, fdr[0]))

    # Define developmentally dynamic connections as any connection that is
    # statistically significant after correction.
    classification = {
        'increase': [],
        'decrease': [],
    }
    for edge in significant_edges:
        syns_normalized = G_pair_normalized.loc[edge]
        juv = np.mean(syns_normalized[:2])
        mat = np.mean(syns_normalized[[d for d in datasets if stage[d] == 'Adult']])
        if mat > juv*5:
            classification['increase'].append(edge)
        elif juv > mat*5:
            classification['decrease'].append(edge)

    # Define stable connections as remaining connections that are present in all
    # datasets or only missing in one dataset.
    classification['stable'] = []
    for edge in G_pair_filtered.index:
        if edge in classification['increase'] or edge in classification['decrease']:
            continue
        if np.count_nonzero(G_pair.loc[edge]) >= len(datasets) - 1:
            classification['stable'].append(edge)

    # Define variable connections the remaining connections and the initially
    # removed very weak connections.
    classification['noise'] = [
        e for e in pair_edges_all if e not in pair_edges_filtered
    ]
    classification['remainder'] = [
        e for e in pair_edges_all if e not in [h for hs in classification.values() for h in hs]
    ]

    # Save classified edges, splitting left-right pooled cells.
    classifications = {}
    pair_classifications = {e: t for t in classification for e in classification[t]}
    for (pre, post) in G.index:
        classification = pair_classifications[npair(pre), npair(post)]
        if (pre, post) not in G_filtered.index:
            classification = 'noise'
        classifications[(pre, post)] = classification

    edge_classifications = pd.Series(classifications, name='edge_classifications')
    edge_classifications.index.set_names(['pre', 'post'], inplace=True)
    edge_pair_classifications = pd.Series(pair_classifications, name='edge_classifications')
    edge_pair_classifications.index.set_names(['pre', 'post'], inplace=True)

    # Count and print classifications.
    G_classifications = G.merge(edge_classifications, left_index=True, right_index=True)
    counts = G_classifications.groupby('edge_classifications').agg(lambda s: s.astype(bool).sum())
    counts.loc['all'] = counts.sum()
    print(counts.mean(axis=1))

    # nemanode_output = {
    #     'increase': sorted([e for e, c in classifications.items() if c == 'increase']),
    #     'decrease': sorted([e for e, c in classifications.items() if c == 'decrease']),
    #     'stable': sorted([e for e, c in classifications.items() if c == 'stable']),
    #     'postembryonic': sorted([e for e in G_postemb if e not in G]),
    #     'variable': sorted([e for e, c in classifications.items() if c in ('noise', 'remainder')]),
    # }
    #
    # with open('../_data/edge_classifications.json', 'w') as f:
    #     json.dump(nemanode_output, f, indent=2)

    return edge_classifications, edge_pair_classifications
