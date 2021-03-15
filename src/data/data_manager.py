import os
import json
import warnings
from functools import lru_cache, wraps
from collections import defaultdict

import numpy as np
import pandas as pd

from src.data.dataset_info import all_datasets, datasets_with_adj, upscale_l3
from src.data.neuron_info import neuron_list, in_brain, class_members, \
    ntype, npair, nclass, is_postemb, contralateral

synapse_data = 'data/synapses'
skeleton_data = 'data/skeletons'
volumetric_data = 'data/physical_contact'
legacy_data = 'data/legacy_data'

valid_cells = set([n for c in neuron_list if in_brain(c) for n in class_members(c)])


def list_to_tuple(func):
    """ Decorator to allow lru_cache to cache list arguments """

    @wraps(func)
    def wrapper(*args, **kwargs):
        args = [tuple(x) if type(x) == list else x for x in args]
        kwargs = {k: tuple(v) if type(v) == list else v for k, v in
                  kwargs.items()}
        return func(*args, **kwargs)

    return wrapper


@list_to_tuple
@lru_cache(maxsize=None)
def get_gapjunctions(datasets=all_datasets):
    synapses = []
    for dataset in datasets:
        with open(os.path.join(synapse_data, dataset + '_gapjunctions.json')) as f:
            dataset_synapses = json.load(f)
        for synapse in dataset_synapses:
            if synapse['n1'] not in valid_cells or synapse['n1'] not in valid_cells:
                continue
            synapse['dataset'] = dataset
            synapses.append(synapse)

    return pd.DataFrame(synapses)


@list_to_tuple
@lru_cache(maxsize=None)
def get_gapjunction_connections(datasets=all_datasets):
    df = get_gapjunctions(datasets=datasets).copy()
    df['count'] = 1
    df = df.groupby(['dataset', 'n1', 'n2'], as_index=False)['count'].sum() \
        .set_index(['n1', 'n2']) \
        .pivot(columns='dataset')

    df = df['count']
    df = df.fillna(0).astype(int)

    return df


@list_to_tuple
@lru_cache(maxsize=None)
def get_synapses(datasets=all_datasets):
    """ Get synapses before splitting into one-to-one.

    Args:
        datasets (tuple of str, optional): Datasets to get synapses for.

    Returns:
        dict: Each key a dataset, and each value a list of synapses.

    """
    synapses = {d: [] for d in datasets}
    for dataset in datasets:

        fpath = os.path.join(synapse_data, dataset + '_synapses.json')
        if dataset == 'white_adult':
            fpath = os.path.join(legacy_data, 'catmaid_N2U_nmjs.json')
        if dataset == 'white_l4':
            fpath = os.path.join(legacy_data, 'catmaid_JSH_nmjs.json')

        with open(fpath) as f:
            dataset_synapses = json.load(f)

        # Filter out synapses between invalid cells.
        for synapse in dataset_synapses:
            if synapse['pre'] not in valid_cells:
                continue

            valid_posts = [n for n in synapse['post'] if n in valid_cells]
            if len(valid_posts) == 0:
                continue

            if 'post_weights' in synapse:
                valid_post_weigths = []
                for valid_post in valid_posts:
                    post_idx = synapse['post'].index(valid_post)
                    valid_post_weigths.append(synapse['post_weights'][post_idx])
                synapse['post_weights'] = valid_post_weigths

            synapse['post'] = valid_posts

            synapses[dataset].append(synapse)

    return synapses


@list_to_tuple
@lru_cache(maxsize=None)
def get_synapses_one_to_one(datasets=all_datasets):
    synapses_one_to_one = []

    synapses = get_synapses(datasets=datasets)
    for dataset in datasets:
        for i, raw_synapse in enumerate(synapses[dataset]):
            pre = raw_synapse['pre']
            for post_idx, (post, post_tid) in enumerate(
                    zip(raw_synapse['post'], raw_synapse['post_tid'])):
                if pre not in valid_cells or post not in valid_cells:
                    continue

                synapse = {
                    'dataset': dataset,
                    'pre': pre, 'post': post,
                    'pre_tid': raw_synapse['pre_tid'], 'post_tid': post_tid,
                    'syn_id': raw_synapse['catmaid_id'][0] if type(
                        raw_synapse['catmaid_id']) == list else raw_synapse[
                        'catmaid_id']
                }

                if 'size' in raw_synapse:
                    synapse['num_post'] = len(raw_synapse['post'])
                    synapse['post_weight'] = raw_synapse['post_weights'][
                        post_idx]
                    synapse['sections'] = len(
                        raw_synapse['size_proportion_per_section'])
                    synapse['size'] = raw_synapse['size'] * synapse[
                        'post_weight']

                synapses_one_to_one.append(synapse)

    synapses_one_to_one = pd.DataFrame(synapses_one_to_one)

    #    if 'size' in synapses_one_to_one:
    #        postemb_classes = [n for n in neuron_list if is_postemb(n)]
    #        pre_postemb = synapses_one_to_one['pre'].map(nclass).isin(postemb_classes)
    #        post_postemb = synapses_one_to_one['post'].map(nclass).isin(postemb_classes)
    #        sums = synapses_one_to_one[~pre_postemb & ~post_postemb].groupby('dataset')['size'].sum()
    #        for d in sums.index:
    #            s = synapses_one_to_one
    #            s.loc[s['dataset'] == d, 'size_norm'] = s[s['dataset'] == d]['size'] / sums[d] * sums.mean()

    return synapses_one_to_one


@list_to_tuple
@lru_cache(maxsize=None)
def get_connections(datasets=all_datasets, synapses=None):

    if synapses is None:
        synapses = get_synapses_one_to_one(datasets=datasets)

    synapses['count'] = 1
    connections = synapses \
        .drop(['post_weight', 'syn_id'], axis=1, errors='ignore') \
        .groupby(['pre', 'post', 'dataset'], as_index=False).sum() \
        .set_index(['pre', 'post']) \
        .pivot(columns='dataset')

    connections['count']  # without this line, fillna() does nothing.
    connections = connections.fillna(0).astype(int)

    if 'Dataset7' in datasets:
        for col in ['sections', 'size']:
            connections[(col, 'Dataset7')] = np.nan

    # Sort columns by dataset age.
    cols = sorted(connections.columns)
    connections = connections[cols]

    return connections


@list_to_tuple
@lru_cache(maxsize=None)
def get_legacy_connections(dataset_to_get):
    # Load Durbin datasets.
    # File is directly from wormatlas.
    # - No duplicate connections are present.
    # - PVT and DVB are switched.
    # - There are 48 connections between neurons that are inconsistent
    #   when comparing Send(_joint) with Receive(_joint) and when
    #   comparing gap junctions. These connections were either missed in
    #   one direction (35) or deviate by one or two synapses (13). Taking
    #   the max should resolve both issues.

    connections = {}
    with open(os.path.join(legacy_data, 'durbin.txt')) as f:
        for line in f:
            pre, post, typ, dataset, synapses = line.strip().split('\t')

            # Skip muscles as they were reannotated.
            if post == 'mu_bod':
                continue

            synapses = int(synapses)
            pre = pre.replace('DVB', 'PVT')
            post = post.replace('DVB', 'PVT')
            if pre not in valid_cells or post not in valid_cells:
                continue

            rev_type = {'Gap_junction': 'Gap_junction',
                        'Send': 'Receive',
                        'Receive': 'Send',
                        'Send_joint': 'Receive_joint',
                        'Receive_joint': 'Send_joint'}

            if dataset_to_get != {'N2U': 'white_adult', 'JSH': 'white_l4'}[dataset]:
                continue

            key = (typ, pre, post)
            rev_key = (rev_type[typ], post, pre)

            connections[key] = max(connections.get(key, 0), synapses)
            connections[rev_key] = max(connections.get(rev_key, 0), synapses)

    connections_joined = defaultdict(int)
    for (typ, pre, post), synapses in connections.items():
        if typ in ('Receive', 'Receive_joint'):
            continue
        if typ == 'Gap_junction':
            continue
        connections_joined[(pre, post)] += synapses

    legacy_connections = pd.Series(connections_joined, name=dataset_to_get)

    catmaid_connections = get_connections([dataset_to_get])[('count', dataset_to_get)]

    legacy_connections.index = legacy_connections.index.set_names(catmaid_connections.index.names)

    return pd.concat([legacy_connections, catmaid_connections]).sort_index()


@lru_cache(maxsize=None)
def get_cook_data():
    edges = defaultdict(int)

    with open(os.path.join(legacy_data, 'wormwiring_N2U.txt')) as f:
        for line in f:
            pre, post, syn = line.split('\t')

            # Clean up labels.
            if pre.endswith('.'):
                pre = pre[:-1]
            if post.endswith('.'):
                post = post[:-1]
            if pre.endswith('_old'):
                pre = pre[:-4]
            if post.endswith('_old'):
                post = post[:-4]
            if post == '[CANL]':
                post = 'CANL'

            if post in ('obj25484', 'exc_gl', 'hyp', 'seamL', 'pharyngealepithelium'):
                continue

            # Remove non-brain.
            if {pre, post}.intersection({
                'SABD', 'AS01', 'VB01', 'DA01', 'DD01', 'DB01', 'VC01', 'VB02', 'VA02'
            }):
                continue

            # Remove fragments.
            if {pre, post}.intersection({
                'DROP-unk', 'unk', 'd-blue-25', '[SMDD.R]_incorrect',
                'dBWMLrename',
                'BWM-DL0r', 'obj25640', 'obj8956', 'obj8957', 'obj25486',
                'obj51960',
                'obj25485', 'obj25487', 'obj41106'
            }):
                continue

            # Fix incorrect muscle assignments.
            if 'BWM' in post:
                post = 'BWM-' + post[0].upper() + post[4] + '0' + post[5]

                muscle_map = {
                    'BWM-VL01': 'BWM-VL02',
                    'BWM-VR01': 'BWM-VR02',
                    'BWM-VL02': 'BWM-VL01',
                    'BWM-VR02': 'BWM-VR01',
                    'BWM-DL02': 'BWM-DL03',
                    'BWM-DR02': 'BWM-DR03',
                    'BWM-DL03': 'BWM-DL02',
                    'BWM-DR03': 'BWM-DR02',
                    'BWM-DL05': 'BWM-DL06',
                    'BWM-DR05': 'BWM-DR06',
                    'BWM-VL05': 'BWM-VL06',
                    'BWM-VR05': 'BWM-VR06',
                    'BWM-DL06': 'BWM-DL05',
                    'BWM-DR06': 'BWM-DR05',
                    'BWM-VL06': 'BWM-VL05',
                    'BWM-VR06': 'BWM-VR05'
                }

                if post in muscle_map:
                    post = muscle_map[post]

            if post == 'BWM-DL07':
                # only traced outside the nerve ring.
                continue

            if pre not in valid_cells:
                print('pre', pre)
            if post not in valid_cells:
                print('pos', post)

            edges[(pre, post)] += int(syn)

    return pd.Series(edges)


@list_to_tuple
@lru_cache(maxsize=None)
def get_adjacency(datasets=all_datasets):
    adjacency = pd.DataFrame(
        index=pd.MultiIndex(levels=[[], []], codes=[[], []],
                            names=['pre', 'post']))

    for dataset in datasets:
        if dataset not in datasets_with_adj:
            adjacency[dataset] = np.nan
            continue

        file_path = os.path.join(volumetric_data, dataset + '_adjacency.csv')

        dataset_adj = pd.read_csv(file_path, index_col=0)
        invalid_columns = [c for c in dataset_adj.columns if
                           c not in valid_cells or ntype(c) == 'other']
        dataset_adj = dataset_adj \
            .drop(invalid_columns, axis=0) \
            .drop(invalid_columns, axis=1) \
            .stack()
        dataset_adj = dataset_adj[dataset_adj > 0]
        dataset_adj.index = dataset_adj.index.set_names(('pre', 'post'))
        dataset_adj.name = dataset
        adjacency = adjacency.merge(dataset_adj, 'outer', left_index=True,
                                    right_index=True, copy=False)

    for dataset in datasets:
        if dataset in datasets_with_adj:
            adjacency[dataset] = adjacency[dataset].fillna(0).astype(int)

    return adjacency


@list_to_tuple
@lru_cache(maxsize=None)
def get_volumes(datasets=all_datasets):
    volumes = pd.DataFrame()
    volumes.index.name = 'neuron'

    for dataset in datasets:
        if dataset not in datasets_with_adj:
            volumes[dataset] = np.nan
            continue

        file_path = os.path.join(volumetric_data, dataset + '_volumes.csv')
        dataset_vol = pd.read_csv(file_path, index_col=0, squeeze=True)

        invalid_cells = [c for c in dataset_vol.index if
                         c not in valid_cells or ntype(c) == 'other']
        dataset_vol = dataset_vol.drop(invalid_cells)

        dataset_vol.name = dataset
        volumes = volumes.merge(dataset_vol, 'outer', left_index=True,
                                right_index=True, copy=False)

    for dataset in datasets:
        if dataset in datasets_with_adj:
            volumes[dataset] = volumes[dataset].fillna(0).astype(int)

    return volumes


@list_to_tuple
@lru_cache(maxsize=None)
def get_skeletons(datasets=all_datasets):
    confirmed_zero_length = {
        'Dataset1': ['ALNL', 'ALNR', 'AQR', 'AVFL', 'AVFR', 'AVM', 'HSNL',
                     'HSNR', 'PLNL', 'PLNR', 'PVNL', 'PVNR', 'RMFL', 'RMFR',
                     'RMHL', 'RMHR', 'SDQL', 'SDQR'],
        'Dataset2': ['ALNL', 'ALNR', 'AQR', 'AVFL', 'AVFR', 'AVM', 'HSNL',
                     'HSNR', 'PLNL', 'PLNR', 'PVNL', 'PVNR', 'RMFL', 'RMFR',
                     'RMHL', 'RMHR', 'SDQL', 'SDQR'],
        'Dataset3': ['ALNL', 'ALNR', 'AQR', 'AVFL', 'AVFR', 'AVM', 'HSNL',
                     'HSNR', 'PLNL', 'PLNR', 'PVNL', 'PVNR', 'RMFL', 'RMFR',
                     'RMHL', 'RMHR', 'SDQL', 'SDQR'],
        'Dataset4': ['AVFL', 'AVFR', 'AVM', 'HSNL', 'HSNR', 'PLNL', 'PLNR',
                     'PVNL', 'PVNR', 'RMFL', 'RMFR', 'SDQL'],
        'Dataset5': ['HSNL', 'HSNR', 'PLNL', 'PLNR', 'PVNL', 'PVNR', 'RMFL',
                     'RMFR'],
        'Dataset6': ['HSNL', 'HSNR', 'PLNL', 'PLNR', 'PVNL', 'PVNR'],
        'Dataset7': [],
        'Dataset8': []
    }
    confirmed_zero_length['Dataset1'].append('BWM-DL07')

    skeletons = {}

    for dataset in datasets:
        with open(os.path.join(skeleton_data, dataset + '_skeletons.json')) as f:
            dataset_skeletons = json.load(f)

        dataset_skeletons = {n: s for n, s in dataset_skeletons.items() if
                             n in valid_cells and ntype(n) != 'other'}

        # Convert JSON keys to int.
        for n, skeleton in dataset_skeletons.items():
            for prop in ['arbor', 'coords', 'dist_to_root']:
                skeleton[prop] = {int(k): v for k, v in skeleton[prop].items()}

        # Apply transform for Dataset6 due to shrinkage.
        if dataset == 'Dataset6':
            for n, skeleton in dataset_skeletons.items():
                skeleton['length'] *= upscale_l3
                for tid in skeleton['dist_to_root']:
                    skeleton['dist_to_root'][tid] *= upscale_l3

        # Check that no cells are missing from export.
        for n in valid_cells:
            if ntype(n) == 'other':
                continue
            if n in confirmed_zero_length[dataset]:
                continue
            if n not in dataset_skeletons or dataset_skeletons[n]['length'] == 0:
                print(f'Error: {n} in {dataset} has a length of 0')

        skeletons[dataset] = dataset_skeletons

    return skeletons


def to_npair(df):
    warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
    df = df.copy()
    index_names = df.index.names
    temp_index_names = []
    for i, idx_name in enumerate(df.index.names):
        temp_index_name = idx_name + '_pair'
        temp_index_names.append(temp_index_name)
        df[temp_index_name] = df.index.get_level_values(i).map(npair)

    df_pair = df.groupby(temp_index_names).sum()
    df_pair.index.names = index_names
    return df_pair


def remove_postemb(df):
    df = df.copy()
    postemb = [n for n in neuron_list if is_postemb(n)]
    n1_postemb = df.index.get_level_values(0).map(nclass).isin(postemb)
    n2_postemb = df.index.get_level_values(1).map(nclass).isin(postemb)
    return df[~n1_postemb & ~n2_postemb]


def series_to_dense_matrix(series):
    valid_neurons = set(
        [n for c in neuron_list if in_brain(c) for n in class_members(c)])

    matrix = series \
        .reset_index() \
        .pivot(columns='pre', index='post', values=series.name) \
        .fillna(0) \
        .astype(int)

    missing_pre = valid_neurons - set(matrix.columns)
    missing_post = valid_neurons - set(matrix.index)
    invalid_pre = set(matrix.columns) - valid_neurons
    invalid_post = set(matrix.index) - valid_neurons

    for n in missing_pre:
        matrix[n] = 0
    for n in missing_post:
        matrix.loc[n] = 0
    for n in invalid_pre:
        matrix = matrix.drop(n, axis=1)
    for n in invalid_post:
        matrix = matrix.drop(n, axis=0)

    return matrix


def remove_noise(df):
    """
        Only use edges with 2 or more synapses in at least two timepoints.
        Or at least one timepoint has 2 or more synapses and the contralateral edge exists.
    """

    assert type(df.columns) == pd.core.indexes.base.Index, 'The colums should be single index'

    df = df.copy()

    noise_edges = []
    for pre, post in df.index:
        edge = pre, post
        edge_c = (contralateral(pre), contralateral(post))
        syns = df.loc[edge]
        if edge != edge_c and edge_c in df.index:
            syns_c = df.loc[edge_c]
        else:
            syns_c = np.zeros(df.shape[1])

        if (syns >= 2).sum() >= 2:
            continue
        if any((syns >= 1) & (syns_c >= 2)):
            continue
        if any((syns >= 2) & (syns_c >= 1)):
            continue

        noise_edges.append(edge)

    return df.drop(noise_edges)


def get_communities():
    communities = {
        'Dataset1': [
            ['ADAL', 'ADAR', 'ADFL', 'ADFR', 'ADLR', 'AIAL', 'AIAR', 'AIBL', 'AIBR', 'AIYL', 'AIYR', 'AIZL', 'AIZR', 'ASER', 'ASKR', 'AUAL', 'AUAR', 'AVAL', 'AVAR', 'AVBL', 'AVBR', 'AVDR', 'AVEL', 'AVER', 'AVHR', 'AVKL', 'AVKR', 'AWCL', 'AWCR', 'BAGL', 'BAGR', 'BWM-DL01', 'BWM-DL02', 'BWM-DR01', 'BWM-DR02', 'BWM-VL01', 'BWM-VL02', 'BWM-VR01', 'BWM-VR02', 'CEPDL', 'CEPVL', 'CEPVR', 'DVA', 'IL1L', 'IL1R', 'IL2L', 'IL2R', 'OLLL', 'OLQVL', 'OLQVR', 'PVR', 'RIAL', 'RIAR', 'RIBL', 'RIBR', 'RICL', 'RICR', 'RIFL', 'RIGL', 'RIGR', 'RIH', 'RIML', 'RIMR', 'RIPL', 'RIPR', 'RIR', 'RIS', 'RMDDL', 'RMDDR', 'RMDL', 'RMDR', 'RMDVL', 'RMDVR', 'RMED', 'RMEL', 'RMER', 'RMEV', 'SAADR', 'SMBVL', 'SMBVR', 'SMDDL', 'SMDDR', 'SMDVL', 'SMDVR', 'URXL', 'URXR', 'URYVR'],
            ['ADEL', 'ADER', 'ADLL', 'AFDL', 'AFDR', 'AIML', 'AIMR', 'AINL', 'AINR', 'ALA', 'ALML', 'ALMR', 'ASEL', 'ASGL', 'ASGR', 'ASHL', 'ASHR', 'ASIL', 'ASIR', 'ASJL', 'ASJR', 'ASKL', 'AVDL', 'AVHL', 'AVJL', 'AVJR', 'AVL', 'AWAL', 'AWAR', 'AWBL', 'AWBR', 'BDUL', 'BDUR', 'BWM-DL03', 'BWM-DL04', 'BWM-DL05', 'BWM-DL06', 'BWM-DL07', 'BWM-DL08', 'BWM-DR03', 'BWM-DR04', 'BWM-DR05', 'BWM-DR06', 'BWM-DR07', 'BWM-DR08', 'BWM-VL03', 'BWM-VL04', 'BWM-VL05', 'BWM-VL06', 'BWM-VL07', 'BWM-VL08', 'BWM-VR03', 'BWM-VR04', 'BWM-VR05', 'BWM-VR06', 'BWM-VR07', 'BWM-VR08', 'CEPDR', 'DVC', 'FLPL', 'FLPR', 'IL1DL', 'IL1DR', 'IL1VL', 'IL1VR', 'IL2DL', 'IL2DR', 'IL2VL', 'IL2VR', 'OLLR', 'OLQDL', 'OLQDR', 'PVCL', 'PVCR', 'PVPL', 'PVPR', 'PVQL', 'PVQR', 'PVT', 'RID', 'RIFR', 'RIVL', 'RIVR', 'RMGL', 'RMGR', 'SAADL', 'SAAVL', 'SAAVR', 'SIADL', 'SIADR', 'SIAVL', 'SIAVR', 'SIBDL', 'SIBDR', 'SIBVL', 'SIBVR', 'SMBDL', 'SMBDR', 'URADL', 'URADR', 'URAVL', 'URAVR', 'URBL', 'URBR', 'URYDL', 'URYDR', 'URYVL'],
            [],
            [],
            [],
            [],
            [],
            [],
        ],
        'Dataset2': [
            [],
            [],
            ['AFDL', 'AFDR', 'ASIL', 'ASIR', 'AVDR', 'AVL', 'BWM-DL01', 'BWM-DL02', 'BWM-DL03', 'BWM-DL04', 'BWM-DL05', 'BWM-DL06', 'BWM-DL07', 'BWM-DL08', 'BWM-DR01', 'BWM-DR02', 'BWM-DR03', 'BWM-DR04', 'BWM-DR05', 'BWM-DR06', 'BWM-DR07', 'BWM-DR08', 'BWM-VL01', 'BWM-VL02', 'BWM-VL03', 'BWM-VL04', 'BWM-VL05', 'BWM-VL06', 'BWM-VL07', 'BWM-VL08', 'BWM-VR01', 'BWM-VR02', 'BWM-VR03', 'BWM-VR04', 'BWM-VR05', 'BWM-VR06', 'BWM-VR07', 'BWM-VR08', 'RID', 'SIADL', 'SIADR', 'SIAVL', 'SIAVR', 'SIBDL', 'SIBDR', 'SIBVL', 'SIBVR'],
            [],
            ['AIAL', 'AIBL', 'AIBR', 'AIZL', 'AIZR', 'AVAL', 'AVAR', 'AVBL', 'AVBR', 'AVDL', 'AVEL', 'AVER', 'IL1L', 'RIAL', 'RIAR', 'RIBL', 'RIBR', 'RICL', 'RICR', 'RIML', 'RIMR', 'RIPL', 'RIPR', 'RMDDL', 'RMDDR', 'RMDL', 'RMDR', 'RMDVL', 'RMDVR', 'RMED', 'RMEL', 'RMER', 'RMEV', 'RMGL', 'SAAVL', 'SAAVR', 'SMBDL', 'SMBDR', 'SMBVL', 'SMBVR', 'SMDDL', 'SMDDR', 'SMDVL', 'SMDVR'],
            [],
            ['ADAL', 'ADAR', 'ADEL', 'ADER', 'ADFL', 'ADFR', 'ADLL', 'ADLR', 'AIAR', 'AIML', 'AIMR', 'AINL', 'AINR', 'AIYL', 'AIYR', 'ALA', 'ALML', 'ALMR', 'ASEL', 'ASER', 'ASGL', 'ASGR', 'ASHL', 'ASHR', 'ASJL', 'ASJR', 'ASKL', 'ASKR', 'AUAL', 'AUAR', 'AVHL', 'AVHR', 'AVJL', 'AVJR', 'AVKL', 'AVKR', 'AWAL', 'AWAR', 'AWBL', 'AWBR', 'AWCL', 'AWCR', 'BAGL', 'BAGR', 'BDUL', 'BDUR', 'CEPDL', 'CEPDR', 'CEPVL', 'CEPVR', 'DVA', 'DVC', 'FLPL', 'FLPR', 'IL1DL', 'IL1DR', 'IL1R', 'IL1VL', 'IL1VR', 'IL2DL', 'IL2DR', 'IL2L', 'IL2R', 'IL2VL', 'IL2VR', 'OLLL', 'OLLR', 'OLQDL', 'OLQDR', 'OLQVL', 'OLQVR', 'PVCL', 'PVCR', 'PVPL', 'PVPR', 'PVQL', 'PVQR', 'PVR', 'PVT', 'RIFL', 'RIFR', 'RIGL', 'RIGR', 'RIH', 'RIR', 'RIS', 'RIVL', 'RIVR', 'RMGR', 'SAADL', 'SAADR', 'URADL', 'URADR', 'URAVL', 'URAVR', 'URBL', 'URBR', 'URXL', 'URXR', 'URYDL', 'URYDR', 'URYVL', 'URYVR'],
            [],
        ],
        'Dataset3': [
            [],
            [],
            ['AVDL', 'AVDR', 'BWM-DL01', 'BWM-DL02', 'BWM-DL03', 'BWM-DL04', 'BWM-DL05', 'BWM-DL06', 'BWM-DL07', 'BWM-DL08', 'BWM-DR01', 'BWM-DR02', 'BWM-DR03', 'BWM-DR04', 'BWM-DR05', 'BWM-DR06', 'BWM-DR07', 'BWM-DR08', 'BWM-VL01', 'BWM-VL02', 'BWM-VL03', 'BWM-VL04', 'BWM-VL05', 'BWM-VL06', 'BWM-VL07', 'BWM-VL08', 'BWM-VR01', 'BWM-VR02', 'BWM-VR03', 'BWM-VR04', 'BWM-VR05', 'BWM-VR06', 'BWM-VR07', 'BWM-VR08', 'RID', 'SIADL', 'SIADR', 'SIAVL', 'SIAVR', 'SIBDL', 'SIBDR', 'SIBVL', 'SIBVR'],
            [],
            ['ADAL', 'ADAR', 'AIAL', 'AIAR', 'AIBL', 'AIBR', 'AIYL', 'AIYR', 'AIZL', 'AIZR', 'ASER', 'AVAL', 'AVAR', 'AVBL', 'AVBR', 'AVEL', 'AVER', 'AVJR', 'AWCL', 'BAGR', 'CEPDL', 'OLLL', 'OLLR', 'RIAL', 'RIAR', 'RIBL', 'RIBR', 'RICL', 'RICR', 'RIGR', 'RIH', 'RIML', 'RIMR', 'RIPL', 'RIPR', 'RIR', 'RIS', 'RMDDL', 'RMDDR', 'RMDL', 'RMDR', 'RMDVL', 'RMDVR', 'RMED', 'RMEL', 'RMER', 'RMEV', 'SAADR', 'SMBDL', 'SMBVL', 'SMBVR', 'SMDDL', 'SMDDR', 'SMDVL', 'SMDVR'],
            [],
            ['ADEL', 'ADER', 'ADFL', 'ADFR', 'ADLL', 'ADLR', 'AFDL', 'AFDR', 'AIML', 'AIMR', 'AINL', 'AINR', 'ALA', 'ALML', 'ALMR', 'ASEL', 'ASGL', 'ASGR', 'ASHL', 'ASHR', 'ASIL', 'ASIR', 'ASJL', 'ASJR', 'ASKL', 'ASKR', 'AUAL', 'AUAR', 'AVHL', 'AVHR', 'AVJL', 'AVKL', 'AVKR', 'AVL', 'AWAL', 'AWAR', 'AWBL', 'AWBR', 'AWCR', 'BAGL', 'BDUL', 'BDUR', 'CEPDR', 'CEPVL', 'CEPVR', 'DVA', 'DVC', 'FLPL', 'FLPR', 'IL1DL', 'IL1DR', 'IL1L', 'IL1R', 'IL1VL', 'IL1VR', 'IL2DL', 'IL2DR', 'IL2L', 'IL2R', 'IL2VL', 'IL2VR', 'OLQDL', 'OLQDR', 'OLQVL', 'OLQVR', 'PVCL', 'PVCR', 'PVPL', 'PVPR', 'PVQL', 'PVQR', 'PVR', 'PVT', 'RIFL', 'RIFR', 'RIGL', 'RIVL', 'RIVR', 'RMGL', 'RMGR', 'SAADL', 'SAAVL', 'SAAVR', 'SMBDR', 'URADL', 'URADR', 'URAVL', 'URAVR', 'URBL', 'URBR', 'URXL', 'URXR', 'URYDL', 'URYDR', 'URYVL', 'URYVR'],
            [],
        ],
        'Dataset4': [
            [],
            [],
            ['BWM-DL01', 'BWM-DL02', 'BWM-DL03', 'BWM-DL04', 'BWM-DL05', 'BWM-DL06', 'BWM-DL07', 'BWM-DL08', 'BWM-DR01', 'BWM-DR02', 'BWM-DR03', 'BWM-DR04', 'BWM-DR05', 'BWM-DR06', 'BWM-DR07', 'BWM-DR08', 'BWM-VL01', 'BWM-VL02', 'BWM-VL03', 'BWM-VL04', 'BWM-VL05', 'BWM-VL06', 'BWM-VL07', 'BWM-VL08', 'BWM-VR01', 'BWM-VR02', 'BWM-VR03', 'BWM-VR04', 'BWM-VR05', 'BWM-VR06', 'BWM-VR07', 'BWM-VR08', 'RID', 'RIPL', 'RIPR', 'SIADL', 'SIADR', 'SIAVL', 'SIAVR', 'SIBDL', 'SIBDR', 'SIBVL', 'SIBVR'],
            [],
            ['AIBL', 'AIBR', 'AVAL', 'AVAR', 'AVBL', 'AVBR', 'AVDR', 'AVEL', 'AVER', 'OLLL', 'RIAL', 'RIAR', 'RIBL', 'RIBR', 'RIML', 'RIMR', 'RIS', 'RMDDL', 'RMDDR', 'RMDL', 'RMDR', 'RMDVL', 'RMDVR', 'SAADL', 'SAADR', 'SMDDL', 'SMDDR', 'SMDVL', 'SMDVR'],
            [],
            ['ADAL', 'ADAR', 'ADEL', 'ADER', 'ALA', 'ALML', 'ALMR', 'AVDL', 'AVJL', 'AVJR', 'AVKL', 'AVKR', 'AVL', 'BDUL', 'BDUR', 'CEPDL', 'CEPDR', 'CEPVL', 'CEPVR', 'FLPL', 'FLPR', 'IL1DL', 'IL1DR', 'IL1L', 'IL1R', 'IL1VL', 'IL1VR', 'IL2DL', 'IL2DR', 'IL2L', 'IL2R', 'IL2VL', 'IL2VR', 'OLLR', 'OLQDL', 'OLQDR', 'OLQVL', 'OLQVR', 'PVCL', 'PVCR', 'PVR', 'RICL', 'RICR', 'RIFL', 'RIFR', 'RIH', 'RIVL', 'RIVR', 'RMED', 'RMEL', 'RMER', 'RMEV', 'RMGL', 'RMGR', 'SAAVL', 'SAAVR', 'SMBDL', 'SMBDR', 'SMBVL', 'SMBVR', 'URADL', 'URADR', 'URAVL', 'URAVR', 'URBL', 'URBR', 'URYDL', 'URYDR', 'URYVL', 'URYVR'],
            ['ADFL', 'ADFR', 'ADLL', 'ADLR', 'AFDL', 'AFDR', 'AIAL', 'AIAR', 'AIML', 'AIMR', 'AINL', 'AINR', 'AIYL', 'AIYR', 'AIZL', 'AIZR', 'ASEL', 'ASER', 'ASGL', 'ASGR', 'ASHL', 'ASHR', 'ASIL', 'ASIR', 'ASJL', 'ASJR', 'ASKL', 'ASKR', 'AUAL', 'AUAR', 'AVHL', 'AVHR', 'AWAL', 'AWAR', 'AWBL', 'AWBR', 'AWCL', 'AWCR', 'BAGL', 'BAGR', 'DVA', 'DVC', 'PVPL', 'PVPR', 'PVQL', 'PVQR', 'PVT', 'RIGL', 'RIGR', 'RIR', 'URXL', 'URXR'],
        ],
        'Dataset5': [
            [],
            [],
            ['BWM-DL01', 'BWM-DL02', 'BWM-DL03', 'BWM-DL04', 'BWM-DL05', 'BWM-DL06', 'BWM-DL07', 'BWM-DL08', 'BWM-DR01', 'BWM-DR02', 'BWM-DR03', 'BWM-DR04', 'BWM-DR05', 'BWM-DR06', 'BWM-DR07', 'BWM-DR08', 'BWM-VL01', 'BWM-VL02', 'BWM-VL03', 'BWM-VL04', 'BWM-VL05', 'BWM-VL06', 'BWM-VL07', 'BWM-VL08', 'BWM-VR01', 'BWM-VR02', 'BWM-VR03', 'BWM-VR04', 'BWM-VR05', 'BWM-VR06', 'BWM-VR07', 'BWM-VR08', 'RID', 'SIADL', 'SIADR', 'SIAVL', 'SIAVR', 'SIBDL', 'SIBDR', 'SIBVL', 'SIBVR'],
            [],
            ['AIBL', 'AIBR', 'AVAL', 'AVAR', 'AVBL', 'AVBR', 'AVEL', 'AVER', 'IL1L', 'RIAL', 'RIAR', 'RIBL', 'RIBR', 'RIML', 'RIMR', 'RIPL', 'RIPR', 'RIVR', 'RMDDL', 'RMDDR', 'RMDL', 'RMDR', 'RMDVL', 'RMDVR', 'RMED', 'RMEL', 'RMER', 'RMEV', 'SAADL', 'SAADR', 'SAAVR', 'SMBDR', 'SMBVL', 'SMBVR', 'SMDDL', 'SMDDR', 'SMDVL', 'SMDVR'],
            [],
            ['ADAL', 'ADAR', 'ADEL', 'ADER', 'AINL', 'AUAL', 'AUAR', 'AVKL', 'AVKR', 'AVL', 'BAGL', 'BAGR', 'CEPDL', 'CEPDR', 'CEPVL', 'CEPVR', 'DVA', 'DVC', 'FLPL', 'FLPR', 'IL1DL', 'IL1DR', 'IL1R', 'IL1VL', 'IL1VR', 'IL2DL', 'IL2DR', 'IL2L', 'IL2R', 'IL2VL', 'IL2VR', 'OLLL', 'OLLR', 'OLQDL', 'OLQDR', 'OLQVL', 'OLQVR', 'PVCL', 'PVPL', 'PVR', 'PVT', 'RICL', 'RICR', 'RIGL', 'RIGR', 'RIH', 'RIR', 'RIS', 'RIVL', 'RMGL', 'RMGR', 'SAAVL', 'SMBDL', 'URADL', 'URADR', 'URAVL', 'URAVR', 'URBL', 'URBR', 'URXL', 'URXR', 'URYDL', 'URYDR', 'URYVL', 'URYVR'],
            ['ADFL', 'ADFR', 'ADLL', 'ADLR', 'AFDL', 'AFDR', 'AIAL', 'AIAR', 'AIML', 'AIMR', 'AINR', 'AIYL', 'AIYR', 'AIZL', 'AIZR', 'ALA', 'ALML', 'ALMR', 'ASEL', 'ASER', 'ASGL', 'ASGR', 'ASHL', 'ASHR', 'ASIL', 'ASIR', 'ASJL', 'ASJR', 'ASKL', 'ASKR', 'AVDL', 'AVDR', 'AVHL', 'AVHR', 'AVJL', 'AVJR', 'AWAL', 'AWAR', 'AWBL', 'AWBR', 'AWCL', 'AWCR', 'BDUL', 'BDUR', 'PVCR', 'PVPR', 'PVQL', 'PVQR', 'RIFL', 'RIFR'],
        ],
        'Dataset6': [
            [],
            [],
            ['AVL', 'BWM-DL01', 'BWM-DL02', 'BWM-DL03', 'BWM-DL04', 'BWM-DL05', 'BWM-DL06', 'BWM-DL07', 'BWM-DL08', 'BWM-DR01', 'BWM-DR02', 'BWM-DR03', 'BWM-DR04', 'BWM-DR05', 'BWM-DR06', 'BWM-DR07', 'BWM-DR08', 'BWM-VL01', 'BWM-VL02', 'BWM-VL03', 'BWM-VL04', 'BWM-VL05', 'BWM-VL06', 'BWM-VL07', 'BWM-VL08', 'BWM-VR01', 'BWM-VR02', 'BWM-VR03', 'BWM-VR04', 'BWM-VR05', 'BWM-VR06', 'BWM-VR07', 'BWM-VR08', 'RID', 'RIPL', 'RIPR', 'SIADL', 'SIADR', 'SIAVL', 'SIAVR', 'SIBDL', 'SIBDR', 'SIBVL', 'SIBVR'],
            ['AVAL', 'AVAR', 'AVBL', 'AVBR', 'AVDL', 'AVDR', 'AVEL', 'AVER'],
            ['IL1R', 'RIAL', 'RIAR', 'RMDDL', 'RMDDR', 'RMDL', 'RMDR', 'RMDVL', 'RMDVR', 'RMED', 'RMEL', 'RMER', 'RMEV', 'SMDDL', 'SMDDR', 'SMDVL', 'SMDVR'],
            ['CEPDL', 'CEPDR', 'CEPVL', 'CEPVR', 'IL1DL', 'IL1DR', 'IL1L', 'IL1VL', 'IL1VR', 'IL2DL', 'IL2DR', 'IL2L', 'IL2R', 'IL2VL', 'IL2VR', 'OLLL', 'OLLR', 'OLQDL', 'OLQDR', 'OLQVL', 'OLQVR', 'PVR', 'RICL', 'RICR', 'RIH', 'RIVL', 'RIVR', 'RMGL', 'SAADL', 'SMBDL', 'SMBDR', 'SMBVL', 'SMBVR', 'URADL', 'URADR', 'URAVL', 'URAVR', 'URBL', 'URBR', 'URYDL', 'URYDR', 'URYVL', 'URYVR'],
            ['ADAL', 'ADAR', 'ADEL', 'ADER', 'ADFL', 'ADFR', 'AIBL', 'AIBR', 'AIZL', 'AIZR', 'ALMR', 'ASHL', 'ASHR', 'AUAL', 'AUAR', 'AVHR', 'AVJL', 'AVJR', 'AVKL', 'AVKR', 'BAGL', 'BAGR', 'BDUR', 'DVA', 'DVC', 'FLPL', 'FLPR', 'PVCL', 'PVCR', 'PVPL', 'PVT', 'RIBL', 'RIBR', 'RIGL', 'RIGR', 'RIML', 'RIMR', 'RIR', 'RIS', 'RMGR', 'SAADR', 'SAAVL', 'SAAVR', 'URXL', 'URXR'],
            ['ADLL', 'ADLR', 'AFDL', 'AFDR', 'AIAL', 'AIAR', 'AIML', 'AIMR', 'AINL', 'AINR', 'AIYL', 'AIYR', 'ALA', 'ALML', 'ASEL', 'ASER', 'ASGL', 'ASGR', 'ASIL', 'ASIR', 'ASJL', 'ASJR', 'ASKL', 'ASKR', 'AVHL', 'AWAL', 'AWAR', 'AWBL', 'AWBR', 'AWCL', 'AWCR', 'BDUL', 'PVPR', 'PVQL', 'PVQR', 'RIFL', 'RIFR'],
        ],
        'Dataset7': [
            [],
            [],
            ['BWM-DL01', 'BWM-DL02', 'BWM-DL03', 'BWM-DL04', 'BWM-DL05', 'BWM-DL06', 'BWM-DL07', 'BWM-DL08', 'BWM-DR01', 'BWM-DR02', 'BWM-DR03', 'BWM-DR04', 'BWM-DR05', 'BWM-DR06', 'BWM-DR07', 'BWM-DR08', 'BWM-VL01', 'BWM-VL02', 'BWM-VL03', 'BWM-VL04', 'BWM-VL05', 'BWM-VL06', 'BWM-VL07', 'BWM-VL08', 'BWM-VR01', 'BWM-VR02', 'BWM-VR03', 'BWM-VR04', 'BWM-VR05', 'BWM-VR06', 'BWM-VR07', 'BWM-VR08', 'RIPL', 'RIPR', 'SIADL', 'SIADR', 'SIAVL', 'SIAVR', 'SIBDL', 'SIBDR', 'SIBVL', 'SIBVR'],
            ['AVAL', 'AVAR', 'AVBL', 'AVBR', 'AVDL', 'AVDR', 'AVEL', 'AVER', 'RID'],
            ['IL1DL', 'IL1DR', 'IL1L', 'IL1R', 'IL1VL', 'IL1VR', 'RIAL', 'RIAR', 'RIVL', 'RIVR', 'RMDDL', 'RMDDR', 'RMDL', 'RMDR', 'RMDVL', 'RMDVR', 'RMED', 'RMEL', 'RMER', 'RMEV', 'SMBDL', 'SMBDR', 'SMBVL', 'SMBVR', 'SMDDL', 'SMDDR', 'SMDVL', 'SMDVR', 'URADL', 'URADR', 'URAVR'],
            ['ADEL', 'ADER', 'AVKL', 'AVKR', 'AVL', 'CEPDL', 'CEPDR', 'CEPVL', 'CEPVR', 'IL2DL', 'IL2DR', 'IL2L', 'IL2R', 'IL2VL', 'IL2VR', 'OLLL', 'OLLR', 'OLQDL', 'OLQDR', 'OLQVL', 'OLQVR', 'PVR', 'PVT', 'RICL', 'RICR', 'RIH', 'RIS', 'RMGL', 'RMGR', 'URAVL', 'URBL', 'URBR', 'URYDL', 'URYDR', 'URYVL', 'URYVR'],
            ['ADAL', 'ADAR', 'ADFL', 'AIBL', 'AIBR', 'AIZL', 'AIZR', 'AUAL', 'AUAR', 'AVJL', 'AVJR', 'BAGL', 'BAGR', 'DVA', 'DVC', 'FLPL', 'FLPR', 'PVCL', 'PVCR', 'PVPL', 'PVPR', 'RIBL', 'RIBR', 'RIGL', 'RIGR', 'RIML', 'RIMR', 'RIR', 'SAADL', 'SAADR', 'SAAVL', 'SAAVR', 'URXL', 'URXR'],
            ['ADFR', 'ADLL', 'ADLR', 'AFDL', 'AFDR', 'AIAL', 'AIAR', 'AIML', 'AIMR', 'AINL', 'AINR', 'AIYL', 'AIYR', 'ALA', 'ALML', 'ALMR', 'ASEL', 'ASER', 'ASGL', 'ASGR', 'ASHL', 'ASHR', 'ASIL', 'ASIR', 'ASJL', 'ASJR', 'ASKL', 'ASKR', 'AVHL', 'AVHR', 'AWAL', 'AWAR', 'AWBL', 'AWBR', 'AWCL', 'AWCR', 'BDUL', 'BDUR', 'PVQL', 'PVQR', 'RIFL', 'RIFR'],
        ],
        'Dataset8': [
            [],
            [],
            ['BWM-DL01', 'BWM-DL02', 'BWM-DL03', 'BWM-DL04', 'BWM-DL05', 'BWM-DL06', 'BWM-DL07', 'BWM-DL08', 'BWM-DR01', 'BWM-DR02', 'BWM-DR03', 'BWM-DR04', 'BWM-DR05', 'BWM-DR06', 'BWM-DR07', 'BWM-DR08', 'BWM-VL01', 'BWM-VL02', 'BWM-VL03', 'BWM-VL04', 'BWM-VL05', 'BWM-VL06', 'BWM-VL07', 'BWM-VL08', 'BWM-VR01', 'BWM-VR02', 'BWM-VR03', 'BWM-VR04', 'BWM-VR05', 'BWM-VR06', 'BWM-VR07', 'BWM-VR08', 'SIADR', 'SIAVL', 'SIAVR'],
            ['AVAL', 'AVAR', 'AVBL', 'AVBR', 'AVDL', 'AVDR', 'AVEL', 'AVER', 'RID', 'RIPL', 'RIPR', 'SIADL', 'SIBDL', 'SIBDR', 'SIBVL', 'SIBVR'],
            ['IL1DL', 'IL1DR', 'IL1L', 'IL1R', 'IL1VL', 'IL1VR', 'RIVL', 'RIVR', 'RMDDL', 'RMDDR', 'RMDL', 'RMDR', 'RMDVL', 'RMDVR', 'RMED', 'RMEL', 'RMER', 'RMEV', 'SAADL', 'SAAVL', 'SAAVR', 'SMBDL', 'SMBDR', 'SMBVL', 'SMBVR', 'SMDDL', 'SMDDR', 'SMDVL', 'SMDVR', 'URADL', 'URADR', 'URAVL', 'URAVR'],
            ['ADAL', 'ADAR', 'ADEL', 'ADER', 'AVKL', 'AVKR', 'AVL', 'CEPDL', 'CEPDR', 'CEPVL', 'CEPVR', 'DVC', 'FLPL', 'FLPR', 'IL2DL', 'IL2DR', 'IL2L', 'IL2R', 'IL2VL', 'IL2VR', 'OLLL', 'OLLR', 'OLQDL', 'OLQDR', 'OLQVL', 'OLQVR', 'PVR', 'PVT', 'RICL', 'RICR', 'RIH', 'RIS', 'RMGL', 'RMGR', 'SAADR', 'URBL', 'URBR', 'URXR', 'URYDL', 'URYDR', 'URYVL', 'URYVR'],
            ['AIBL', 'AIBR', 'AIZL', 'AIZR', 'DVA', 'RIAL', 'RIAR', 'RIBL', 'RIBR', 'RIGL', 'RIGR', 'RIML', 'RIMR'],
            ['ADFL', 'ADFR', 'ADLL', 'ADLR', 'AFDL', 'AFDR', 'AIAL', 'AIAR', 'AIML', 'AIMR', 'AINL', 'AINR', 'AIYL', 'AIYR', 'ALA', 'ALML', 'ALMR', 'ASEL', 'ASER', 'ASGL', 'ASGR', 'ASHL', 'ASHR', 'ASIL', 'ASIR', 'ASJL', 'ASJR', 'ASKL', 'ASKR', 'AUAL', 'AUAR', 'AVHL', 'AVHR', 'AVJL', 'AVJR', 'AWAL', 'AWAR', 'AWBL', 'AWBR', 'AWCL', 'AWCR', 'BAGL', 'BAGR', 'BDUL', 'BDUR', 'PVCL', 'PVCR', 'PVPL', 'PVPR', 'PVQL', 'PVQR', 'RIFL', 'RIFR', 'RIR', 'URXL'],
        ]
    }

    community_names = ['', '', 'Motor output', 'Body movement', 'Head movement',
                       'Anterior sensory', 'Medial sensory/\ninterneuron',
                       'Posterior sensory']

    return communities, community_names
