import itertools
from collections import defaultdict

import numpy as np
import scipy as sc
import pandas as pd
from statsmodels.stats.proportion import proportions_ztest
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

from src.data import data_manager
from src.data.neuron_info import neuron_list, class_members, in_brain, is_postemb, ntype, dark_colors
from src.data.dataset_info import all_datasets, timepoint

from src.plotting import plotter


class Figure(object):

    def __init__(self, output_path, page_size=7.20472):
        self.plt = plotter.Plotter(output_path=output_path, page_size=page_size)

    def _get_branches(self, exclude_postemb=True):

        skeletons = data_manager.get_skeletons()
        cells_with_skeletons = [n for c in neuron_list for n in class_members(c) if (in_brain(c) and ntype(c) != 'other')]
        if exclude_postemb:
            cells_with_skeletons = [n for n in cells_with_skeletons if not is_postemb(n)]
        mean_length = {d: np.mean([s['length'] for n, s in skeletons[d].items() if n in cells_with_skeletons]) for d in all_datasets}

        # Separate branches into long and short.
        long_branches = defaultdict(list)
        short_branches = defaultdict(list)
        for d in all_datasets:
            for n in cells_with_skeletons:
                if n not in skeletons[d]:
                    continue
                skeleton = skeletons[d][n]
                for branch in skeleton['branches']:
                    branch_dists = [skeleton['dist_to_root'][tid] for tid in branch]
                    l = max(branch_dists) - min(branch_dists)
                    orphan_tids = [  # nerve_ring_starts and nerve_ring_ends tags on the same treenode
                        6373282, 6373249, 6373247,   # Dataset 3
                        7060635,  # Dataset6
                        8152020, 8152012,  # Dataset 7
                        8242523, 8340868,  # Dataset 8
                    ]
                    if l/mean_length[d] > 0.1 or branch[0] in orphan_tids:
                        long_branches[d].append((n, branch))
                    else:
                        short_branches[d].append((n, branch))
                        
        # Add excluded treenodes (due to 'not a branch' tags) to branches.
        for d in all_datasets:
            tid_to_branch_map = {}
            for idx, (n, branch) in enumerate(short_branches[d]):
                for tid in branch:
                    tid_to_branch_map[tid] = ('short', idx)
            for idx, (n, branch) in enumerate(long_branches[d]):
                for tid in branch:
                    tid_to_branch_map[tid] = ('long', idx)
                    
            for n in cells_with_skeletons:
                if n not in skeletons[d]:
                    continue
                arbor = skeletons[d][n]['arbor']
                for tid in arbor:
                    if tid in tid_to_branch_map:
                        continue
                    parent = arbor[tid]
                    while parent is not None:
                        if parent in tid_to_branch_map:
                            branch_type, branch_idx = tid_to_branch_map[parent]
                            if branch_type == 'short':
                                short_branches[d][branch_idx][1].append(tid)
                            if branch_type == 'long':
                                long_branches[d][branch_idx][1].append(tid)
                            break
                        parent = arbor.get(parent)
                    else:
                        print(f'Error: for {n} in {d}, {tid} was not assigned to any branch')

        return long_branches, short_branches

    def _treenodes_with_synapses(self):
        synapses = data_manager.get_synapses()
        pre_tids = {d: [] for d in all_datasets}
        post_tids = {d: [] for d in all_datasets}
        for dataset in all_datasets:
            for synapse in synapses[dataset]:
                pre_tids[dataset].append(synapse['pre_tid'])
                post_tids[dataset].extend(synapse['post_tid'])
        return pre_tids, post_tids

    def number_of_branches(self, f):
        long_branches, short_branches = self._get_branches()
        x = [timepoint[d] for d in all_datasets]
        y = (
#            ('Long branches', [len(long_branches[d]) for d in all_datasets]),
            ('Short branches', [len(short_branches[d]) for d in all_datasets]),
        )

        c = {'Long branches': 'grey', 'Long branches_edge': 'grey',
             'Short branches': 'k', 'Short branches_edge': 'k'}
        self.plt.plot(
            'xy_graph', (x, y),
            colors=c, ylim=(0, 1250), y_label='Spine-like protrusions', no_legend=True,
            yticks=range(0, 1500, 250),
            save=f+'_number_of_branches', linkpoints=False, stats='spearmanr',
            size=0.22, x_label='Developmental age',
            margin={'left': 0.05, 'right': 0.01, 'top': 0.03, 'bottom': 0.03},
        )

    def _short_branches_with_synapses(self, synapse_type='input'):
        long_branches, short_branches = self._get_branches()
        pre_tids, post_tids = self._treenodes_with_synapses()
        
        short_branches_with_inputs = {}
        short_branches_without_inputs = {}
        
        if synapse_type == 'input':
            tids_with_synapse = post_tids
        elif synapse_type == 'output':
            tids_with_synapse = pre_tids
        else:
            tids_with_synapse = {d: [] for d in all_datasets}
            for d, tids in post_tids.items():
                tids_with_synapse[d].extend(tids)
            for d, tids in pre_tids.items():
                tids_with_synapse[d].extend(tids)
        
        for d in all_datasets:
            short_branches_with_inputs[d] = defaultdict(list)
            short_branches_without_inputs[d] = defaultdict(list)
            for n, branch in short_branches[d]:
                number_of_inputs = len([tid for tid in branch[1:] if tid in tids_with_synapse[d]])
                if number_of_inputs > 0:
                    short_branches_with_inputs[d][n].append(branch)
                else:
                    short_branches_without_inputs[d][n].append(branch)
        
        return short_branches_with_inputs, short_branches_without_inputs

    def proportion_of_branches_with_synapses(self, f):
        short_branches_with_inputs, short_branches_without_inputs = self._short_branches_with_synapses()

        proportion_with_inputs = []

        for d in all_datasets:
            with_input = 0
            without_inputs = 0
            for n in list(short_branches_with_inputs[d].keys()) + list(short_branches_without_inputs[d].keys()):
                with_input += len(short_branches_with_inputs[d][n])
                without_inputs += len(short_branches_without_inputs[d][n])

            proportion_with_inputs.append(with_input/(with_input+without_inputs))

        x = [timepoint[d] for d in all_datasets]
        y = proportion_with_inputs

        self.plt.plot(
            'xy_graph', (x, (('', y), )), ylim=(0, 0.6),
            y_label='Proportion of spine-like\nprotrusions with synaptic input',
            save=f+'_proportion_of_branches_with_input', stats='spearmanr', linkpoints=False,
            size=0.22, x_label='Developmental age',
            margin={'left': 0.05, 'right': 0.01, 'top': 0.03, 'bottom': 0.03},
        )

    def _synapses_on_branches(self, remove_unknown=True, exclude_postemb=True):
        long_branches, short_branches = self._get_branches(exclude_postemb=exclude_postemb)
        synapses = data_manager.get_synapses_one_to_one().copy()
        
        if remove_unknown:
            valid_cells = set([n for c in neuron_list if in_brain(c) for n in class_members(c) if ntype(c) not in ('other', )])
            synapses = synapses[synapses['post'].isin(valid_cells)]
        
        if exclude_postemb:
            synapses = synapses[
                ~synapses['pre'].apply(is_postemb) &
                ~synapses['post'].apply(is_postemb)
            ]
        
        branch_tids = set()
        for d in all_datasets:
            for n, branch in short_branches[d]:
                branch_tids = branch_tids.union(branch[1:])
                
        synapses['post_is_branch'] = synapses['post_tid'].isin(branch_tids)

        return synapses

    def proportion_of_synapses_on_branches(self, f):

        synapses = self._synapses_on_branches()

        x = [timepoint[d] for d in all_datasets]
        y = synapses.groupby('dataset')['post_is_branch'].apply(lambda x: x.sum()/x.size)

        self.plt.plot(
            'xy_graph', (x, (('', y), )),
            y_label='Proportion of synapses\nonto spine-like protrusions', ylim=(0, 0.2),
            yticks=np.arange(0, 0.21, 0.05),
            save=f+'_proportion_of_synapses_on_branches',
            stats='spearmanr', linkpoints=False,
            size=0.22, x_label='Developmental age',
            margin={'left': 0.05, 'right': 0.01, 'top': 0.03, 'bottom': 0.03},
        )

    def branches_per_neuron(self, f):

        short_branches_with_inputs, short_branches_without_inputs = self._short_branches_with_synapses()
        valid_neurons = set([n for c in neuron_list if in_brain(c) for n in class_members(c) if ntype(n) not in ('other', )])
        
        number_of_branches = pd.DataFrame(index=valid_neurons, columns=all_datasets).fillna(0)
        number_of_branches_with_inputs = pd.DataFrame(index=valid_neurons, columns=all_datasets).fillna(0)
        proportion_of_branches_with_inputs = pd.DataFrame(index=valid_neurons, columns=all_datasets)
        for d in all_datasets:
            for n in list(short_branches_with_inputs[d].keys()) + list(short_branches_without_inputs[d].keys()):
                with_input = len(short_branches_with_inputs[d][n])
                without_inputs = len(short_branches_without_inputs[d][n])
                total = with_input + without_inputs
                if total == 0:
                    continue
                number_of_branches.loc[n, d] = total
                number_of_branches_with_inputs.loc[n, d] = with_input
                proportion_of_branches_with_inputs.loc[n, d] = with_input/total
                
        data_to_plot = (
            ('At birth', number_of_branches_with_inputs[['Dataset1', 'Dataset2']].mean(axis=1, skipna=False)),
            ('In adulthood', number_of_branches_with_inputs[['Dataset7', 'Dataset8']].mean(axis=1, skipna=False)),
        )
        
        self.plt.plot(
            'hist', data_to_plot, colors=[(0.5, 0.5, 0.5, 0.5), (0.8, 1, 0.0, .5)],
            size=0.22, x_label='Number of spine-like protrusions with inputs', 
            y_label='Number of neurons',
            save=f+'_branches_per_neuron',
            hist_range=(0, 16), bins=16,
            xlim=(0, 16), yticks=range(0, 151, 50),
            xticks=range(0, 20, 5),
            margin={'left': 0.05, 'right': 0.05, 'top': 0.05, 'bottom': 0.03},
        )

    def branches_by_type(self, f):
        synapses_branches = self._synapses_on_branches(exclude_postemb=False)
        
        synapses_branches['post_type'] = synapses_branches['post'].apply(ntype)
    
        df = synapses_branches.groupby(['dataset', 'post_type', 'post'])['post_is_branch'].apply(lambda x: x.sum()/x.size)
        df = df[['Dataset7', 'Dataset8']]
        df = df.reset_index().groupby(['post_type', 'post'], as_index=False).mean()
        
        cell_types = ['sensory', 'modulatory', 'inter', 'motor', 'muscle']
    
        self.plt.plot(
            'simple_violin', 
            [df.loc[df['post_type'] == cell_type, 'post_is_branch'] for cell_type in cell_types],
            colors=[dark_colors[t] for t in cell_types],
            labels=cell_types,
            ylim=(0, 1),
            size=0.22, 
            x_label='Postsynaptic cell type',
            y_label='Proportion of inputs on\nspine-like protrusions',
            save=f+'_branches_by_type',
        )

    def branch_locations(self, f):

        skeletons = data_manager.get_skeletons()

        long_branches, short_branches = self._get_branches()
        tids_long = defaultdict(dict)
        for d in long_branches:
            for n, branch in long_branches[d]:
                tids_long[d][n] = tids_long[d].get(n, set()).union(branch)
        
        short_branches_with_inputs, short_branches_without_inputs = self._short_branches_with_synapses(synapse_type='input')
        short_branches_with_outputs, short_branches_without_outputs = self._short_branches_with_synapses(synapse_type='output')
        short_branches_with_synapses, short_branches_without_synapses = self._short_branches_with_synapses(synapse_type='all')    
        
        location_of_branches_with_inputs = []
        location_of_branches_without_inputs = []
        for short_branches, location_of_branches in (
            (short_branches_with_inputs, location_of_branches_with_inputs),
            (short_branches_without_inputs, location_of_branches_without_inputs),
    #        (short_branches_with_outputs, location_of_branches_with_outputs),
    #        (short_branches_without_synapses, location_of_branches_without_synapses),
        ):
            for d in short_branches:
                if d not in ('Dataset7', 'Dataset8'):
                    continue
                for n, branches in short_branches[d].items():
                    for branch in branches:
                        if ntype(n) in ('other', 'muscle'):
                            continue
            
                        skel = skeletons[d][n]
                        tid = branch[0]
                        backbone_tids = tids_long[d][n]
            
                        while tid not in backbone_tids:
                            tid = skel['arbor'][tid]
                        
                        location_of_branches.append(skel['dist_to_root'][tid]/max(skel['dist_to_root'].values()))

        data_to_plot = (
            ('With inputs', location_of_branches_with_inputs),
            ('Without inputs', location_of_branches_without_inputs),
        )
        
        self.plt.plot(
            'kde', data_to_plot, colors=[(0.25, 0, 0.5), (0.5, 0.5, 0.5)],
            size=0.22, #x_label='Number of spine-like protrusions with inputs', 
            x_label='Location of spine-like protrusions',
            xticks=[0.16, 0.9], xtickpad=2,
            #no_y=True, 
            yticks=[], 
            no_x=True,
            y_label='Frequency', ypad=5,
            xticklabels=['Proximal', 'Distal'],
            save=f+'_branch_locations',
            xlim=(0, 1.1),
            #xlim=(0, 16), yticks=range(0, 151, 50),
            margin={'left': 0.05, 'right': 0.01, 'top': 0.05, 'bottom': 0.03},
        )
        
        print(sc.stats.mannwhitneyu(location_of_branches_with_inputs, location_of_branches_without_inputs))

    def branch_synapses_by_edge_classification(self, f, edge_classifications):

        edge_classifications = edge_classifications.copy()

        synapses_branches = self._synapses_on_branches(exclude_postemb=True)
        
        categories = {
            'increase': 'Dev.\ndynamic',
            'decrease': 'Dev.\ndynamic',
            'stable': 'Stable',
            'noise': 'Variable',
            'remainder': 'Variable'
        }
    
        synapses_with_classifications = synapses_branches.merge(edge_classifications, 'left', left_on=('pre', 'post'), right_index=True)
    
        assert synapses_with_classifications.isna().sum()['edge_classifications'] == 0
    
        synapses_with_classifications['edge_classifications'] = synapses_with_classifications['edge_classifications'].replace(categories)
    
        df = synapses_with_classifications[synapses_with_classifications['dataset'].isin(['Dataset7', 'Dataset8'])] \
            .groupby(['edge_classifications'])['post_is_branch'].agg(['sum', 'size'])
        df['proportion'] = df['sum']/df['size']
    
        classifications = ['Stable', 'Variable', 'Dev.\ndynamic']

        pairs = [
            ('Stable', 'Variable'), ('Stable', 'Dev.\ndynamic'), ('Variable', 'Dev.\ndynamic'),
        ]
        p_values = []
        for pair in pairs:
            _, p = proportions_ztest(
                [
                    df.loc[pair[0], 'sum'], 
                    df.loc[pair[1], 'sum']
                ], [
                    df.loc[pair[0], 'size'], 
                    df.loc[pair[1], 'size']
                ]
            )
            p_values.append(p)
            
        fdr = fdrcorrection0(p_values)
    
        adjusted_p_values = fdr[1]
    
        significant_pairs = list(itertools.compress(pairs, fdr[0]))
        significant_p_values = list(itertools.compress(adjusted_p_values, fdr[0]))
        
        # for c in classifications:
        #     print(c, df.loc[c])
        print('Significant pairs:')
        for pair, p in zip(significant_pairs, significant_p_values):
            print(' - '.join(pair).replace('\n', ' '), p)
            
        xy = [('', [
            [1, 2, 3],
            df.loc[classifications, 'proportion']
        ])]
        colors = ['#000000', '#eeeeee', '#007eff', '#007eff']
        stats = (None, None, 2.68, 0.3, '***')
        self.plt.plot(
            'bar_graph',
            xy,
            colors=colors,
            xticks=[1, 2, 3],
            xticklabels=classifications,
            y_label='Proportion of synapses onto\nspine-like protrusions',
            x_label='Connection classification',
            ylim=(0, 0.3),
            customstats=stats,
            legend=False,
            size=(0.11, 0.15),
            save=f+'_branch_synapses_by_edge_classification',
        )
