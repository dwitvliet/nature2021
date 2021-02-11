import os
from itertools import compress
from collections import defaultdict

import numpy as np
import pandas as pd
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
from statsmodels.stats.proportion import proportions_ztest

from src.data import data_manager
from src.data.neuron_info import ntype

from src.plotting import plotter


categories = {
    'increase': 'Developmental change',
    'decrease': 'Developmental change',
    'stable': 'Stable',
    'noise': 'Variable',
    'remainder': 'Variable'
}


class Figure(object):

    def __init__(self, output_path, page_size=7.20472):
        self.path = output_path
        self.plt = plotter.Plotter(output_path=output_path, page_size=page_size)


    def changes_type_overrepresentation(self, fname, edge_classifications, color='red'):

        G = data_manager.get_connections()['count'].copy()
        G = data_manager.remove_postemb(G)

        classifcations = edge_classifications.copy()

        edges_sample = [e for e, t in classifcations.items() if t in ('increase', 'decrease')]
        edges_population = [e for e, t in classifcations.items() if t in ('increase', 'decrease', 'stable')]

        G_type_sample = defaultdict(int)
        G_type_population = defaultdict(int)
        for (pre, post) in G.index:
            edge = (pre, post)
            edge_t = (ntype(pre), ntype(post))
            if edge in edges_sample:
                G_type_sample[edge_t] += 1
            if edge in edges_population:
                G_type_population[edge_t] += 1

        p_values = []
        labels = []
        signs = []
        totals = []

        type_count_change = defaultdict(int)
        type_count_stable = defaultdict(int)

        for edge_t in G_type_population:
#            n = G_type_population[edge_t]
#            if n < 3:
#                continue
            if edge_t[1] == 'other':
                continue

            for typ in [edge_t[0]]:
                type_count_change[typ] += G_type_sample[edge_t]
                type_count_stable[typ] += G_type_population[edge_t]-G_type_sample[edge_t]

            count = G_type_sample[edge_t]
            average = G_type_population[edge_t]*float(sum(G_type_sample.values()))/sum(G_type_population.values())
            total = G_type_population[edge_t]

            _, p = proportions_ztest([count, average], [total, total])

            totals.append(total)
            signs.append('under' if count < average else 'over')
            p_values.append(p)
            labels.append(edge_t)


        fdr = fdrcorrection0(p_values)

        adjusted_p_values = fdr[1]

        significant_totals = list(compress(totals, fdr[0]))
        significant_labels = list(compress(labels, fdr[0]))
        significant_p_values = list(compress(adjusted_p_values, fdr[0]))
        significant_signs = list(compress(signs, fdr[0]))

        print('Statistically significant over- or underrepresentation:')
        for l, p, s, t in zip(significant_labels, significant_p_values, significant_signs, significant_totals):
            print('_'.join(l), s, 'p =', p, 'n =', t)

        data = []
        types = ['sensory', 'modulatory', 'inter', 'motor', 'muscle']

        output = []
        for i, post in enumerate(types):
            x, y = [], []
            for j, pre in enumerate(types):
                edge = (pre, post)

                sample = G_type_sample[edge]
                population = G_type_population[edge]

                if population < 10:
                    continue

                proportion = 0
                if population:
                    proportion = sample/(population+0.0)
#                    proportion = sample/(population+0.0)*3-1.3 for variable
                red = int(min(255, 255*proportion*2))
                
                if color == 'red':
                    c = '#%02x%02x%02x' % (red, 0, 0)
                if color == 'orange':
                    c = '#%02x%02x%02x' % (red, red // 2, 0)
                if color == 'pink':
                    c = '#%02x%02x%02x' % (red, 0, red // 1.5)
#                color = '#%02x%02x%02x' % (255-red, 255-red, 255-red)

                output.append((pre,  post, population, c, proportion))

                if pre == 'muscle':
                    continue
                if population == 0:
                    continue

                x.append(j*(len(types)+1)+i)
                y.append(sample/float(population))

            data.append((post, (x, y)))

        edge_to_p = dict(zip(labels, adjusted_p_values))
        output = sorted(output, key=lambda x: (edge_to_p[(x[0], x[1])] < 0.05, x[2]))

        fpath = os.path.join(self.path, fname + '_changes_type_overrepresentation.txt')

        with open(fpath, 'w') as f:
            f.write('"rowname"\t"key"\t"value"\t"color"\t"percentage_changing"\n')
            for (pre,  post, population, c, proportion) in output:
                if pre in ('sensory', 'modulatory', 'inter', 'motor'):
                    if pre != 'inter':
                        pre += ' '
                    pre += 'neuron'
                if post in ('sensory', 'modulatory', 'inter', 'motor'):
                    if post != 'inter':
                        post += ' '
                    post += 'neuron'
                pre = pre.capitalize()
                post = post.capitalize()

                f.write('"{}"\t"{}"\t{}\t"{}"\t"{}"\n'.format(pre,  post, population, c, proportion*100))

        print(f'Saved to `{fpath}`')

    def relative_synapse_increase_by_type(self, fname, edge_classifications):

        G = data_manager.get_connections()['count'].copy()
        G = data_manager.remove_postemb(G)

        edge_classifications = edge_classifications.copy()
        
        G['pre_type'] = [ntype(n) for n in G.index.get_level_values('pre')]
        G['post_type'] = [ntype(n) for n in G.index.get_level_values('post')]
        G['type'] = G['pre_type'] + G['post_type']
        G = G.join(edge_classifications.rename('classification'))
        
        G['rel_synapse_addition'] = G[['Dataset7', 'Dataset8']].mean(axis=1) / G[['Dataset1', 'Dataset2']].mean(axis=1)
        
        G_stable = G[G['classification'] == 'stable']
        G_stable = G_stable[G_stable['rel_synapse_addition'] != np.inf]
        stable_synapse_addition = G_stable.groupby('type')['rel_synapse_addition'].mean()

        stable_synapse_addition_percentage = (stable_synapse_addition - 1) * 100
        stable_synapse_addition_percentage.name = 'percentage_synapse_number_increase_for_stable_connections'

        fpath = os.path.join(self.path, fname + '_synapse_increase_by_type.txt')
        stable_synapse_addition_percentage.to_csv(fpath)
        print(f'Saved to `{fpath}`')
