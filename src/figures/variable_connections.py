from collections import defaultdict

import numpy as np

from src.data import data_manager
from src.data.neuron_info import ntype, dark_colors

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
        self.plt = plotter.Plotter(output_path=output_path, page_size=page_size)

    def variable_edge_prevalence_quantified(self, f, edge_classifications,
                                            inputs=False, variable=None,
                                            nonvariable=None):

        G = data_manager.get_connections()['count'].copy()
        G = data_manager.remove_postemb(G)
        
        classifications = edge_classifications.copy()

        variable_defined_by_adult = bool(variable or nonvariable)

        if not variable_defined_by_adult:
            nonvariable = defaultdict(float)
            variable = defaultdict(float)

            for (pre, post), syns in G.iterrows():
                if np.max(syns[-2:]) <= 0:
                    continue

                s = np.count_nonzero(syns[-2:])
                if s == 0:
                    continue

                classification = classifications[(pre, post)]
                if classification in ('stable', 'increase', 'decrease'):
                    nonvariable[post if inputs else pre] += s

                if classification in ('remainder', 'noise'):
                    variable[post if inputs else pre] += s



        y_label='Proportion of variable connection'

        if inputs:
            y_label += ' inputs'
            x_label='Cell type'
            size = (0.2*1.25, 0.12)
            margin_top = 0.1
            ylim = (0, 1.01)
        else:
            y_label += ' outputs'
            x_label='Neuron type'
            size = (0.2, 0.12)
            ylim = (0, 1)
            margin_top = 0.05

        ns = list(set(list(nonvariable.keys()) + list(variable.keys())))


        data, l, c = [], [], []
        types = ['sensory', 'modulatory', 'inter', 'motor',]
        if inputs:
            types += ['muscle']
        for typ in types:
            data.append([variable[n]/(nonvariable[n]+variable[n]) for n in ns if ntype(n) == typ])
            c.append(dark_colors[typ])
            l.append(typ.capitalize())

        # for d in data:
        #     print(np.mean(d), len(d))

        stats = ((1,2), (1,3), (4,3), (4, 5), (3,5), (2,5), (5, 1), ) if inputs else ((3, 2), (2, 1), (3, 4), (2, 4), (1, 4))
        self.plt.plot(
            'box_plot', data,
            margin={'left': 0.04, 'right': 0.01, 'top': margin_top, 'bottom': 0.02},
            colors=c, xticklabels=l, y_label=y_label, x_label=x_label, ylim=ylim,
            xtickpad=2, xpad=3,
            show_outliers=False, stats=stats, size=size,
            save=f+'_variable_edge_prevalence_quantified'+('_in' if inputs else ''), darkmedian=True
        )
