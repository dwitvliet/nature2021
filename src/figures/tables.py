import os

from src.data import data_manager
from src.data.neuron_info import neuron_list, class_members, in_brain, is_postemb, ntype, dark_colors
from src.data.dataset_info import all_datasets

from src.plotting import plotter


class Figure(object):

    def __init__(self, output_path, page_size=7.20472):
        self.path = output_path
        self.plt = plotter.Plotter(output_path=output_path, page_size=page_size)

    def neuron_table(self, f):
        fpath = os.path.join(self.plt.output_path, f + '_neurons.csv')
        with open(fpath, 'w') as f:
            f.write(','.join(['class', 'members', 'type', 'integration']) + '\n')
            for c in sorted(neuron_list):
                if not in_brain(c):
                    continue
                if c in ('CAN', ):
                    continue

                n = str(len(class_members(c)))
                t = 'glia' if c in ('GLR', 'CEPsh') else ntype(c)
                p = 'post-embryonic' if is_postemb(c) else 'embryonic'
                f.write(','.join([c, n, t, p]) + '\n')

        print(f'Saved to `{fpath}`')

    def cells_in_modules_table(self, f):

        communities, community_names = data_manager.get_communities()

        datasets = [
            'L1 (dataset 1)', 'L1 (dataset 2)', 'L1 (dataset 3)', 'L1 (dataset 4)',
            'L2 (dataset 5)', 'L3 (dataset 6)', 'Adult (dataset 7)', 'Adult (dataset 8)'
        ]
        columns = []
        for d in all_datasets:
            rows = []
            for c in communities[d]:
                if not c:
                    continue
                rows.append(c)
            columns.append(rows)
#        f.write(','.join(datasets) + '\n')
#        for i, d in enumerate(all_datasets):

        colors = {n: dark_colors[ntype(n)] for c in neuron_list for n in class_members(c)}

        self.plt.plot(
            'table', columns, colors=colors, header=datasets,
            row_names=[c for c in community_names if c],
            margin={'left': 0.005, 'right': 0.095, 'top': 0.005, 'bottom': 0.005},
            xticks='remove', yticks='remove',
            size=(0.9, 0.75),
            width=900, height=780,
            no_x=True, no_y=True,
            save=f+'_cells_in_modules',
        )
