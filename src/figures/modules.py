import numpy as np

from src.data import data_manager
from src.data.dataset_info import all_datasets, stage

from src.plotting import plotter

community_colors = {
    'Posterior sensory': '#C891FF',
    'Medial sensory/\ninterneuron': '#6DB6FF',
    'Anterior sensory': '#9DDC8A',
    'Head movement': '#A68259',
    'Body movement': '#FFD265',
    'Motor output': '#FF8B8B',
    '': '#eeeeee',
}


class Figure(object):

    def __init__(self, output_path, page_size=7.20472):
        self.plt = plotter.Plotter(output_path=output_path, page_size=page_size)

    def modules_across_development(self, f):

        communities, community_names = data_manager.get_communities()

        community_counts = [(c, [len(communities[d][i]) for d in all_datasets]) for i, c in enumerate(community_names)]

        colors = community_colors

        time_labels = [stage[d] for d in all_datasets]

        lines = []
        dataset_counts = [len([n for c in communities[d] for n in c]) for d in all_datasets]
        cumulative_counts = np.zeros(len(all_datasets))
        for i, (c, counts) in enumerate(community_counts):

            for j, count in enumerate(counts):
                if j > 1 and count > 0:
                    prev_count = counts[j-1]
                    if prev_count == 0:
                        if c == 'Body movement':
                            prev_count += [counts2 for c, counts2 in community_counts if c == 'Head movement'][0][j-1]
                        elif c == 'Anterior sensory':
                            prev_count += [counts2 for c, counts2 in community_counts if c == 'Medial sensory/\ninterneuron'][0][j-1]
                        elif c == 'Posterior sensory':
                            prev_count -= [counts2 for c, counts2 in community_counts if c == 'Medial sensory/\ninterneuron'][0][j-1]
                    lines.append((
                        (j-1, (prev_count/2.0 + cumulative_counts[j-1])/dataset_counts[j-1]),
                        (j, (count/2.0 + cumulative_counts[j])/dataset_counts[j])
                    ))

            cumulative_counts += counts

        ds0_points = [community_counts[0][1][0]*0.5/dataset_counts[0],
                      (community_counts[0][1][0]+community_counts[1][1][0]*0.5)/dataset_counts[0]]
        ds1_points = [c1[1] for c1, c2 in lines if c1[0] == 1]

        lines.append((
            (0, ds0_points[1]),
            (1, ds1_points[1])
        ))
        lines.append((
            (0, ds0_points[1]),
            (1, ds1_points[2])
        ))
        lines.append((
            (0, ds0_points[0]),
            (1, ds1_points[2])
        ))
        lines.append((
            (0, ds0_points[0]),
            (1, ds1_points[0])
        ))

        for i, (name, counts) in enumerate(community_counts):
            community_counts[i] = (name, counts)

        self.plt.plot(
            'stacked_bar_graph', community_counts, size=0.25,
            margin={'left': 0.02, 'right': 0.08, 'top': 0.01, 'bottom': 0.03},
            nospacing=True, directlegend=True, colors=colors, xtickpad=2,
            y_label='Cells per module', x_label='Dataset', xlabels=time_labels,
            lines=lines,
            save=f+'_modules_across_development',
        )
