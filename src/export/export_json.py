import os
import json

import pandas as pd

from src.data import data_manager
from src.data.dataset_info import all_datasets


def export_nemanode_connections(path):
    synapses = data_manager.get_synapses_one_to_one().copy().assign(typ=0, syn=1)
    gapjunctions = data_manager.get_gapjunctions().copy().assign(typ=2, syn=1)

    cols_map = {
        'n1': 'pre', 'n2': 'post',
        'n1_tid': 'pre_tid', 'n2_tid': 'post_tid',
        'catmaid_id': 'ids', 'syn_id': 'ids',
    }
    cols = ['dataset', 'pre', 'post', 'pre_tid', 'post_tid', 'ids', 'typ', 'syn']

    df = pd.concat([
        synapses.rename(cols_map, axis=1)[cols],
        gapjunctions.rename(cols_map, axis=1)[cols]
    ])

    for dataset in all_datasets:
        df_dataset = df[df['dataset'] == dataset].drop('dataset', axis=1)
        df_dataset_grouped = df_dataset.groupby(['pre', 'post', 'typ']).agg(list)

        json_output = []
        for (pre, post, typ), series in df_dataset_grouped.iterrows():
            json_output.append({
                'pre': pre,
                'post': post,
                'typ': typ,
                **series.to_dict()
            })

        fpath = os.path.join(path, f'witvliet_2020_{dataset[-1]}.json')
        with open(fpath, 'w') as f:
            f.write(json.dumps(json_output, indent=2, sort_keys=True))

        print('Exported', fpath)
