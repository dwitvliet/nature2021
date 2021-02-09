all_datasets = (
    'Dataset1', 'Dataset2', 'Dataset3', 'Dataset4', 'Dataset5', 'Dataset6', 'Dataset7', 'Dataset8',
)

datasets_with_adj = (
    'Dataset1', 'Dataset2', 'Dataset3', 'Dataset4', 'Dataset5', 'Dataset6', 'Dataset8',
)

timepoint = {
    'Dataset1': 0,
    'Dataset2': 5,
    'Dataset3': 8,
    'Dataset4': 16,
    'Dataset5': 23,
    'Dataset6': 27,
    'Dataset7': 50,
    'Dataset8': 50,
}

stage = {
    'Dataset1': 'L1',
    'Dataset2': 'L1',
    'Dataset3': 'L1',
    'Dataset4': 'L1',
    'Dataset5': 'L2',
    'Dataset6': 'L3',
    'Dataset7': 'Adult',
    'Dataset8': 'Adult',
    'white_l4': 'L4',
    'white_adult': 'Adult',
}

# Dataset 6 coordinates need to be scaled due to shrinkage.
upscale_l3 = 1.1
