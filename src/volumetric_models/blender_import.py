import bpy

import os
import math


path = '/path/to/smoothed/objs'


obj_paths = [f for f in sorted(os.listdir(path)) if f.lower().endswith('.obj')]


for f in obj_paths:
    bpy.ops.import_scene.obj(filepath=os.path.join(path, f))


for o in bpy.data.objects:
    if 'Dataset5' in obj_path:
        o.rotation_euler.y += math.pi/2
    if 'Dataset2' in obj_path:
        o.rotation_euler.x = 4.71238899230957
        o.rotation_euler.y = 0.62831842899322
    if 'Dataset8' in obj_path or 'Dataset6' in obj_path \
    or 'Dataset1' in obj_path or 'Dataset3' in obj_path \
    or 'Dataset5' in obj_path or 'Dataset4' in obj_path:
        o.rotation_euler.y += math.pi/2

for o in bpy.data.objects:
    o.name = o.name.split(' (')[0].split('.')[-1].split('_')[-1]


for o in bpy.data.objects:
    if o.name.startswith('BWM'):
        o.name ='B' + o.name[-2:] + o.name[-4:-2]

for o in bpy.data.objects:
    if o.name[:3] not in (
        'ADF', 'ADL', 'AFD', 'ALM', 'ALN', 'AQR', 'ASE', 'ASG', 'ASH', 'ASI',
        'ASJ', 'ASK', 'AUA', 'AVM', 'AWA', 'AWB', 'AWC', 'BAG', 'DVA', 'FLP',
        'IL2', 'OLL', 'OLQ', 'PLM', 'PLN', 'PVD', 'SAA', 'SDQ', 'URB', 'URX',
        'URY', 'unk', 'ADA', 'AIA', 'AIB', 'AIN', 'AIY', 'AIZ', 'AVA', 'AVB',
        'AVD', 'AVE', 'AVG', 'BDU', 'DVC', 'PVC', 'PVN', 'PVP', 'PVR', 'PVT',
        'RIA', 'RIB', 'RIF', 'RIG', 'RIH', 'RIM', 'RIP', 'RIR', 'IL1', 'RIV',
        'RMD', 'RME', 'RMF', 'RMH', 'SIA', 'SIB', 'SMB', 'SMD', 'URA', 'ADE',
        'AIM', 'ALA', 'AVF', 'AVH', 'AVJ', 'AVK', 'AVL', 'CEP', 'HSN', 'PDE',
        'PVQ', 'RIC', 'RID', 'RIS', 'RMG',
        'B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'fra'
    ):
        o.hide = True
        print('Hidden', o.name)
        continue

    if not o.name in bpy.data.materials:
        bpy.data.materials.new(o.name)
    o.data.materials.append(bpy.data.materials[o.name])
