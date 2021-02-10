import re


neuron_list = (
    'ADA', 'ADE', 'ADF', 'ADL', 'AFD', 'AIA', 'AIB', 'AIM', 'AIN', 'AIY', 'AIZ',
    'ALA', 'ALM', 'ALN', 'AQR', 'ASE', 'ASG', 'ASH', 'ASI', 'ASJ', 'ASK', 'ASn',
    'AUA', 'AVA', 'AVB', 'AVD', 'AVE', 'AVF', 'AVG', 'AVH', 'AVJ', 'AVK', 'AVL',
    'AVM', 'AWA', 'AWB', 'AWC', 'BAG', 'BDU', 'CEP', 'DAn', 'DBn', 'DDn', 'DVA',
    'DVB', 'DVC', 'FLP', 'HSN', 'IL1', 'IL2', 'LUA', 'OLL', 'OLQ', 'PDA',
    'PDB', 'PDE', 'PHA', 'PHB', 'PHC', 'PLM', 'PLN', 'PQR', 'PVC', 'PVD', 'PVM',
    'PVN', 'PVP', 'PVQ', 'PVR', 'PVT', 'PVW', 'RIA', 'RIB', 'RIC', 'RID', 'RIF',
    'RIG', 'RIH', 'RIM', 'RIP', 'RIR', 'RIS', 'RIV', 'RMD', 'RME', 'RMF', 'RMG',
    'RMH', 'SAA', 'SAB', 'SDQ', 'SIA', 'SIB', 'SMB', 'SMD', 'URA', 'URB', 'URX',
    'URY', 'VAn', 'VBn', 'VCn', 'VDn',
    'BWM01', 'BWM02', 'BWM03', 'BWM04', 'BWM05', 'BWM06', 'BWM07', 'BWM08',
    'CAN', 'CEPsh', 'GLR', 'excgl', 'hyp'
)


def class_members(cls):
    if cls in (
        'AVG', 'DVC', 'PVR', 'PVT', 'RIH', 'RIR', 'DVA', 'AQR', 'AVM', 'PQR',
        'PVM', 'DVB', 'PDA', 'PDB', 'ALA', 'AVL', 'RID', 'RIS',
        'I3', 'I4', 'I5', 'I5', 'M1', 'M4', 'M5', 'MI', 'SABD'
    ):
        return [cls]
    if cls in (
        'ADA', 'AIA', 'AIB', 'AIN', 'AIY', 'AIZ', 'BDU', 'LUA', 'PVN', 'PVP',
        'PVW', 'RIA', 'RIB', 'RIF', 'RIG', 'RIM', 'RIP', 'AVA', 'AVD', 'AVE',
        'AVB', 'PVC', 'ADL', 'AFD', 'ASE', 'ASG', 'ASH', 'ASI', 'ASJ', 'ASK',
        'AUA', 'AWA', 'AWB', 'AWC', 'BAG', 'FLP', 'OLL', 'URB', 'RMG', 'PDE',
        'ALM', 'ALN', 'PHA', 'PHB', 'PHC', 'PLM', 'PLN', 'PVD', 'SDQ', 'RIV',
        'RMF', 'RMH', 'AIM', 'AVF', 'AVH', 'AVJ', 'AVK', 'PVQ', 'RIC', 'ADE',
        'ADF', 'HSN', 'URX',
        'I1', 'I2', 'M2', 'M3', 'MC', 'NSM',
        'CAN',
        'SAAD', 'SAAV', 'URYD', 'URYV', 'SMBD', 'SMBV', 'SMDD', 'SMDV', 'URAD',
        'URAV', 'SIBD', 'SIBV', 'SIAD', 'SIAV', 'CEPD', 'CEPV', 'OLQD', 'OLQV',
        'IL1D', 'IL1V', 'IL2D', 'IL2V', 'RMDD', 'RMDV', 'GLRD', 'GLRV', 'CEPshD',
        'CEPshV', 'SABV'
    ):
        return [cls + n for n in ('L', 'R')]
    if cls in (
        'SAA', 'URY', 'SMB', 'SMD', 'URA', 'SIB', 'SIA', 'CEP', 'OLQ', 'CEPsh'
    ):
        return [cls + n for n in ('DL', 'DR', 'VL', 'VR')]
    if cls in ('IL1', 'IL2', 'RMD', 'GLR'):
        return [cls + n for n in ('DL', 'DR', 'L', 'R', 'VL', 'VR')]
    if cls in ('GLRL/R', 'RMDL/R', 'RMEL/R', 'IL1L/R', 'RMED/V', 'IL2L/R'):
        return [cls[:3] + cls[3], cls[:3] + cls[5]]

    if cls == 'SAB':
        return ['SABD', 'SABVL', 'SABVR']
    if cls == 'RME':
        return ['RMED', 'RMEL', 'RMER', 'RMEV']
    if cls == 'DAn':
        return ['DA' + str(i + 1) for i in range(9)]
    if cls == 'DBn':
        return ['DB' + str(i + 1) for i in range(7)]
    if cls in ('DDn', 'VCn'):
        return [cls[:2] + str(i + 1) for i in range(6)]
    if cls == 'VAn':
        return ['VA' + str(i + 1) for i in range(12)]
    if cls in ('VBn', 'ASn'):
        return [cls[:2] + str(i + 1) for i in range(11)]
    if cls == 'VDn':
        return ['VD' + str(i + 1) for i in range(13)]

    if cls == 'muscle':
        return ['muscle'] + \
           ['BWM-DL' + str(i + 1).zfill(2) for i in range(24)] + \
           ['BWM-DR' + str(i + 1).zfill(2) for i in range(24)] + \
           ['BWM-VL' + str(i + 1).zfill(2) for i in range(23)] + \
           ['BWM-VR' + str(i + 1).zfill(2) for i in range(24)]

    if cls in (
        'BWM01', 'BWM02', 'BWM03', 'BWM04', 'BWM05', 'BWM06', 'BWM07', 'BWM08',
    ):
        return ['BWM-' + n + cls[-2:] for n in ('DL', 'DR', 'VL', 'VR')]
    if cls in (
        'BWM01D', 'BWM02D', 'BWM03D', 'BWM04D', 'BWM05D', 'BWM06D', 'BWM07D', 'BWM08D',
        'BWM01V', 'BWM02V', 'BWM03V', 'BWM04V', 'BWM05V', 'BWM06V', 'BWM07V', 'BWM08V',
    ):
        return ['BWM-' + cls[-1] + n + cls[-3:-1] for n in ('L', 'R')]

    return [cls]


def ntype(n):
    n = nclass(n)
    if n not in neuron_list:
        n = nclass(class_members(n)[0])

    if n in (
        'ADF', 'ADL', 'AFD', 'ALM', 'ALN', 'AQR', 'ASE', 'ASG', 'ASH', 'ASI',
        'ASJ', 'ASK', 'AUA', 'AVM', 'AWA', 'AWB', 'AWC', 'BAG', 'DVA', 'FLP',
        'IL2', 'OLL', 'OLQ', 'PHA', 'PHB', 'PHC', 'PLM', 'PLN', 'PQR', 'PVD',
        'PVM', 'SAA', 'SDQ', 'URB', 'URX', 'URY'
    ):
        return 'sensory'
    if n in (
        'ADA', 'AIA', 'AIB', 'AIN', 'AIY', 'AIZ', 'AVA', 'AVB', 'AVD', 'AVE',
        'AVG', 'BDU', 'LUA', 'PVC', 'PVP', 'PVR', 'PVT', 'PVW',
        'RIA', 'RIB', 'RIF', 'RIG', 'RIH', 'RIM', 'RIR', 'RIP', 'AVJ',
    ):
        return 'inter'
    if n in (
        'ASn', 'DAn', 'DBn', 'DDn', 'DVB', 'IL1', 'PDA', 'PDB', 'RIV', 'RMD',
        'RME', 'RMF', 'RMH', 'SAB', 'SIA', 'SIB', 'SMB', 'SMD', 'URA', 'VAn',
        'VBn', 'VCn', 'VDn',
    ):
        return 'motor'
    if n in (
        'ADE', 'AIM', 'ALA', 'AVF', 'AVH', 'AVK', 'AVL', 'CEP', 'HSN',
        'PDE', 'PVQ', 'PVN', 'RIC', 'RID', 'RIS', 'RMG', 'DVC',
    ):
        return 'modulatory'
    if n in (
        'BWM01', 'BWM02', 'BWM03', 'BWM04', 'BWM05', 'BWM06', 'BWM07', 'BWM08'
    ):
        return 'muscle'
    if n in ('CAN', 'CEPsh', 'GLR', 'excgl', 'hyp'):
        return 'other'
    print(n, 'is not a valid neuron')
    return 'nonvalid'


def is_neuron(n):
    n = nclass(n)
    return n in neuron_list and ntype(n) not in ('muscle', 'other')


def nclass(n):
    if n in (
        'AVG', 'DVC', 'PVR', 'PVT', 'RIH', 'RIR', 'DVA', 'AQR', 'AVM', 'PQR',
        'PVM', 'DVB', 'PDA', 'PDB', 'ALA', 'AVL', 'RID', 'RIS',
        'I3', 'I4', 'I5', 'I5', 'M1', 'M4', 'M5', 'MI'
    ):
        return n
    if len(n) == 4 and n[-1] in 'LR' and n[:3] in (
        'ADA', 'ADE', 'ADF', 'ADL', 'AFD', 'AIA', 'AIB', 'AIM', 'AIN', 'AIY',
        'AIZ', 'ALM', 'ALN', 'ASE', 'ASG', 'ASH', 'ASI', 'ASJ', 'ASK', 'AUA',
        'AVA', 'AVB', 'AVD', 'AVE', 'AVF', 'AVH', 'AVJ', 'AVK', 'AWA', 'AWB',
        'AWC', 'BAG', 'BDU', 'CAN', 'FLP', 'GLR', 'HSN', 'IL1', 'IL2', 'LUA',
        'OLL', 'PDE', 'PHA', 'PHB', 'PHC', 'PLM', 'PLN', 'PVC', 'PVD', 'PVN',
        'PVP', 'PVQ', 'PVW', 'RIA', 'RIB', 'RIC', 'RIF', 'RIG', 'RIM', 'RIP',
        'RIV', 'RMD', 'RMF', 'RMG', 'RMH', 'SDQ', 'URB', 'URX'
    ):
        return n[:3]
    if len(n) == 5 and n[-2:] in ('DL', 'DR', 'VL', 'VR') and n[:3] in (
        'CEP', 'GLR', 'IL1', 'IL2', 'OLQ', 'RMD', 'SAA', 'SIA', 'SIB', 'SMB',
        'SMD', 'URA', 'URY'
    ):
        return n[:3]
    if len(n) == 8 and re.match('BWM-[DV][LR]0[0-8]', n):
        return 'BWM' + n[-2:]
    if n in (
        'RMED', 'RMEL', 'RMER', 'RMEV', 'SABD', 'SABVL', 'SABVR',
    ):
        return n[:3]
    if n in (
        'CEPshDL', 'CEPshDR', 'CEPshVL', 'CEPshVR'
    ):
        return n[:5]
    if n[:2] in ('AS', 'VB', 'VA', 'VD') and n[2:] in map(str, range(12)):
        return n[:2] + 'n'
    if n in ('VA12', 'VD12', 'VD13'):
        return n[:2] + 'n'
    if re.match('^(DA[1-9])|(DB[1-7])|(DD[1-6])|(VC[1-6])$', n):
        return n[:2] + 'n'
    return n


def in_brain(n):
    n = nclass(n)
    assert n in neuron_list, n + ' is not a valid neuron'
    return n in (
        'ADA', 'ADE', 'ADF', 'ADL', 'AFD', 'AIA', 'AIB', 'AIM', 'AIN', 'AIY',
        'AIZ', 'ALA', 'ALM', 'ALN', 'AQR', 'ASE', 'ASG', 'ASH', 'ASI', 'ASJ',
        'ASK', 'AUA', 'AVA', 'AVB', 'AVD', 'AVE', 'AVF', 'AVH', 'AVJ', 'AVK',
        'AVL', 'AVM', 'AWA', 'AWB', 'AWC', 'BAG', 'BDU', 'BWM01', 'BWM02',
        'BWM03', 'BWM04', 'BWM05', 'BWM06', 'BWM07', 'BWM08', 'CAN', 'CEP',
        'CEPsh', 'DVA', 'DVC', 'FLP', 'GLR', 'HSN', 'IL1', 'IL2', 'OLL', 'OLQ',
        'PLN', 'PVC', 'PVN', 'PVP', 'PVQ', 'PVR', 'PVT', 'RIA', 'RIB', 'RIC',
        'RID', 'RIF', 'RIG', 'RIH', 'RIM', 'RIP', 'RIR', 'RIS', 'RIV', 'RMD',
        'RME', 'RMF', 'RMG', 'RMH', 'SAA', 'SDQ', 'SIA', 'SIB', 'SMB', 'SMD',
        'URA', 'URB', 'URX', 'URY'
    )


def npair(n):
    if n in (
        'AVG', 'DVC', 'PVR', 'PVT', 'RIH', 'RIR', 'DVA', 'AQR', 'AVM',
        'PQR',
        'PVM', 'DVB', 'PDA', 'PDB', 'ALA', 'AVL', 'RID', 'RIS',
        'I3', 'I4', 'I5', 'I5', 'M1', 'M4', 'M5', 'MI',
        'SABD', 'excgl'
    ):
        return n
    cls = nclass(n)
    if cls in (
        'ADA', 'AIA', 'AIB', 'AIN', 'AIY', 'AIZ', 'BDU', 'LUA', 'PVN', 'PVP',
        'PVW', 'RIA', 'RIB', 'RIF', 'RIG', 'RIM', 'RIP', 'AVA', 'AVD', 'AVE',
        'AVB', 'PVC', 'ADL', 'AFD', 'ASE', 'ASG', 'ASH', 'ASI', 'ASJ', 'ASK',
        'AUA', 'AWA', 'AWB', 'AWC', 'BAG', 'FLP', 'OLL', 'URB', 'RMG', 'PDE',
        'ALM', 'ALN', 'PHA', 'PHB', 'PHC', 'PLM', 'PLN', 'PVD', 'SDQ', 'RIV',
        'RMF', 'RMH', 'AIM', 'AVF', 'AVH', 'AVJ', 'AVK', 'PVQ', 'RIC', 'ADE',
        'ADF', 'HSN', 'URX',
        'I1', 'I2', 'M2', 'M3', 'MC', 'NSM',
        'CAN'
    ):
        return cls
    if cls in (
        'ASn', 'DAn', 'DBn', 'DDn', 'VAn', 'VBn', 'VCn', 'VDn'
    ):
        return n
    if cls in (
        'SAA', 'URY', 'SMB', 'SMD', 'URA', 'SIB', 'SIA', 'CEP', 'OLQ',
        'CEPsh'
    ):
        return n[:-1]
    if n[:-1] in (
        'SABV', 'IL1D', 'IL1V', 'IL2D', 'IL2V', 'RMDD', 'RMDV', 'GLRD',
        'GLRV',
    ):
        return n[:-1]
    if n in (
        'IL1L', 'IL1R', 'IL2L', 'IL2R', 'RMDL', 'RMDR', 'GLRL', 'GLRR',
        'RMEL', 'RMER'
    ):
        return n[:3] + 'L/R'
    if n in ('RMED', 'RMEV'):
        return 'RMED/V'

    if len(n) == 8 and re.match('BWM-[DV][LR]0[0-8]', n):
        return 'BWM' + n[-2:] + n[4]

    print(n, 'is not a valid cell?')

    return n


def contralateral(n):
    if nclass(n) == n:
        return n
    if n == 'RMED':
        return 'RMEV'
    if n == 'RMEV':
        return 'RMED'
    c = {'L': 'R', 'R': 'L'}
    if n.endswith('R') or n.endswith('L'):
        return n[:-1] + c[n[-1]]
    if n.startswith('BWM-'):
        return 'BWM-' + n[4] + c[n[5]] + n[6:]  # BWM-[VD][LR][0-9]+
    return n


def is_postemb(n):
    n = nclass(n)
    if n not in neuron_list:
        n = nclass(class_members(n)[0])
    assert n in neuron_list, n + ' is not a valid neuron'
    return n in (
        'ALN', 'AQR', 'ASn', 'AVF', 'AVM', 'DVB', 'HSN', 'PDA', 'PDB', 'PDE',
        'PHC', 'PLN', 'PQR', 'PVD', 'PVM', 'PVN', 'PVW', 'RMF', 'RMH', 'SDQ',
        'VAn', 'VBn', 'VCn', 'VDn'
    )


dark_colors = {
    'sensory': '#fda0fd',
    'inter': '#ff442f',
    'motor': '#5cafff',
    'modulatory': '#ffc000',
    'muscle': '#7ec95a',
    'other': '#d9d9d9',
    'nonvalid': '#000000',
}