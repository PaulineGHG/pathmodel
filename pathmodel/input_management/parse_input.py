import yaml

from pathmodel.input_management.utils import *

MOL = 'Molecules'
DOM = 'Domain'
RXN = 'Reactions'
PWY = 'Pathways'
SPC = 'Species'
MZR = 'MZ'
ABS = 'Absent'
SRC = 'Source'
TRG = 'Target'


def load_from_yaml(input_file, draw=False):
    lp_input = input_file.replace('.yaml', '.lp')
    with open(input_file) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)
    print(data)
    if draw:
        smiles_to_2d_structure(data[MOL], input_file.replace('.yaml', ''))
        for rxn, v in data[RXN].items():
            draw_rxn(data[MOL][v[0]], data[MOL][v[1]], rxn)
    with open(lp_input, 'w') as f:
        # MOLECULES
        f.write('%*\nMOLECULES\n=========\n*%\n')
        for mol, smiles in data[MOL].items():
            f.write(f'\n%{mol}\n')
            asp_atoms = mol_to_asp(mol_name=mol, mol_code=smiles, encoding=SMILES)
            for atom in asp_atoms:
                f.write(f'{atom}\n')
        # DOMAIN
        # REACTIONS
        f.write('%*\REACTIONS\n=========\n*%\n')
        for rxn, react in data[RXN].items():
            pass
    


INPUT_YML = '/home/phamongi/Documents/Dev/pathmodel/Files/Inputs/oxylipins_input.yaml'

load_from_yaml(INPUT_YML)
