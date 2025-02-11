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
        f.write(f'%*\nMOLECULES\n{100*"="}\n*%\n')
        for mol, smiles in data[MOL].items():
            f.write(f'\n%{mol}\n')
            asp_atoms = mol_to_asp(mol_name=mol, mol_code=smiles, encoding=SMILES)
            for atom in asp_atoms:
                f.write(f'{atom}\n')
        # DOMAIN
        f.write(f'\n%*\nSHARED DOMAIN\n{100 * "="}\n*%\n')
        for domain, smiles in data[DOM].items():
            asp_atoms = mol_to_asp(mol_name=domain, mol_code=smiles, encoding=SMILES, domain=True)
            for atom in asp_atoms:
                f.write(f'{atom}\n')
        # REACTIONS
        f.write(f'\n%*\nREACTIONS\n{100*"="}\n*%\n')
        for rxn, react in data[RXN].items():
            f.write(f'\n%{rxn}\n')
            f.write(rxn_to_asp(rxn, react[0], react[1]) + '\n')
        # PATHWAYS
        # SPECIES
        # MZ
        if data[MZR] is not None:
            f.write(f'\n%*\nMZ RATIO\n{100 * "="}\n*%\n\n')
            for mz in data[MZR]:
                f.write(f'mzfiltering({mz}).\n')
        # ABSENT
        if data[ABS] is not None:
            f.write(f'\n%*\nABSENT MOLECULES\n{100 * "="}\n*%\n\n')
            for abs_mol in data[ABS]:
                f.write(f'absentmolecules("{abs_mol}").\n')
        # GOAL
        f.write(f'\n%*\nSOURCE - GOAL\n{100*"="}\n*%\n')
        f.write(f'\n%SOURCE\nsource("{data[SRC]}").\n')
        f.write(f'\n%GOAL\ngoal(pathway("{data[SRC]}","{data[TRG]}"))')


INPUT_YML = '../../Files/Inputs/oxylipins_input.yaml'

load_from_yaml(INPUT_YML)
