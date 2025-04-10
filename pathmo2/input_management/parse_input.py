import os.path
import csv
from utils import *


# CONSTANTS ----------------------------------------------------------------------------------------
INPUT_DIR = 'Inputs'
OUTPUT_DIR = 'Outputs'
CHEM_INPUT = 'Chemicals_input.tsv'
RXN_INPUT = 'Reactions_input.tsv'
LP_INPUT = 'input.lp'

CHEM = 'Chemical'
RXN = 'Reaction'
ROLE = 'ROLE'
MZ = 'MZ'
ORIGIN = 'ORIGIN'
LINKS = 'LINKS'
REACTANT = 'REACTANT'
PRODUCT = 'PRODUCT'
EC = 'EC'
TAXON = 'TAXON'

DOMAIN = 'DOMAIN'
ABSENT = 'ABSENT'
SOURCE = 'SOURCE'
TARGET = 'TARGET'

CHEM_COLUMNS = [CHEM, SMILES, ROLE, MZ, ORIGIN, LINKS]
RXN_COLUMNS = [RXN, REACTANT, PRODUCT, ORIGIN, LINKS, EC, TAXON]


#     if draw:
#         smiles_to_2d_structure(data[MOL], input_file.replace('.yaml', ''))
#         for rxn, v in data[RXN].items():
#             draw_rxn(data[MOL][v[0]], data[MOL][v[1]], rxn)

#         # MZ
#         if data[MZR] is not None:
#             f.write(f'\n%*\nMZ RATIO\n{100 * "="}\n*%\n\n')
#             for mz in data[MZR]:
#                 f.write(f'mzfiltering({mz}).\n')

#         # ABSENT
#         if data[ABS] is not None:
#             f.write(f'\n%*\nABSENT MOLECULES\n{100 * "="}\n*%\n\n')
#             for abs_mol in data[ABS]:
#                 f.write(f'absentmolecules("{abs_mol}").\n')


def init_run(name, output_path):
    run_dir = os.path.join(output_path, name)
    if not os.path.exists(run_dir):
        if os.path.exists(output_path):
            os.mkdir(run_dir)
            os.mkdir(os.path.join(run_dir, INPUT_DIR))
            os.mkdir(os.path.join(run_dir, OUTPUT_DIR))
            with open(os.path.join(run_dir, INPUT_DIR, CHEM_INPUT), 'w') as f:
                csvwriter = csv.writer(f, delimiter='\t')
                csvwriter.writerow(CHEM_COLUMNS)
                csvwriter.writerow(['', '', SOURCE, '', '', ''])
                csvwriter.writerow(['', '', TARGET, '', '', ''])
            with open(os.path.join(run_dir, INPUT_DIR, RXN_INPUT), 'w') as f:
                csvwriter = csv.writer(f, delimiter='\t')
                csvwriter.writerow(RXN_COLUMNS)
        else:
            raise FileNotFoundError(f'Directory {output_path} not found')
    else:
        raise FileExistsError(f'Directory {run_dir} already exists, not overwriting')
    return run_dir


def generate_lp_input(run_dir, draw=False):
    chem_input = os.path.join(run_dir, INPUT_DIR, CHEM_INPUT)
    rxn_input = os.path.join(run_dir, INPUT_DIR, RXN_INPUT)
    lp_input = os.path.join(run_dir, INPUT_DIR, LP_INPUT)

    with open(chem_input, 'r') as cf, open(rxn_input, 'r') as rf, open(lp_input, 'w') as lf:
        chemicals_lst = list(csv.DictReader(cf, delimiter='\t'))
        reactions_lst = list(csv.DictReader(rf, delimiter='\t'))
        write_chemicals(chemicals_lst, lf)
        write_reactions(reactions_lst, [ch[CHEM] for ch in chemicals_lst], lf)


def write_chemicals(chemicals_lst, lp_f):
    sources = [chem for chem in chemicals_lst if chem[ROLE] == SOURCE]
    targets = [chem for chem in chemicals_lst if chem[ROLE] == TARGET]
    domains = [chem for chem in chemicals_lst if chem[ROLE] == DOMAIN]
    env_chem = [chem for chem in chemicals_lst if chem[ROLE] == '' or chem[ROLE] is None]
    if len(sources) != 1:
        print(f'Exactly 1 source must be defined, {len(sources)} found.')
    if len(targets) < 1:
        print(f'At least 1 target must be specified, {len(targets)} found.')
    lp_f.write(f'%*\nCHEMICALS\n{100 * "="}\n*%\n')
    for chem in sources + targets + env_chem:
        lp_f.write(f'\n% {chem[CHEM]}\n')
        asp_atoms = mol_to_asp(mol_name=chem[CHEM], mol_code=chem[SMILES], encoding=SMILES)
        for atom in asp_atoms:
            lp_f.write(f'{atom}\n')
    lp_f.write(f'\n%*\nSHARED DOMAIN\n{100 * "="}\n*%\n')
    for dom in domains:
        lp_f.write(f'\n% {dom[CHEM]}\n')
        asp_atoms = mol_to_asp(mol_name=dom[CHEM], mol_code=dom[SMILES], encoding=SMILES,
                               domain=True)
        for atom in asp_atoms:
            lp_f.write(f'{atom}\n')

    lp_f.write(f'\n%*\nSOURCE - GOAL\n{100*"="}\n*%\n')
    source = sources[0][CHEM]
    lp_f.write(f'\n% SOURCE\nsource("{source}").\n')
    for trg in targets:
        lp_f.write(f'\n% GOAL\ngoal(pathway("{source}","{trg[CHEM]}"))')


def write_reactions(reactions_lst, chemicals_ids, lp_f):
    lp_f.write(f'\n\n%*\nREACTIONS\n{100 * "="}\n*%\n')
    for rxn in reactions_lst:
        if rxn[REACTANT] not in chemicals_ids:
            print(f'Reactant {rxn[REACTANT]} of reaction {rxn[RXN]} not defined in chemicals.')
        elif rxn[PRODUCT] not in chemicals_ids:
            print(f'Product {rxn[PRODUCT]} of reaction {rxn[RXN]} not defined in chemicals.')
        lp_f.write(f'\n% {rxn[RXN]}\n')
        lp_f.write(rxn_to_asp(rxn[RXN], rxn[REACTANT], rxn[PRODUCT]) + '\n')


# ==================================================================================================

RUN_PATH = '/home/phamongi/Documents/Dev/pathmodel/Files'
RUN_NAME = 'Oxylipins'

generate_lp_input(os.path.join(RUN_PATH, RUN_NAME))
