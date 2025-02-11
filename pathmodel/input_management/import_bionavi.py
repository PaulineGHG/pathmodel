import os.path
import shutil

from pathmodel.plotting import *
from utils import *


PROD = 'product'
SUB = 'substrates'


def draw_bionavi(bionavi_output, output_path):
    with open(bionavi_output, 'r') as f:
        f.__next__()
        for l in f:
            l = l.strip().split('\t')
            rxn = l[0]
            score = round(float(l[3]), 2)
            sub_name = l[2]
            prod_name = l[1]
            output_sub = os.path.join(output_path, f'RXN{rxn}_{score}_{sub_name}_sub.png')
            output_prod = os.path.join(output_path, f'RXN{rxn}_{score}_{prod_name}_prod.png')
            smiles_to_2d_structure(smiles=l[8], output=output_sub)
            smiles_to_2d_structure(smiles=l[7], output=output_prod)


def extract_bionavi_rxn(bionavi_output):
    bio_rxn_info = {}
    bio_mol = {}
    bio_rxn = []
    with open(bionavi_output, 'r') as f:
        f.__next__()
        for l in f:
            l = l.strip().split('\t')
            rxn = l[0]
            if rxn not in bio_rxn_info:
                bio_rxn_info[rxn] = {SUB: {}, PROD: {}}

            sub_name = l[2]
            prod_name = l[1]

            bio_rxn_info[rxn][SUB][sub_name] = l[8]
            bio_rxn_info[rxn][PROD][prod_name] = l[7]

            if sub_name in bio_mol:
                print(bio_mol[sub_name] == mol_to_asp(sub_name, l[8], SMILES))
            if prod_name in bio_mol:
                print(bio_mol[prod_name] == mol_to_asp(prod_name, l[7], SMILES))

            bio_mol[sub_name] = mol_to_asp(sub_name, l[8], SMILES)
            bio_mol[prod_name] = mol_to_asp(prod_name, l[7], SMILES)

            bio_rxn.append(f'("bionavi_{rxn}", "{sub_name}", "{prod_name}").')

        return bio_rxn, bio_mol


def import_bionavi(bionavi_output, pathmodel_input):
    pm_file = pathmodel_input.replace('.lp', '_bionavi.lp')
    shutil.copyfile(pathmodel_input, pm_file)
    bio_rxn, bio_mol = extract_bionavi_rxn(bionavi_output)
    with open(pm_file, 'a') as f:
        f.write('\n%Bionavi mol\n%%%%%%%%%%%%\n')
        for mol_name, mol in bio_mol.items():
            f.write(f'\n%{mol_name}\n')
            for mol_info in mol:
                f.write(mol_info + '\n')
        f.write('\n%Bionavi rxn\n%%%%%%%%%%%%\n')
        for rxn in bio_rxn:
            print(rxn)


BN_OUTPUT_FILE = '/home/phamongi/Documents/Analysis_files/BioNavi/Oxylipins/pathway.txt'
PM_INPUT_FILE = 'Inputs/oxylipins_pwy.lp'
# draw_bionavi(BN_OUTPUT_FILE, 'Outputs/Bionavi/oxy')
import_bionavi(BN_OUTPUT_FILE, PM_INPUT_FILE)
