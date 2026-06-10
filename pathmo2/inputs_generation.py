# IMPORTS
# --------------------------------------------------------------------------------------------------
import json
import os.path
from utils import mol_to_asp, rxn_to_asp

from rdkit.Chem import AllChem
from rdkit.Chem import Draw


# CONSTANTS
# --------------------------------------------------------------------------------------------------
CUR_PATH = os.path.dirname(os.path.realpath(__file__))
METACYC_MAPPING_FILE = os.path.join(CUR_PATH, 'data', 'atom-mappings-smarts.json')
METACYC_REACTIONS_FILE = os.path.join(CUR_PATH, 'data', 'reactions.json')
METACYC_COMPOUNDS_FILE = os.path.join(CUR_PATH, 'data', 'compounds.json')

INPUTS_DIR = 'Inputs'
OUTPUTS_DIR = 'Outputs'


# FUNCTIONS
# --------------------------------------------------------------------------------------------------

# Main Function
def generate_input(run_name, output_path, source, target, metacyc_ref_reactions):

    # MANAGE DIRECTORIES
    run_path = os.path.join(output_path, run_name)
    mapping_dir = os.path.join(run_path, INPUTS_DIR, 'Mappings')
    directories = [run_path,
                   os.path.join(run_path, INPUTS_DIR),
                   os.path.join(run_path, OUTPUTS_DIR),
                   mapping_dir]
    for directory in directories:
        if not os.path.exists(directory):
            os.mkdir(directory)

    # MANAGE FILES
    reactions_file = os.path.join(run_path, INPUTS_DIR, 'Reactions_references.tsv')
    with open(reactions_file, 'w') as f_rxn:
        f_rxn.write('\t'.join(['DataBase', 'Rxn ID', 'Reactants', 'Products', 'Direction', 'Mapping SMARTS']))
    input_lp_file = os.path.join(run_path, INPUTS_DIR, 'input.lp')

    # EXTRACT DATA
    rxn_data = {}
    cpd_data = {}
    map_data = {}
    if metacyc_ref_reactions:
        rxn_data, cpd_data, map_data = import_metacyc(metacyc_ref_reactions, rxn_data, cpd_data, map_data)

    # WRITE INPUT FILES
    with open(input_lp_file, 'w') as flp:
        write_reactions(rxn_data, flp)
        write_chemicals(source, target, cpd_data, flp)
        write_mappings(map_data, mapping_dir, flp)

    with open(reactions_file, 'a') as f_rxn:
        for rxn_id, rxn_info in rxn_data.items():
            f_rxn.write('\n' + '\t'.join([str(rxn_info[3]), rxn_id, str(rxn_info[0]),
                                          str(rxn_info[1]), str(rxn_info[2]), map_data[rxn_id]]))


# Source data extraction functions
def import_metacyc(metacyc_ref_reactions, rxn_data, cpd_data, map_data):
    # EXTRACT RXN DATA
    with open(METACYC_REACTIONS_FILE, 'r') as f:
        r_f = json.load(f)
        for rxn in metacyc_ref_reactions:
            data = r_f[rxn]
            if data['direction'] == 'LEFT-TO-RIGHT':
                reactants = data['left']
                products = data['right']
                direction = 'uni'
            if data['direction'] == 'RIGHT-TO-LEFT':
                reactants = data['right']
                products = data['left']
                direction = 'uni'
            if data['direction'] == 'REVERSIBLE':
                reactants = data['left']
                products = data['right']
                direction = 'bi'
            rxn_data[rxn] = (reactants, products, direction, 'MetaCyc')
    # EXTRACT CPD DATA
    all_cpd = set()
    for rxn, data in rxn_data.items():
        for c in data[0]:
            all_cpd.add(c)
        for c in data[1]:
            all_cpd.add(c)
    with open(METACYC_COMPOUNDS_FILE, 'r') as f:
        c_f = json.load(f)
        for cpd in all_cpd:
            cpd_data[cpd] = c_f[cpd]
    # EXTRACT MAP DATA
    with open(METACYC_MAPPING_FILE, 'r') as f:
        mc_map = json.load(f)
        for r in metacyc_ref_reactions:
            map_data[r] = mc_map[r]

    return rxn_data, cpd_data, map_data


# LP writing functions
def write_chemicals(source, target, rxn_chemicals_lst, lp_f):
    lp_f.write(f'\n%*\nCHEMICALS\n{100 * "="}\n*%\n')
    for chem_id, smile in {**source, **target, **rxn_chemicals_lst}.items():
        lp_f.write(f'\n% {chem_id}\n')
        asp_atoms = mol_to_asp(mol_name=chem_id, mol_code=smile, encoding='SMILES')
        for atom in asp_atoms:
            lp_f.write(f'{atom}\n')

    lp_f.write(f'\n%*\nSOURCE - GOAL\n{100*"="}\n*%\n')
    source = list(source)[0]
    lp_f.write(f'\n% SOURCE\nsource("{source}").\n')
    lp_f.write(f'\n% GOAL\ngoal(pathway("{source}","{list(target)[0]}")).\n')


def write_reactions(reactions_data, lp_f):
    lp_f.write(f'%*\nREACTIONS\n{100 * "="}\n*%\n')
    for rxn, data in reactions_data.items():
        reactants = data[0]
        products = data[1]
        direction = data[2]
        lp_f.write(f'\n% {rxn}\n')
        rxn_to_asp(rxn, reactants, products, direction, lp_f)


def write_mappings(mappings_data, mapping_dir, lp_f):
    lp_f.write(f'\n%*\nMAPPINGS\n{100 * "="}\n*%\n')
    for rxn_id, mapping in mappings_data.items():
        rxn = AllChem.ReactionFromSmarts(mapping)
        reactants = rxn.GetReactants()
        products = rxn.GetProducts()

        lp_f.write(f'\n% {rxn_id}\n')
        for r in reactants:
            asso_if = dict()
            for atom in r.GetAtoms():
                idx = atom.GetIdx()
                symbol = atom.GetSymbol()
                amap = atom.GetAtomMapNum()
                lp_f.write(f'atomMappingReactant("{rxn_id}","{symbol}",{amap}).\n')
                asso_if[idx] = amap
            for bond in r.GetBonds():
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                bond_type = str(bond.GetBondType()).lower()
                lp_f.write(f'bondMappingReactant("{rxn_id}","{asso_if[a1]}","{asso_if[a2]}","{bond_type}").\n')

        for p in products:
            asso_if = dict()
            for atom in p.GetAtoms():
                idx = atom.GetIdx()
                symbol = atom.GetSymbol()
                amap = atom.GetAtomMapNum()
                lp_f.write(f'atomMappingProduct("{rxn_id}","{symbol}",{amap}).\n')
                asso_if[idx] = amap
            for bond in p.GetBonds():
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                bond_type = str(bond.GetBondType()).lower()
                lp_f.write(
                    f'bondMappingProduct("{rxn_id}","{asso_if[a1]}","{asso_if[a2]}","{bond_type}").\n')

        draw_rxn(rxn, os.path.join(mapping_dir, rxn_id + '.png'))


# Mol Draw functions
def draw_rxn(rxn, output):
    d2d = Draw.MolDraw2DCairo(1600, 600)
    d2d.DrawReaction(rxn, highlightByReactant=False)
    png = d2d.GetDrawingText()
    open(output, 'wb+').write(png)


# --------------------------------------------------------------------------------------------------

RUN_NAME = 'OXY'
RUN_PATH = '/home/phamongi/Documents/Dev/pathmodel/Files'
SOURCE = {'linoleate': 'CCCCC\C=C/C\C=C/CCCCCCCC([O-])=O'}
TARGET = {'12_hydroxy_13_glutation_OME': 'CCCCCC(SCC(NC(=O)CCC(N)C(=O)O)C(=O)NCC(=O)O)C(O)C/C=C/CCCCCCCC(=O)O'}
MC_REF_RXN = ['LEUKOTRIENE-C4-SYNTHASE-RXN', 'RXN-8495', 'LIPOXYGENASE-RXN']

# generate_input(RUN_NAME, RUN_PATH, SOURCE, TARGET, MC_REF_RXN)