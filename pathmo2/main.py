import os
import csv
import clyngor
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from pathmo2.inputs_generation import *
from pathmo2.utils import *
from pathmo2.rxn_mapper_mapping import generate_input_transformations

ROOT = os.path.dirname(__file__)


def generate_transformations(run_path):
    """
    Detect reaction sites by comparing molecules implied in a reaction.
    Return the result as a string.

    """
    print('~~~~~Creation of Reaction~~~~~')
    transformations_script = os.path.join(*[ROOT, 'asp', 'TransformationsExtraction.lp'])
    input_file = os.path.join(run_path, INPUTS_DIR, 'input.lp')

    reaction_solver = clyngor.solve([input_file, transformations_script], use_clingo_module=False)
    result_atoms = []
    transformations = {}
    reactant_transformation_centers = {}
    product_transformation_centers = {}
    for atom in next(reaction_solver.parse_args.int_not_parsed.sorted):
        atom_str = atom[0] + '(' + ', '.join(atom[1]) + ')'
        result_atoms.append(atom_str)
        # Extract transformations
        if atom[0].endswith('Difference'):
            if atom[1][0] not in transformations:
                transformations[atom[1][0]] = []
            transformations[atom[1][0]].append(atom)
        # Extract Reactant Transformation Centers
        if atom[0].startswith('reactantTransformationCenter'):
            if atom[1][0] not in reactant_transformation_centers:
                reactant_transformation_centers[atom[1][0]] = []
            reactant_transformation_centers[atom[1][0]].append(atom)
        # Extract Product Transformation Centers
        if atom[0].startswith('productTransformationCenter'):
            if atom[1][0] not in product_transformation_centers:
                product_transformation_centers[atom[1][0]] = []
            product_transformation_centers[atom[1][0]].append(atom)

    print(result_atoms)
    print(transformations)
    print(product_transformation_centers)
    print(reactant_transformation_centers)

    product_mol, reactant_mol = export_transformations_patterns(product_transformation_centers, reactant_transformation_centers)
    for r, m in reactant_mol.items():
        draw_mol(m, os.path.join(run_path, INPUTS_DIR, r[1:-1] + '_reactant_RC'))
    for r, m in product_mol.items():
        draw_mol(m, os.path.join(run_path, INPUTS_DIR, r[1:-1] + '_product_RC'))

    source_test = 'C=CC1OC1CCCC(=O)O'
    target_test = 'C=CC(O)C(CCCC(=O)O)SC(=O)O'
    hypothetical_mol = Chem.MolFromSmiles(target_test)
    substructure_search(hypothetical_mol, product_mol, reactant_mol)

    return result_atoms, transformations, product_transformation_centers, reactant_transformation_centers


    # output_folder = os.path.join(run_path, OUTPUT_DIR)
    # pathmodel_output_transformation_path = os.path.join(output_folder, 'pathmodel_data_transformations.tsv')
    # with open(pathmodel_output_transformation_path, 'w') as transformation_file:
    #     csvwriter = csv.writer(transformation_file, delimiter = '\t')
    #     csvwriter.writerow(['reaction_id', 'reactant_sbustructure', 'product_substructure'])
    #     for reaction in reactions:
    #         if reaction in transformation_reactants:
    #             reactant = transformation_reactants[reaction]
    #         else:
    #             reactant = []
    #         if reaction in transformation_products:
    #             product = transformation_products[reaction]
    #         else:
    #             product = []
    #         csvwriter.writerow([reaction, reactant, product])
    #
    #
    # reaction_result = '\n'.join([atom+'.' for atom in reaction_results])
    #
    # return reaction_result


RDKIT_BONDS = {"single": Chem.BondType.SINGLE,
               "double": Chem.BondType.DOUBLE,
               "triple": Chem.BondType.TRIPLE,
               "aromatic": Chem.BondType.AROMATIC}


def export_transformations_patterns(product_transformation_centers, reactant_transformation_centers):
    product_mol = {}
    reactant_mol = {}
    for trans_name, trans_lst in product_transformation_centers.items():
        trans_atoms = [x for x in trans_lst if x[0] == 'productTransformationCenterAtom']
        trans_bonds = [x for x in trans_lst if x[0] == 'productTransformationCenterBond']
        mol = asp_to_mol(trans_atoms, trans_bonds)
        product_mol[trans_name] = mol

    for trans_name, trans_lst in reactant_transformation_centers.items():
        trans_atoms = [x for x in trans_lst if x[0] == 'reactantTransformationCenterAtom']
        trans_bonds = [x for x in trans_lst if x[0] == 'reactantTransformationCenterBond']
        mol = asp_to_mol(trans_atoms, trans_bonds)
        reactant_mol[trans_name] = mol
    return product_mol, reactant_mol


def asp_to_mol(atoms, bonds):
    re_numbering = {}
    mol = Chem.RWMol()
    mol_atom_pos = 0
    for asp_atom in atoms:
        atom_symbol = asp_atom[1][1].upper()
        atom_position = asp_atom[1][2]
        re_numbering[atom_position] = mol_atom_pos
        mol.AddAtom(Chem.Atom(atom_symbol))
        mol_atom_pos += 1
    for asp_atom in bonds:
        bond_atom1 = re_numbering[asp_atom[1][1]]
        bond_atom2 = re_numbering[asp_atom[1][2]]
        bond_order = RDKIT_BONDS[asp_atom[1][3]]
        if mol.GetBondBetweenAtoms(bond_atom1, bond_atom2) is None:
            mol.AddBond(bond_atom1, bond_atom2, bond_order)
    mol = mol.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def substructure_search(hypothetical_mol, product_pattern, reactant_pattern):
    for trans, pattern in reactant_pattern.items():
        if hypothetical_mol.HasSubstructMatch(pattern):
            substructure_matches = hypothetical_mol.GetSubstructMatch(pattern)
            for p_idx, m_idx in enumerate(substructure_matches):
                print(f"Pattern atom {p_idx+1} → Mol atom {m_idx+1}")
    for trans, pattern in product_pattern.items():
        if hypothetical_mol.HasSubstructMatch(pattern):
            substructure_matches = hypothetical_mol.GetSubstructMatch(pattern)
            for p_idx, m_idx in enumerate(substructure_matches):
                print(f"Pattern atom {p_idx + 1} → Mol atom {m_idx + 1}")

        for atom in hypothetical_mol.GetAtoms():
            atom.SetProp('atomNote', str(atom.GetIdx() + 1))
        Draw.MolToFile(hypothetical_mol, 'hyp_mol.svg', size=(300, 300), imageType='svg')

        for atom in pattern.GetAtoms():
            atom.SetProp('atomNote', str(atom.GetIdx() + 1))
        Draw.MolToFile(pattern, 'pat_mol.svg', size=(300, 300), imageType='svg')

        hit_bonds = []
        for bond in pattern.GetBonds():
            aid1 = substructure_matches[bond.GetBeginAtomIdx()]
            aid2 = substructure_matches[bond.GetEndAtomIdx()]
            hit_bonds.append(hypothetical_mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())
        d = rdMolDraw2D.MolDraw2DSVG(500, 500)
        rdMolDraw2D.PrepareAndDrawMolecule(d, hypothetical_mol, highlightAtoms=substructure_matches,
                                           highlightBonds=hit_bonds)
        d.FinishDrawing()
        svg = d.GetDrawingText()
        with open("highlight.svg", "w") as f:
            f.write(svg)

        # from skfp.fingerprints import MACCSFingerprint, PubChemFingerprint
        # fp_maccs = MACCSFingerprint(n_jobs=-1)
        # fp_maccs_count = MACCSFingerprint(count=True, n_jobs=-1)
        #
        # fp_pubchem = PubChemFingerprint(n_jobs=-1)
        # fp_pubchem_count = PubChemFingerprint(count=True, n_jobs=-1)
        #
        # mols = [hypothetical_mol, pattern]
        # X_maccs = fp_maccs.transform(mols)
        # X_maccs_count = fp_maccs_count.transform(mols)
        #
        # X_pubchem = fp_pubchem.transform(mols)
        # X_pubchem_count = fp_pubchem_count.transform(mols)
        # print("Binary MACCS:")
        # print(f"Shape: {X_maccs.shape}")
        # print(f"Example values: {X_maccs[0, -10:]}")
        # print()
        # print("Count MACCS:")
        # print(f"Shape: {X_maccs_count.shape}")
        # print(f"Example values: {X_maccs_count[0, -10:]}")
        # print()
        # print("Binary PubChem:")
        # print(f"Shape: {X_pubchem.shape}")
        # print(f"Example values: {X_pubchem[0, :10]}")
        # print()
        # print("Count PubChem:")
        # print(f"Shape: {X_pubchem_count.shape}")
        # print(f"Example values: {X_pubchem_count[0, :10]}")
        # print()


# ==================================================================================================

RUN_PATH = '/home/phamongi/Documents/Dev/pathmodel/Files'
# RUN_PATH = 'C:\\Users\\Octav\\PycharmProjects\\pathmodel\\Files'
RUN_NAME = 'ToyExemple'

SOURCE = {'Source': 'C=CC1OC1CCCC(=O)O'}
TARGET = {'Target': 'C=CC(O)C(CCCC(=O)O)SC(=O)O'}
MC_REF_RXN = []
MAPPINGS_REF = {'Reaction1': (['C(=O)(O)CC1OC1/C=C/C', 'C(S)(=O)O'],
                              ['C(=O)(O)CC(SC(O)=O)C(O)/C=C/C'])}


generate_input_transformations(MAPPINGS_REF)

# generate_input(RUN_NAME, RUN_PATH, SOURCE, TARGET, MC_REF_RXN)
# generate_transformations(os.path.join(RUN_PATH, RUN_NAME))


# SOURCE = {'linoleate': 'CCCCC\C=C/C\C=C/CCCCCCCC([O-])=O'}
# TARGET = {'12_hydroxy_13_glutation_OME':
#           'CCCCCC(SCC(NC(=O)CCC(N)C(=O)O)C(=O)NCC(=O)O)C(O)C/C=C/CCCCCCCC(=O)O'}
# MC_REF_RXN = ['LEUKOTRIENE-C4-SYNTHASE-RXN', 'RXN-8495', 'LIPOXYGENASE-RXN']

