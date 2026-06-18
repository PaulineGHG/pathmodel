from rdkit.Chem import inchi
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions

SMILES = 'SMILES'
INCHI = 'InChI'
ATOM_CORRESP = {'C': 'c',
                'O': 'o',
                'S': 's',
                'N': 'n'}


def get_canonical_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    sm = Chem.MolToSmiles(mol, canonical=True)
    mol = Chem.MolFromSmiles(sm)
    # inchi = Chem.MolToInchi(mol)
    # mol = Chem.MolFromInchi(inchi)
    return mol


def simplify_list_with_index(lst):
    result = {}
    start = 0
    for i, a in enumerate(lst):
        if a != lst[start]:
            value = lst[start]
            if value not in result:
                result[value] = []
            result[value].append((start + 1, i))
            start = i
    value = lst[start]
    if value not in result:
        result[value] = []
    result[value].append((start + 1, len(lst)))
    return result


def mol_to_asp(mol_name, mol_code, encoding=INCHI, domain=False):
    atom_str = 'atom'
    bond_str = 'bond'
    if domain:
        atom_str += 'Domain'
        bond_str += 'Domain'
    if encoding == SMILES:
        mol = get_canonical_mol(mol_code)
    elif encoding == INCHI:
        mol = inchi.MolFromInchi(mol_code)
    else:
        raise ValueError(f'{encoding} not a valid encoding name, must be in : {[SMILES, INCHI]}')
    atom_lst = []
    atom_deg = []
    atom_n_deg = []
    for a in mol.GetAtoms():
        atom_lst.append(a.GetSymbol())
        atom_deg.append(len(a.GetNeighbors()))
        atom_n_deg.append('-'.join(sorted([str(len(n.GetNeighbors())) for n in a.GetNeighbors()])))
    # atom_pos = simplify_list_with_index(atom_lst)
    asp_val = []
    for i in range(len(atom_lst)):
        asp_val.append(f'{atom_str}("{mol_name}",'
                       f'{i+1},'
                       f'{ATOM_CORRESP[atom_lst[i]]}).')
    # for a, pos_lst in atom_pos.items():
    #     for pos in pos_lst:
    #         if pos[0] == pos[1]:
    #             asp_val.append(f'{atom_str}("{mol_name}",{pos[0]},{ATOM_CORRESP[a]}).')
    #         else:
    #             asp_val.append(f'{atom_str}("{mol_name}",{pos[0]}..{pos[1]},{ATOM_CORRESP[a]}).')
    for b in mol.GetBonds():
        asp_val.append(f'{bond_str}("{mol_name}",'
                       f'{str(b.GetBondType()).lower()},'
                       f'{b.GetBeginAtom().GetIdx() + 1},'
                       f'{b.GetEndAtom().GetIdx() + 1}).')
    return asp_val


def rxn_to_asp(rxn, reactants, products, direction, lp_f):
    for r in reactants:
        lp_f.write(f'reactant({rxn},"{r}","{direction}").\n')
    for p in products:
        lp_f.write(f'product({rxn},"{p}","{direction}").\n')


def smiles_to_2d_structure(smiles, output):
    if type(smiles) == str:
        mol = get_canonical_mol(smiles)
        for atom in mol.GetAtoms():
            atom.SetProp('atomNote', str(atom.GetIdx() + 1))
        Draw.MolToFile(mol, output, size=(300, 300), imageType='svg')
    elif type(smiles) == dict:
        mol_lst = []
        names_lst = []
        for n_mol, smiles in smiles.items():
            mol = get_canonical_mol(smiles)
            mol_lst.append(mol)
            names_lst.append(n_mol)
        x = Draw.MolsToGridImage(mol_lst, molsPerRow=3, legends=names_lst, useSVG=True,
                                 subImgSize=(350, 350))
        with open(output + '.svg', 'w') as f:
            f.write(x)


def draw_mol(mol, output):
    for atom in mol.GetAtoms():
        atom.SetProp('atomNote', str(atom.GetIdx() + 1))
    Draw.MolToFile(mol, output + '.svg', size=(300, 300), imageType='svg')


def draw_rxn(sub_smiles, prod_smiles, output):
    sub = get_canonical_mol(sub_smiles)
    prod = get_canonical_mol(prod_smiles)
    rxn = rdChemReactions.ChemicalReaction()
    rxn.AddReactantTemplate(sub)
    rxn.AddProductTemplate(prod)

    x = Draw.ReactionToImage(rxn, useSVG=True)
    with open(output + '.svg', 'w') as f:
        f.write(x)


