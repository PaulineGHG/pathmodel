from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem import rdMolDescriptors


SMILES = 'SMILES'
INCHI = 'InChI'
ATOM_CORRESP = {'C': 'carb',
                'O': 'oxyg',
                'S': 'sulf',
                'N': 'nitr'}


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


def generate_pathmodel_mol_input(mol_name, mol_code, encoding=INCHI):
    mol_name = mol_name.replace(' ', '_')
    if encoding == SMILES:
        mol = Chem.MolFromSmiles(mol_code)
        sm = Chem.MolToSmiles(mol, canonical=True)
        mol = Chem.MolFromSmiles(sm)
    elif encoding == INCHI:
        mol = inchi.MolFromInchi(mol_code)
    else:
        raise ValueError(f'{encoding} not a valid encoding name, must be in : {[SMILES, INCHI]}')
    atom_lst = []
    for a in mol.GetAtoms():
        atom_lst.append(a.GetSymbol())
    atom_pos = simplify_list_with_index(atom_lst)
    for a, pos_lst in atom_pos.items():
        for pos in pos_lst:
            if pos[0] == pos[1]:
                print(f'atom("{mol_name}",{pos[0]},{ATOM_CORRESP[a]}).')
            else:
                print(f'atom("{mol_name}",{pos[0]}..{pos[1]},{ATOM_CORRESP[a]}).')
    for b in mol.GetBonds():
        print(f'bond("{mol_name}",'
              f'{str(b.GetBondType()).lower()},'
              f'{b.GetBeginAtom().GetIdx() + 1},'
              f'{b.GetEndAtom().GetIdx() + 1}).')

    functional_groups = rdMolDescriptors.Properties()
    props = functional_groups.ComputeProperties(mol)
    print("Groupements fonctionnels identifiés :", props)


LA4_N = 'leukotriene A4'
LA4_IK = 'UFPQIRYSPUYQHK-WAQVJNLQSA-M'
LA4_I = 'InChI=1S/C20H30O3/c1-2-3-4-5-6-7-8-9-10-11-12-13-15-18-19(23-18)16-14-17-20(21)22/h6-7,9' \
        '-13,15,18-19H,2-5,8,14,16-17H2,1H3,(H,21,22)/p-1/b7-6-,10-9-,12-11+,15-13+/t18-,19-/m0/s1'

# generate_pathmodel_mol_input(LA4_N, LA4_I)

EPOME_N = '12(13)-EpOME'
EPOME_I = 'InChI=1S/C18H32O3/c1-2-3-10-13-16-17(21-16)14-11-8-6-4-5-7-9-12-15-18(19)20/h8,11,16-' \
          '17H,2-7,9-10,12-15H2,1H3,(H,19,20)/b11-8-'
EPOME_S = 'CCCCCC1OC1C/C=C/CCCCCCCC(=O)O'

generate_pathmodel_mol_input(EPOME_N, EPOME_S, SMILES)

LC4_N = 'leukotriene C4'

CY_N = 'Cycloartenol'
CY_I = 'InChI=1S/C30H50O/c1-20(2)9-8-10-21(3)22-13-15-28(7)24-12-11-23-26(4,5)25(31)14-16-29(23)' \
       '19-30(24,29)18-17-27(22,28)6/h9,21-25,31H,8,10-19H2,1-7H3/t21-,22-,23+,24+,25+,27-,28+,2' \
       '9-,30+/m1/s1'

