import csv
from rxnmapper import RXNMapper
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
rxn_mapper = RXNMapper()


def extract_chemicals(chem_input):
    chem_d = dict()
    with open(chem_input, 'r') as f:
        dat = csv.reader(f, delimiter='\t')
        dat.__next__()
        for l in dat:
            chem_d[l[0]] = l[1]
    return chem_d


def extract_rxn(rxn_input):
    r_dict = dict()
    with open(rxn_input, 'r') as f:
        dat = csv.reader(f, delimiter='\t')
        dat.__next__()
        for l in dat:
            r_dict[l[0]] = (l[1], l[2])
    return r_dict


def generate_input_transformations(r_dict):
    rxn_smiles = []
    for r_id, r in r_dict.items():
        reactant_smiles = '.'.join(r[0])
        product_smiles = '.'.join(r[1])
        rxn_smiles.append(f'{reactant_smiles}>>{product_smiles}')
    results = rxn_mapper.get_attention_guided_atom_maps(rxn_smiles)
    i = 0
    for res in results:
        i += 1
        mapped_rxn = res['mapped_rxn']
        rxn = AllChem.ReactionFromSmarts(mapped_rxn)
        d2d = Draw.MolDraw2DCairo(1600, 600)
        d2d.DrawReaction(rxn)
        png = d2d.GetDrawingText()
        open('rxn_mapped' + str(i) + '.png', 'wb+').write(png)
