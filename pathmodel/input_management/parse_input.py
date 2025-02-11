import yaml
import mols2grid

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


def load_from_yaml(input):
    with open(input) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)
    print(data)
    smiles_to_2d_structure(data[MOL], input.replace('.yaml', ''))
    for rxn, v in data[RXN].items():
        draw_rxn(data[MOL][v[0]], data[MOL][v[1]], rxn)


INPUT_YML = 'oxylipins_input.yaml'

load_from_yaml(INPUT_YML)
