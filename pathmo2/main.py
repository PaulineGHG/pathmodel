import os
import csv
import clyngor


ROOT = os.path.dirname(__file__)


def reaction_creation(input_file, output_folder):
    '''
    Detect reaction sites by comparing molecules implied in a reaction.
    Return the result as a string.

    Args:
        input_file (str): path to the input data file
        output_folder (str): path to the output folder
    Returns:
        reaction_result (str): ASP answer as str
    '''
    print('~~~~~Creation of Reaction~~~~~')
    reaction_site_extraction_script = os.path.join(*[ROOT, 'asp', 'ReactionSiteExtraction.lp'])
    reaction_solver = clyngor.solve([input_file, reaction_site_extraction_script], use_clingo_module=False)
    reaction_results = []
    transformation_reactants = {}
    transformation_products = {}
    reactions = []
    for atom in next(reaction_solver.parse_args.int_not_parsed.sorted):
        reaction_results.append(atom[0] + '(' + ','.join(atom[1]) + ')')
        if 'diff' in atom[0]:
            reaction_id = atom[1][0]
            substructures = atom[1][1:]
            reactions.append(reaction_id)
            if 'Before' in atom[0]:
                if reaction_id not in transformation_reactants:
                    transformation_reactants[reaction_id] = [substructures]
                else:
                    transformation_reactants[reaction_id].append([substructures])
            elif 'After' in atom[0]:
                if reaction_id not in transformation_products:
                    transformation_products[reaction_id] = [substructures]
                else:
                    transformation_products[reaction_id].append([substructures])

    reactions = set(reactions)
    pathmodel_output_transformation_path = os.path.join(output_folder, 'pathmodel_data_transformations.tsv')
    with open(pathmodel_output_transformation_path, 'w') as transformation_file:
        csvwriter = csv.writer(transformation_file, delimiter = '\t')
        csvwriter.writerow(['reaction_id', 'reactant_sbustructure', 'product_substructure'])
        for reaction in reactions:
            if reaction in transformation_reactants:
                reactant = transformation_reactants[reaction]
            else:
                reactant = []
            if reaction in transformation_products:
                product = transformation_products[reaction]
            else:
                product = []
            csvwriter.writerow([reaction, reactant, product])


    reaction_result = '\n'.join([atom+'.' for atom in reaction_results])

    return reaction_result