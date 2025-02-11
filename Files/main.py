import pathmodel

# pathmodel infer -i data.lp -o output_folder -s 100

INPUT_FILE = '/home/phamongi/Documents/Packages/pathmodel/pathmodel/data/oxylipins_pwy.lp'
OUTPUT_DIR = 'Outputs/oxylipins'

pathmodel.pathmodel_analysis(INPUT_FILE, OUTPUT_DIR, step_limit=100)
# pathmodel_plot -i OUT