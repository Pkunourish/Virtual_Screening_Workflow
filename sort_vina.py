import pandas as pd

# Process: Sort vina score and extract top molecules
# Configuration
input_file = './sort_1443_vinascore.csv'
output_file = 'sort_1443_vinascore_top10.csv'

# Vina results are extracted and sorted by vina score elsewhere.
res = pd.read_csv(input_file, sep = ' ')
top10_mol = res[res['score']<-7.6].copy()
top10_mol.to_csv(output_file, index=False)