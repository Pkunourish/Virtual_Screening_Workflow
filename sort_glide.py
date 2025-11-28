import pandas as pd

# Process: Sort glide score and extract top molecules
# Configuration
input_file = "sort_1443_glidescore.csv"
output_file = "sort_1443_glidescore_top10.csv"

# Glide results are extracted glide score elsewhere.
# Molecules separated by chirality in glide are firsted gouped by title.
df = pd.read_csv(input_file)
result = df.groupby('title', as_index=False).apply(lambda x: x.loc[x['r_i_glide_gscore'].idxmin()])
result = result.sort_values('r_i_glide_gscore')
top10_mol = result[result['r_i_glide_gscore']<-5.85].copy()
top10_mol[['SMILES', 'title', 'r_i_glide_gscore']].copy().to_csv(output_file, index=False)