import pandas as pd

# Process: Using PharmFit scores to fliter molecules.
# PharmFit scores are calculted elsewhere with PocketV4 and PharmFit.
# Configuration
input_file = "sort_1443_PharmFit.csv"
output_file = "pharm_1443_select.csv"

# Load RMSD mapping file. Title in PharmFit results are ended with '.pdbqt'.
df_map = pd.read_csv('./rmsd_1443_select.csv')
df_map['title_vina'] = df_map['title_vina'] + '.pdbqt'

pharm = pd.read_csv(input_file)

sort_pharm = df_map.merge(
    pharm[['title_vina', 'PharmFit_score']].copy(), 
    on = 'title_vina',
    how='left',
    suffixes=('', '_dup')
)

# PharmFit results is added to the map file.
# A 0.50 cutoff is applied. You can adjust it based on your PharmFit results.
final_pharm = sort_pharm.sort_values(by = 'PharmFit_score', ascending=False).copy()
final_pharm['PharmFit_score'] = pd.to_numeric(final_pharm['PharmFit_score'], errors='coerce')
final_pharm[final_pharm['PharmFit_score']>=0.50].copy().to_csv(output_file, index=False)

# Optional: extract molecule names for binding pose extractions. dos format.
final_pharm[final_pharm['PharmFit_score']>=0.50]['title_vina'].copy().to_csv('./1443_pharm_vinalist.txt', index=False, header=False)
final_pharm[final_pharm['PharmFit_score']>=0.50]['title_glide'].copy().to_csv('./1443_pharm_glidelist.txt', index=False, header=False)