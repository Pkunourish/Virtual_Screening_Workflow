import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import Draw
from rdkit.Chem import FilterCatalog
from rdkit.Chem import Crippen

# Process: Remove PAINS and fliter with admet properties
# Configuration
output_file = "./1443_admet_select.csv"

df_map = pd.read_csv('./rsa_1443_select.csv')

#Calculate PAINS and clogP from rdkit
def detect_pains(smiles):
    mol = Chem.MolFromSmiles(smiles)
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)
    if mol is None:
        return False
    else:
        return catalog.HasMatch(mol)

def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    clogp = Crippen.MolLogP(mol)
    return clogp

df_map['PAINS'] = df_map['standardized_smiles'].apply(detect_pains)
df_map['clogP'] = df_map['standardized_smiles'].apply(calculate_properties)

df_map = df_map[df_map['PAINS']==False]
df_map = df_map[df_map['clogP']<=5.0]

df_map.to_csv(output_file, index=False)