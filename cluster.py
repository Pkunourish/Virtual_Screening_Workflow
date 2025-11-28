import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, rdFingerprintGenerator
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
import matplotlib.pyplot as plt
from sklearn.manifold import MDS

# Process: 2D-clustering
# Configuration
input_file = 'admet_1443_select.csv'
output_file = './clsuter_1443_select.csv'

df_map = pd.read_csv('./1443_admet_select.csv')

def clusterfps(fps, cutoff=0.6):
    #Cluster function
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1-x for x in sims])
    
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    cs = sorted(cs, key=len, reverse=True)
    return cs

compound = []
smile = []
compound_df = pd.read_csv(input_file)
compound_df.insert(0, 'index', range(1, len(compound_df) + 1))

for _, number, smiles in compound_df[["index", "standardized_smiles"]].itertuples():
    compound.append((number, Chem.MolFromSmiles(smiles)))
    smile.append((number, smiles))

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
fps = [gen.GetFingerprint(x) for id, x in compound]

cluster = clusterfps(fps, cutoff=0.6)
print("Number of clusters:", len(cluster))
print("Number of molecules in largest cluster:", len(cluster[0]))

# Export cluster representative molecules(best PharmFit score)
cluster_first = []
for i in range(len(cluster)):
    cluster[i] = np.sort(cluster[i])
    list2 = [smile[x][1] for x in cluster[i]]
    cluster_first.append(list2[0])
cluster_first_df = pd.DataFrame(cluster_first, columns=['standardized_smiles'])

cluster_select = cluster_first_df.merge(df_map, on = 'standardized_smiles', how = 'left')
cluster_select.to_csv(output_file, index=False)