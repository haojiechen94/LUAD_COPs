import pandas as pd
import numpy as np
from sklearn.utils import shuffle

cell_type_proportions_df=pd.read_csv('../cell_types_composition.txt',sep='\t')
COPs_df=pd.read_csv('../meta_data.txt',sep='\t')


data=[]
celltypes=cell_type_proportions_df.index
C_IDs=list(set(COPs_df['seurat_clusters']))
for celltype in celltypes:
    temp=[]
    for C_ID in C_IDs:
        COP=COPs_df[COPs_df['seurat_clusters']==C_ID].index
        observed_mean=cell_type_proportions_df.loc[COP,celltype].mean()
        differences=[]
        for i in range(1000):
            shuffled_cell_type_proportions_df=shuffle(cell_type_proportions_df)
            shuffled_cell_type_proportions_df.index=cell_type_proportions_df.index
            permutated_mean=shuffled_cell_type_proportions_df.loc[COP,celltype].mean()
            difference=ture_mean-permutated_mean
            differences.append(difference)
        temp.append(np.mean(differences)/np.std(differences))
    data.append(temp)

relative_enrichment_scores=pd.DataFrame(data,index=celltypes,columns=C_IDs)
relative_enrichment_scores.to_csv('../relative_enrichment_or_depletion_scores.txt',sep='\t')