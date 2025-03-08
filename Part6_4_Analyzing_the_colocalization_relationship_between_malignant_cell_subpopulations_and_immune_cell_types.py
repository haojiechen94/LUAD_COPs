# co-location of malignant subtypes and immune subtypes
# author: Xinyue Zhang & Haojie Chen
# date: 2025-03-08
# email: zhangxinyue2021@sinh.ac.cn & chenhaojie2017@sinh.ac.cn



import math
import re
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams, gridspec
import matplotlib as mpl
import pandas as pd
import PyComplexHeatmap as pch
import scipy.stats

rcParams['font.family'] = 'Arial'
rcParams['pdf.fonttype'] = 42
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 20)



if __name__=='__main__':

	# P12 for example
	donor_id = 'P12' 

	# ++++++++++++++++++++++++++++++++++++
	# load data

    # composition of immune subtypes of each spot
	immune_prop = pd.read_csv(
        f'/data/ST/CARD_deconvolution/immune_cell_types/{donor_id}_immune_celltypes_composition.txt',
        sep='\t', header=0, index_col=0)

	# composition of malignant subtypes of each spot
    mal_prop_f = f'data/inte_scRNA_ST/{donor_id}_malignant_subtypes_composition_RCTD.txt'
    mal_prop = pd.read_csv(mal_prop_f, sep='\t', header=0, index_col=0)

	# in tumor region
    T_spots_df = pd.read_csv(f'/data/ST/hand_labeled_tumor_barcodes/{donor_id}_T.txt',
                             sep='\t', header=None, index_col=False, names=['T_spots'])
    T_spots = set(T_spots_df['T_spots']) & set(immune_prop.index)

    tumor_spots = list(T_spots & set(mal_prop.index))


    # ++++++++++++++++++++++++++++++++++++++++++++++++++
    # calculate the PCC between immune cell types and subtypes proportion in each spot ##########################################

    immune_subtypes = ['CD4+ Naive T cells', 'CD4+ Treg', 'GZMK+ NK',
    					'GZMK+ CD8+ effector T cells', 'NKT', 'CD4+ Tem',
                        'GIMAP7+ CD8+ effector T cells', 'HSP1A1+ CD8+ effector T cells',
                        'XCL+ CD8+ T cells', 'CXCL13+ CD4+ T cells', 'FCGR3A+ AREG+ NK',
                        'GZMK+ CD4+ CTL', 'CD4+ Tcm', 'LTB+ NK', 'KLRB1+ CD8+ T cells',
                        'GIMAP7+ NK', 'Exhausted CD8+ T cells', 'Plasma cells',
                        'B cells', 'Langerhans cells',
                        'Monocyte derived macrophage', 'Alveolar macrophage',
                        'CLEC10A+ dendritic cell(cDC)', 'CCR7+ mature dendritic cells',
                        'Neutrophil', 'CD14+ monocytes', 'Plasmacytoid dendritic cell(pDC)',
                        'CD16+ monocytes', 'CLEC9A+ dendritic cells',
                        'Lung basophil mast cell']


    malignant_subtypes = [celltype for celltype in mal_prop.columns if 'Malignant' in celltype]

    scc_matrix = pd.DataFrame(index=immune_subtypes, columns=malignant_subtypes)
    p_matrix = pd.DataFrame(index=immune_subtypes, columns=malignant_subtypes)

    for sub_type in malignant_subtypes:
        temp1 = []
        temp_p = []
        for immune_type in immune_subtypes:
            subtype_prop = tumor_epi_prop[sub_type].tolist()
            immune_prop = tumor_donor_immune_prop[immune_type].tolist()
            scc, p = scipy.stats.spearmanr(subtype_prop, immune_prop)
            scc_matrix.at[immune_type, sub_type] = scc
            p_matrix.at[immune_type, sub_type] = p



    # +++++++++++++++++++++++++++++++++++++++++++++++
    # select immune subtypes to plt
    disfunctional_immune = ['CCR7+ mature dendritic cells', 'GZMK+ NK', 'LTB+ NK', 'XCL+ CD8+ T cells','GZMK+ CD8+ effector T cells',
                            'KLRB1+ CD8+ T cells', 'Exhausted CD8+ T cells']
    functional_immune = ['GIMAP7+ NK', 'FCGR3A+ AREG+ NK', 'HSP1A1+ CD8+ effector T cells', 'NKT','GIMAP7+ CD8+ effector T cells']
    
    scc_matrix.loc[disfunctional_immune, 'fun'] = 'disfunctional'
    scc_matrix.loc[functional_immune, 'fun'] = 'functional'
    scc_matrix = scc_matrix[pd.isnull(scc_matrix['fun'])==False]

    # sort
    scc_matrix.sort_values(by=['fun','Malignant_C3'], ascending=[True, False], inplace=True)
    malignant_subtypes = ['Malignant_C1', 'Malignant_C2', 'Malignant_C3']


    fig = plt.subplots(figsize=(4, 8))
    gs = mpl.gridspec.GridSpec(8, 2, width_ratios=[10, 1], height_ratios=[2, 2, 2, 2, 2, 2, 2, 2], hspace=0.1)
    ax = plt.subplot(gs[1:8, 0])
    plt.imshow(scc_matrix, cmap='coolwarm', vmin=vmin, vmax=vmax, aspect='auto')

    plt.yticks(list(range(len(scc_matrix))), scc_matrix.index, size=15)
    plt.xticks(list(range(len(malignant_subtypes))), malignant_subtypes, size=15, rotation=90)

    ax = plt.subplot(gs[2, 1])
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb1 = mpl.colorbar.ColorbarBase(ax, cmap=mpl.cm.coolwarm, norm=norm)
    cb1.set_label('SCC', size=20)
    cb1.set_ticks([-0.5, 0.5])
    cb1.ax.tick_params(labelsize=15)
    ax = plt.subplot(gs[3, 1])
    plt.annotate(' ', xy=(1, 0), ha='center', va='center', size=18)
    plt.axis('off')

    plt.show()






    



    