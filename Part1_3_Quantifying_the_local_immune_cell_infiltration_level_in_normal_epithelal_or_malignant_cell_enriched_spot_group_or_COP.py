# quantifying the local immune cell infiltration level in normal epithelial/malignant cell enriched spot group/COP
# author: Xinyue Zhang & Haojie Chen
# date: 2025-03-08
# email: zhangxinyue2021@sinh.ac.cn & chenhaojie2017@sinh.ac.cn


import re

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
import scipy.stats
import glob



rcParams['font.family'] = 'Arial'
rcParams['pdf.fonttype'] = 42
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows',20)


if __name__ == '__main__':

    work_dir = '/data/ST/CARD_deconvolution/'

    # ++++++++++++++++++++++++++++++++++++++++++++++++
    # composition of each celltype for each spot
    # immune cell types
    immune_composition_dict = {}
    for path in glob.glob(work_dir + '/immune_cell_types/*_immune_celltypes_composition.txt'):
        ID = '_'.join(path.split('\\')[-1].split('_')[:2])
        if ID != 'P12':
            temp_df = pd.read_csv(path, sep='\t')
            temp_df.index = [ID + '_' + i for i in temp_df.index]
            immune_composition_dict[ID] = temp_df
    immune_composition_df = pd.concat([immune_composition_dict[i] for i in immune_composition_dict])

    # major cell types
    major_composition_dict = {}
    for path in glob.glob(work_dir + '/major_cell_types/*_major_celltypes_composition.txt'):
        ID = '_'.join(path.split('\\')[-1].split('_')[:2])
        if ID != 'P12':
            temp_df = pd.read_csv(path, sep='\t')
            temp_df.index = [ID + '_' + i for i in temp_df.index]
            major_composition_dict[ID] = temp_df
    major_composition_df = pd.concat([major_composition_dict[i] for i in major_composition_dict])


    # +++++++++++++++++++++++++++++++++++++++++++++++++
    # example
    # COP label
    cop_meta_f = work_dir + 'CARD_deconvolution/major_cell_types/COPs_meta_data.txt'
    clustering_result = pd.read_csv(cop_meta_f, sep='\t', header=0, index_col=0)
    print(f"clusters or COPs: {set(clustering_result['COPs'])}")
    print('\ncovert cluster_id to COP_id:')
    stC_dic = {}
    for i in clustering_result.index:
        # C_ID = 'COP%d' % (clustering_result.loc[i, 'seurat_clusters'])
        C_ID = clustering_result.loc[i, 'COPs']
        if C_ID in stC_dic:
            stC_dic[C_ID].add(i)
        else:
            stC_dic[C_ID] = set([i])


	# ++++++++++++++++++++++++++++++++++++++++++++++++++++
	# composition of MIA and IAC spots for each COP
    print('\ncount num of MIA and IAC spots per COP')
    count_dic = {}
    for COP in stC_dic:
        count_dic[COP] = {'MIA': 0, 'IAC': 0}
        for i in stC_dic[COP]:
            ID = '_'.join(i.split('_')[:2])
            if ID in ['P01', 'P02', 'P03', 'P04']:
                count_dic[COP]['MIA'] += 1
            elif ID in ['P11', 'P12', 'P14']:
                count_dic[COP]['IAC'] += 1



	# ++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # composition of cell types in each COP

    # true mean_composition
    cell_types = ['Endothelia.cells', 'AT1', 'AT2', 'Fibroblast',
                  'TorNK.cells', 'B.cells',
                  'Myeloid.cells', 'Mast.cells',
                  'Ciliated.epithelia.cells', 'Club_cell', 'Malignent_cells']

    ordered_clusters = [0, 1, 2, 3, 4, 5, 6, 7, 8]  # COPs
    print('for the composition matrix:')
    print(f'cols: {ordered_clusters}')
    print(f'rows: {cell_types}')
    cell_type_proportions_df = major_composition_df.loc[:, cell_types]
    major_data1 = []
    for cell_type in cell_types:
        temp_rmpermute = []
        temp1 = []
        for C_ID in ordered_clusters:
            # spot in each COP
            spatial_cluster = [j for j in stC_dic['COP%d' % C_ID] if j in major_composition_df.index]
            true_mean = cell_type_proportions_df.loc[spatial_cluster, cell_type].mean()
            temp1.append(true_mean)
            differences = []
        major_data1.append(scipy.stats.zscore(temp1))

    # remove the random permutation(background)
    data_rmpermute = []
    for cell_type in cell_types:
        temp = []
        for C_ID in clusters:
            spatial_cluster = [j for j in stC_dic['COP%d' % C_ID] if j in immune_composition_df.index]
            ture_mean = cell_type_proportions_df.loc[spatial_cluster, cell_type].mean()
            differences = []
            for i in range(1000):
                shuffled_cell_type_proportions_df = shuffle(cell_type_proportions_df)
                shuffled_cell_type_proportions_df.index = cell_type_proportions_df.index
                permutated_mean = shuffled_cell_type_proportions_df.loc[spatial_cluster, cell_type].mean()
                difference = ture_mean - permutated_mean
                differences.append(difference)
            temp.append(np.mean(differences) / np.std(differences))
        data_rmpermute.append(temp)
    df_rmpermute = pd.DataFrame(data_rmpermute)
    df_rmpermute.columns = clusters
    df_rmpermute.index = cell_types



    major_data1 = major_data1     # or df_rmpermute.to_numpy()


    # +++++++++++++++++++++++++++++++++++++++++++++++++++++
    # plt
    clusters = [0, 1, 2, 3, 4, 5, 6, 7, 8]   # COPs
    
    clusters_color_map = {0: 'gold', 1: 'yellowgreen', 2: 'darkgreen', 3: 'orange', 4: 'blue', 5: 'deepskyblue',
                          6: 'salmon', 7: 'blueviolet', 8: 'lime', 9: 'turquoise', 10: 'peru', 11: 'green',
                          12: 'darkblue', 13: 'cyan', 14: 'indigo', 15: 'purple', 16: 'violet', 17: 'black',
                          18: 'fuchsia', 19: 'pink', 20: 'red'}
    
    plt.figure(figsize=(8, 10))
    
    gs = mpl.gridspec.GridSpec(10, 2, width_ratios=[15, 1], height_ratios=[2, 2, 2, 2, 2, 2, 2, 2, 2, 2], hspace=0.1)
    
    # COP and Pathological subtype
    ax = plt.subplot(gs[0, 0])
    left = 0
    for i in clusters:
        r = count_dic['COP%s' % i]['MIA'] / (count_dic['COP%s' % i]['MIA'] + count_dic['COP%s' % i]['IAC'])
        plt.barh(-0.5, r - 0.01, height=1, left=left + 0.01, color='blue')
        plt.barh(-0.5, 1 - r - 0.01, height=1, left=left + r, color='red')
        plt.barh(1, 1, height=1, left=left, color=clusters_color_map[i], edgecolor='black')
        left = left + 1
    
    plt.xlim(0, len(clusters))
    plt.yticks([-0.5, 1], ['Pathological subtype', 'COP'], size=20)
    plt.xticks([], [])
    
    # heatmap of mean composition
    vmin = -2
    vmax = 2
    cell_types_plt = ['Endothelial cells', 'AT1', 'AT2', 'Fibroblast',
                  'T/NK cells', 'B cells',
                  'Myeloid cells', 'Mast cells',
                  'Ciliated epithelial cells', 'Club cell', 'Malignent cells']
    ax = plt.subplot(gs[1:10, 0])
    
    plt.imshow(major_data1, cmap='coolwarm', vmin=vmin, vmax=vmax, aspect='auto')
    
    plt.yticks(list(range(len(cell_types_plt))), cell_types_plt, size=20)
    plt.xticks(list(range(len(clusters))), ['COP%d' % i for i in clusters], size=20, rotation=90)
    
    # legend
    ax = plt.subplot(gs[2, 1])
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb1 = mpl.colorbar.ColorbarBase(ax, cmap=mpl.cm.coolwarm, norm=norm)
    cb1.set_label('Relative abundance', size=20)
    cb1.set_ticks([-2, 2])
    cb1.set_ticklabels(['High', 'Low'])
    cb1.ax.tick_params(labelsize=15)
    
    ax = plt.subplot(gs[3, 1])
    
    patches = [plt.scatter([], [], marker='s', color='red', edgecolors='black', label='IAC'),
               plt.scatter([], [], marker='s', color='blue', edgecolors='black', label='MIA')]
    
    plt.legend(bbox_to_anchor=(3, -0.2), handles=patches, markerscale=2, fontsize=18, ncol=1)
    plt.axis('off')
    
    plt.show()







