# linking COPs to pathologically annotated areas by using fisher_exact to determine the degree of enrichment or depletion of annotated areas in each COP
# author: Xinyue Zhang & Haojie Chen
# date: 2025-03-08
# email: zhangxinyue2021@sinh.ac.cn & chenhaojie2017@sinh.ac.cn

import re
import pandas as pd
import glob
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib as mpl

from matplotlib import rcParams

rcParams['font.family']='Arial'
rcParams['pdf.fonttype']=42
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows',20)


def pvalue_to_zvalue(oddsratio,pvalue):
    z_value=0
    if oddsratio>1:
        z_value=scipy.stats.norm.ppf(1-pvalue/2)
        z_value=z_value if z_value!=np.inf else scipy.stats.norm.ppf(0.9999999999999999)
    else:
        z_value=scipy.stats.norm.ppf(pvalue/2)
        z_value=z_value if z_value!=np.inf else scipy.stats.norm.ppf(1-0.9999999999999999)
        z_value=z_value if z_value>scipy.stats.norm.ppf(1-0.9999999999999999) else scipy.stats.norm.ppf(1-0.9999999999999999)
    return z_value



if __name__ == '__main__':

    work_dir = 'data/ST/'


    # ++++++++++++++++++++++++++
    # COP label
    cop_meta_f = work_dir + 'CARD_deconvolution/major_cell_types/COPs_meta_data.txt'
    cop_meta_df = pd.read_csv(cop_meta_f, sep='\t')
    cluster_dic = {}
    for i in meta_data_df.index:
        ID = '_'.join(i.split('_')[:2])
        if ID != 'P12_GHY':
            barcode = i
            C_ID = meta_data_df.loc[i, 'seurat_clusters']
            if ID in cluster_dic:
                if C_ID in cluster_dic[ID]:
                    cluster_dic[ID][C_ID].add(barcode)
                else:
                    cluster_dic[ID][C_ID] = set([barcode])
            else:
                cluster_dic[ID] = {}
                cluster_dic[ID][C_ID] = set([barcode])


    # ++++++++++++++++++++++++++
    # T and NAT label
	dic = {'T': set()}
    for path in glob.glob(work_dir + 'hand_labeled_tumor_barcodes/*.txt'):
        # ID to ID_patient
        # P01 to P01_WYX
        ID = path.split('\\')[-1].split('_')[0]
        patient_ID = ID2patient_ID_map[ID]
        # whether to remove P13 with low quality
        if patient_ID != 'P13':
            with open(path) as infile:
                for line in infile:
                    # ID_patient_sbc
                    barcode = patient_ID + '_' + line.strip()
                    if barcode in meta_data_df.index:
                        dic['T'].add(barcode)
    dic['TAT'] = set([i for i in meta_data_df.index if 'P13' not in i]) - dic['T']



    # ++++++++++++++++++++++++++++++++++++++++
    # pathological annotation label
    dic1 = {}
    for i in dic['T']:
        patient_ID = '_'.join(i.split('_')[:2])
        if patient_ID in dic1:
            dic1[patient_ID]['T'].add(i)
        else:
            dic1[patient_ID] = {}
            dic1[patient_ID]['T'] = set([i])
            dic1[patient_ID]['TAT'] = set()
    for i in dic['TAT']:
        patient_ID = '_'.join(i.split('_')[:2])
        dic1[patient_ID]['TAT'].add(i)

	pathological_annotation_dic = {}
    for patient_ID in patientID_list:
        annotated_type_df = pd.read_csv(work_dir+'hand_labeled_histopathological_regions/%s.txt' % patient_ID,
            sep='\t')
        annotated_type_df.index = [patient_ID + '_' + i for i in annotated_type_df['barcode']]
        temp_dic = {}
        temp_dic['T'] = dic1[patient_ID]['T'] - set(annotated_type_df.index)
        temp_dic['TAT'] = dic1[patient_ID]['TAT'] - set(annotated_type_df.index)
        for i in set(annotated_type_df['annotated_type']):
            temp_dic[i] = set(annotated_type_df[annotated_type_df['annotated_type'] == i].index)
        pathological_annotation_dic[patient_ID] = temp_dic





	# ++++++++++++++++++++++++++++++++++++
	# patient label
	spots_in_each_patient_dic = {}
    for i in meta_data_df.index:
        ID = '_'.join(i.split('_')[:2])
        barcode = i
        if ID in spots_in_each_patient_dic:
            spots_in_each_patient_dic[ID].add(barcode)
        else:
            spots_in_each_patient_dic[ID] = set([barcode])




	# +++++++++++++++++++++++++++++++++++++
	# calculate the enrichment/depletion score of histopathological type for each COP
	# ['T', 'TAT', 'LA', 'coal_dust', 'vascular', 'TLSs', 'bronchiole']
	patientID_list
    C_IDs = list(set(meta_data_df['seurat_clusters'].tolist()))
    results_dic = {}
    for C_ID in C_IDs:
        results_dic[C_ID] = {}
        for ID in patientID_list:
            if C_ID in cluster_dic[ID]:
                results_dic[C_ID][ID] = {}
                for pathological_type in ['T', 'TAT', 'LA', 'coal_dust', 'vascular', 'TLSs', 'bronchiole']:
                    if pathological_type in pathological_annotation_dic[ID]:
                        oddsratio, pvalue = scipy.stats.fisher_exact([[
                            len(cluster_dic[ID][C_ID] & pathological_annotation_dic[ID][pathological_type]),
                            len((spots_in_each_patient_dic[ID] - cluster_dic[ID][C_ID]) &
                                pathological_annotation_dic[ID][pathological_type])],
                            [len((spots_in_each_patient_dic[ID] - pathological_annotation_dic[ID][pathological_type]) &
                                 cluster_dic[ID][C_ID]),
                             len((spots_in_each_patient_dic[ID] - cluster_dic[ID][C_ID]) &
                                 (spots_in_each_patient_dic[ID] - pathological_annotation_dic[ID][
                                     pathological_type]))]])
                        results_dic[C_ID][ID][pathological_type] = pvalue_to_zvalue(oddsratio, pvalue)
                    else:
                        results_dic[C_ID][ID][pathological_type] = 0
            else:
                results_dic[C_ID][ID] = {}
                for pathological_type in ['T', 'TAT', 'LA', 'coal_dust', 'vascular', 'TLSs', 'bronchiole']:
                    results_dic[C_ID][ID][pathological_type] = 0


	# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	# prop of each slide/patient for each COP
    proportions = []
    for C_ID in C_IDs:
        temp = []
        for patient_ID in patientID_list:
            if C_ID in cluster_dic[patient_ID]:
                temp.append(len(cluster_dic[patient_ID][C_ID]) /
                            len(spots_in_each_patient_dic[patient_ID]))
            else:
                temp.append(0)
        proportions.append(temp)





    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # zscore matrix 
    # for each pathological_type: patient x cluster
    data = []

    temp2 = []
    for C_ID in C_IDs:
        temp1 = []
        for i in patientID_list:
            temp1.append(results_dic[C_ID][i]['T'])
        temp2.append(temp1)
    data.append(temp2)

    temp2 = []
    for C_ID in C_IDs:
        temp1 = []
        for i in patientID_list:
            temp1.append(results_dic[C_ID][i]['TAT'])
        temp2.append(temp1)
    data.append(temp2)

    temp2 = []
    for C_ID in C_IDs:
        temp1 = []
        for i in patientID_list:
            if 'LA' in results_dic[C_ID][i]:
                temp1.append(results_dic[C_ID][i]['LA'])
            else:
                temp1.append(0)
        temp2.append(temp1)
    data.append(temp2)

    temp2 = []
    for C_ID in C_IDs:
        temp1 = []
        for i in patientID_list:
            if 'TLSs' in results_dic[C_ID][i]:
                temp1.append(results_dic[C_ID][i]['TLSs'])
            else:
                temp1.append(0)
        temp2.append(temp1)
    data.append(temp2)

    temp2 = []
    for C_ID in C_IDs:
        temp1 = []
        for i in patientID_list:
            if 'vascular' in results_dic[C_ID][i]:
                temp1.append(results_dic[C_ID][i]['vascular'])
            else:
                temp1.append(0)
        temp2.append(temp1)
    data.append(temp2)

    temp2 = []
    for C_ID in C_IDs:
        temp1 = []
        for i in patientID_list:
            if 'coal_dust' in results_dic[C_ID][i]:
                temp1.append(results_dic[C_ID][i]['coal_dust'])
            else:
                temp1.append(0)
        temp2.append(temp1)
    data.append(temp2)

    temp2 = []
    for C_ID in C_IDs:
        temp1 = []
        for i in patientID_list:
            if 'bronchiole' in results_dic[C_ID][i]:
                temp1.append(results_dic[C_ID][i]['bronchiole'])
            else:
                temp1.append(0)
        temp2.append(temp1)
    data.append(temp2)


    # plot 
    clusters_color_map = {0: 'gold', 1: 'yellowgreen', 2: 'darkgreen', 3: 'orange', 4: 'blue', 5: 'deepskyblue',
                          6: 'salmon', 7: 'blueviolet', 8: 'lime', 9: 'turquoise', 10: 'peru', 11: 'green',
                          12: 'darkblue', 13: 'cyan', 14: 'indigo', 15: 'purple', 16: 'violet', 17: 'black',
                          18: 'fuchsia', 19: 'pink', 20: 'red'}
    plt.figure(figsize=(10, 8))

    gs = mpl.gridspec.GridSpec(len(C_IDs), 9, width_ratios=[1] + [6] * 8, hspace=0, wspace=0)
    vmin = -8
    vmax = 8
    for i in range(len(C_IDs)):
        ax = plt.subplot(gs[i, 0])
        plt.barh(0, 1, left=-1.5, height=1, color=clusters_color_map[i])
        plt.yticks([0], ['COP%s' % C_IDs[i]], size=20)
        plt.xticks([], [])
        plt.ylim(-0.5, 0.5)

        ax = plt.subplot(gs[i, 1])
        if i == 0:
            plt.title('    T', size=20, rotation=90)
        plt.imshow([data[0][i]], cmap='bwr', vmin=vmin, vmax=vmax, aspect='auto')

        plt.xticks([], [])
        plt.yticks([], [])

        ax = plt.subplot(gs[i, 2])
        if i == 0:
            plt.title('    TAT', size=20, rotation=90)
        plt.imshow([data[1][i]], cmap='bwr', vmin=vmin, vmax=vmax, aspect='auto')
        plt.xticks([], [])
        plt.yticks([], [])

        ax = plt.subplot(gs[i, 3])
        if i == 0:
            plt.title('    LA', size=20, rotation=90)
        plt.imshow([data[2][i]], cmap='bwr', vmin=vmin, vmax=vmax, aspect='auto')
        plt.xticks([], [])
        plt.yticks([], [])

        ax = plt.subplot(gs[i, 4])
        if i == 0:
            plt.title('    TLSs', size=20, rotation=90)
        plt.imshow([data[3][i]], cmap='bwr', vmin=vmin, vmax=vmax, aspect='auto')
        plt.xticks([], [])
        plt.yticks([], [])

        ax = plt.subplot(gs[i, 5])
        if i == 0:
            plt.title('   Vascular', size=20, rotation=90)
        plt.imshow([data[4][i]], cmap='bwr', vmin=vmin, vmax=vmax, aspect='auto')
        plt.xticks([], [])
        plt.yticks([], [])

        ax = plt.subplot(gs[i, 6])
        if i == 0:
            plt.title('   Coal dust', size=20, rotation=90)
        plt.imshow([data[5][i]], cmap='bwr', vmin=vmin, vmax=vmax, aspect='auto')
        plt.xticks([], [])
        plt.yticks([], [])

        ax = plt.subplot(gs[i, 7])
        if i == 0:
            plt.title('  Bronchiole', size=20, rotation=90)
        plt.imshow([data[6][i]], cmap='bwr', vmin=vmin, vmax=vmax, aspect='auto')
        plt.xticks([], [])
        plt.yticks([], [])

        ax = plt.subplot(gs[i, 8])
        if i == 0:
            plt.title('  Whole slide', size=20, rotation=90)
        plt.imshow([proportions[i]], cmap='Reds', vmin=0, vmax=0.5, aspect='auto')
        plt.xticks([], [])
        plt.yticks([], [])


    cax = plt.axes([0.95, 0.6, 0.05, 0.1])
    norm = mpl.colors.Normalize(vmin=-6, vmax=6)
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap=mpl.cm.bwr, norm=norm)
    cb1.set_label('Relative enrichment', size=20)
    cb1.ax.tick_params(labelsize=15)

    cb1.set_ticks([-6, 6])
    cb1.ax.set_yticklabels(['Depleted', 'Enriched'], size=15)

    cax = plt.axes([0.95, 0.3, 0.05, 0.1])
    norm = mpl.colors.Normalize(vmin=0, vmax=0.5)
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap=mpl.cm.Reds, norm=norm)
    cb1.set_label('Proprotion', size=20)
    cb1.ax.tick_params(labelsize=15)

    cb1.set_ticks([0, 0.5])
    cb1.ax.set_yticklabels(['Low', 'High'], size=15)

    plt.show()




