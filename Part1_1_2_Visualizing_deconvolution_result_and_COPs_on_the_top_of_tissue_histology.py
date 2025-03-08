# visualizing the spots with the composition of celltypes on the top of tissue histology (HE image)
# visualizing the spots of each COPs on the top of tissue histology (HE image)
# author: Xinyue Zhang & Haojie Chen
# date: 2025-03-08
# email: zhangxinyue2021@sinh.ac.cn & chenhaojie2017@sinh.ac.cn


import glob
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec
from PIL import Image

rcParams['font.family'] = 'Arial'
rcParams['pdf.fonttype'] = 42
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows',20)


def plt_prop_onHE(
	deconvolution_f,
	tissue_positions_list_df,
	HE_image_path,
	border_spots):

    """Plot the spots with the composition of celltypes on the HE image 
    Args:
        deconvolution_f: the file storing the deconvolution result
        tissue_positions_list_df: the dataframe storing the position of spots
        HE_image_path: the path storing the HE image
        border_spots: the border spots
    """

    donor_major_prop = pd.read_csv(deconvolution_f, sep='\t', header=0, index_col=0)
    vmax = max(donor_major_prop.max())
    vmin = min(donor_major_prop.min())


    fig = plt.figure(figsize=(24, 18))
    gs = gridspec.GridSpec(3, 4)

    # spots to plt
    spots = [i.split('_')[-1] for i in list(donor_major_prop.index)]
    print(f'total num of spots: {len(spots)}')

    # celltype to plt
    anno = plt_anno  # 'AT1' for example

    grid_x = []
    grid_y = []
    color_str = []
    size = []
    for spot in spots:
        if spot in tissue_positions_list_df[tissue_positions_list_df['under_or_out_of_tissue'] == 1].index:
            ax, ay, px, py = tissue_positions_list_df.loc[spot, ['array_position_x', 'array_position_y',
                                                              'image_pixel_position_x',
                                                              'image_pixel_position_y']]

            grid_x.append(py)
            grid_y.append(px)
            # color_str.append(donor_major_prop.at[donor+'_'+spot, anno])
            color_str.append(donor_major_prop.at[spot, anno])

    plt_ax = plt.subplot(gs[row, col])
    # plt image
    img = Image.open(image_path)
    plt_ax.imshow(img)
    # plt the prop of each spot
    im = plt_ax.scatter(grid_x, grid_y, s=7, alpha=1, c=color_str, cmap='viridis',
                        # edgecolors='black',
                        linewidths=0.03)
    # plot border
    color = 'white'
    r = 2
    pxs = []
    pys = []
    for i in border_spots:
        if i in tissue_positions_list_df[tissue_positions_list_df['under_or_out_of_tissue'] == 1].index:
            ax, ay, px, py = tissue_positions_list_df.loc[i, ['array_position_x', 'array_position_y',
                                                              'image_pixel_position_x', 'image_pixel_position_y']]
            pxs.append(px)
            pys.append(py)

    im = plt_ax.scatter(pys, pxs, c=color, s=r,
                        # edgecolor='black',
                        marker='o', linewidth=0.03)

    plt_ax.axis('off')
    plt_ax.set_title(anno, fontsize=15)
    axes.append(plt_ax)

    # colorbar
    cb = fig.colorbar(im, ax=plt_ax, fraction=0.2, shrink=0.2, pad=0.05)
    cb.mappable.set_clim(vmin, vmax)
    colorbarfontdict = {"size": 15, 'family': 'Arial'}
    cb.ax.set_title('Proportion', fontdict=colorbarfontdict, pad=8)
    cb.ax.tick_params(labelsize=11, direction='in')

    plt.show()


if __name__ == '__main__':

    donor = plt_donor     # 'P01' for example

    st_dir = '/data/ST/'
    dic = {}
    scale_factor = 0.6666667
    
    # position
    tissue_positions_list_df = pd.read_csv(st_dir + f'{donor}/spatial/tissue_positions_list.csv',
                                               names=['under_or_out_of_tissue',
                                                      'array_position_x', 'array_position_y',
                                                      'image_pixel_position_x', 'image_pixel_position_y'])
    tissue_positions_list_df['image_pixel_position_x'] = [i * scale_factor for i in
                                                              tissue_positions_list_df['image_pixel_position_x']]
    tissue_positions_list_df['image_pixel_position_y'] = [i * scale_factor for i in
                                                              tissue_positions_list_df['image_pixel_position_y']]

    # HE image
    HE_image_path = st_dir + f'{donor}/spatial/tissue_hires_image.png'
    
    # border
    border_spots_f = st_dir + f'{donor}/spatial/border_spots/{donor}_border_barcodes.txt'
    border_spots_df = pd.read_csv(border_spots_f, sep='\t')
    border_spots = border_spots_df['barcode'].tolist()



    # plt prop on HE
    prop_plt_dir = '/data/ST/CARD_deconvolution/major_celltypes/HE_plt/pltProp_onHE/'
    if not os.path.isdir(prop_plt_dir):
        os.makedirs(prop_plt_dir)

    deconvolution_f = f'/data/ST/CARD_deconvolution/major_celltypes/{donor}_major_celltypes_composition.txt'

    plt_prop_onHE(
    	deconvolution_f,
    	tissue_positions_list_df,
    	HE_image_path,
    	border_spots)

    