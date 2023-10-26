# select_spots_on_HE.py
# 2023-10-25
# Haojie Chen

"""
python plot_spots_on_HE.py --HE_image=HE_image_path --tissue_position_list=tissue_position_list_csv_file 
                           --spots=spots_file


--HE_image=<str>                    HE image in PNG format.

--tissue_position_list=<str>        Tissue position list in CSV format.

--spots=<str>                       Selected spot list genereated from select_spot_on_HE.py tool (TXT format).

--help/-h                           print this page.

Description: This tool is used for showing spots on HE image.
"""

import numpy as np
from PIL import ImageDraw
from PIL import Image
import matplotlib.pyplot as plt
import pandas as pd
from sys import argv, stderr, stdout
from getopt import getopt


if __name__ == '__main__':

    HE_image=''
    tissue_position_list=''
    spots=''
    try:
        opts,args=getopt(argv[1:],'h',['HE_image=','tissue_position_list=','spots=','help'])
        for i,j in opts:   
            if i=="-h" or i=="--help":
                stdout.write(__doc__)
                exit(0)
            elif i=='--HE_image':
                HE_image=j
            elif i=='--spots':
                spots=j
            elif i=='--tissue_position_list':
                tissue_position_list=j               
            else:
                raise Exception("Internal errors occur when parsing command line arguments.")
    except Exception as e:
        stderr.write("%s\n" % e)
        stderr.write("Type 'python plot_spots_on_HE.py --help' for more information.\n")
        exit(1)


    fig, axx = plt.subplots()
    fig.set_size_inches((100,100))

    fig, axx = plt.subplots()
    fig.set_size_inches((20,20))

    selected_spots=[]
    with open(spots) as infile:
        for line in infile:
            selected_spots.append(line.strip())    
    #'C:\\Users\\chenhaojie\\Desktop\\STs_and_scRNA_seq\\Datasets\\spatial_infos\\P02_JGZ.jpg'
    img=Image.open(HE_image)
    plt.imshow(img)

    #'C:\\Users\\chenhaojie\\Desktop\\STs_and_scRNA_seq\\Datasets\\spatial_infos\\P02_JGZ_tissue_positions_list.csv'
    tissue_positions_list_df=pd.read_csv(tissue_position_list,
                                     names=['under_or_out_of_tissue',
                                            'array_position_x','array_position_y',
                                            'image_pixel_position_x',
                                            'image_pixel_position_y'])
    grid_x=[]    
    grid_y=[]                                    
    for i in selected_spots:
        ax,ay,px,py,under_or_out_of_tissue=tissue_positions_list_df.loc[i,['array_position_x','array_position_y',
                                     'image_pixel_position_x','image_pixel_position_y','under_or_out_of_tissue']]
        if under_or_out_of_tissue==1:
            grid_x.append(py)
            grid_y.append(px)

    pts = axx.scatter(grid_x, grid_y,s=5,alpha=1)

    plt.axis('off')
    plt.show()
