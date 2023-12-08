# select_spots_on_HE.py
# 2023-10-25
# Haojie Chen

"""
python select_spots_on_HE.py --HE_image=HE_image_path --tissue_position_list=tissue_position_list_csv_file 
                             [--scale_factor=<float>] [--outdir=output_directory]


--HE_image=<str>                    HE image in PNG format.

--tissue_position_list=<str>        Tissue position list in CSV format.

[--scale_factor=<float>]            Scale factor for HE image.
                                    Default: 1.0

[--outdir=<str>]                    Output directory for the processed result.
                                    Default: current directory

--help/-h                           print this page.

Description: This tool is used for selecting spots on HE image and save the spots list into files in TXT format.
"""

import numpy as np
from PIL import ImageDraw
from PIL import Image
from matplotlib.widgets import PolygonSelector
from matplotlib.path import Path
import matplotlib.pyplot as plt
import pandas as pd
from sys import argv, stderr, stdout
from getopt import getopt
import time


out_dir=False

class SelectFromCollection:
    def __init__(self, ax, collection, alpha_other=0.1):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Npts, 1))

        self.poly = PolygonSelector(ax, self.onselect)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.poly.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()


def toggle_selector(event):
    if event.key in ['A', 'a']:
        selected_spots=[]
        print('\nSelected points:...')
        print(len(selector.ind))
        for py,px in selector.xys[selector.ind]:
            selected_spots.append(tissue_positions_list_df[
                                           (tissue_positions_list_df['image_pixel_position_x']==px)&
                                           (tissue_positions_list_df['image_pixel_position_y']==py)].index[0])
        print(selected_spots[0],'...')
        with open('%s\\spots_%s.txt'%(out_dir,time.time()),'w') as outfile:
            for i in selected_spots:
                outfile.write(i+'\n')
        selected_spots=[]



if __name__ == '__main__':

    HE_image=''
    tissue_position_list=''
    scale_factor=1.0
    try:
        opts,args=getopt(argv[1:],'h',['HE_image=','tissue_position_list=','scale_factor=','outdir=','help'])
        for i,j in opts:   
            if i=="-h" or i=="--help":
                stdout.write(__doc__)
                exit(0)
            elif i=='--HE_image':
                HE_image=j
            elif i=='--outdir':
                out_dir=j
            elif i=='--tissue_position_list':
                tissue_position_list=j               
            elif i=='--scale_factor':
                scale_factor=float(j)
            else:
                raise Exception("Internal errors occur when parsing command line arguments.")
    except Exception as e:
        stderr.write("%s\n" % e)
        stderr.write("Type 'python select_spots_on_HE.py --help' for more information.\n")
        exit(1)

    if not out_dir:
        out_dir=os.getcwd()

    fig, axx = plt.subplots()
    fig.set_size_inches((100,100))
    img=Image.open(HE_image)
    plt.imshow(img)

    tissue_positions_list_df=pd.read_csv(tissue_position_list,
                                     names=['under_or_out_of_tissue',
                                            'array_position_x','array_position_y',
                                            'image_pixel_position_x',
                                            'image_pixel_position_y'])
    tissue_positions_list_df['image_pixel_position_x']=[i*scale_factor for i in tissue_positions_list_df['image_pixel_position_x']]
    tissue_positions_list_df['image_pixel_position_y']=[i*scale_factor for i in tissue_positions_list_df['image_pixel_position_y']]
    grid_x=[]    
    grid_y=[]                                    
    for i in tissue_positions_list_df.index:
        ax,ay,px,py,under_or_out_of_tissue=tissue_positions_list_df.loc[i,['array_position_x','array_position_y',
                                     'image_pixel_position_x','image_pixel_position_y','under_or_out_of_tissue']]
        if under_or_out_of_tissue==1:
            grid_x.append(py)
            grid_y.append(px)

    pts = axx.scatter(grid_x, grid_y,s=5,alpha=0)

    selector = SelectFromCollection(axx, pts)

    plt.connect('key_press_event',toggle_selector) 

    print("Select points in the figure by enclosing them within a polygon.")
    print("Press the 'esc' key to start a new polygon.")
    print("Try holding the 'shift' key to move all of the vertices.")
    print("Try holding the 'ctrl' key to move a single vertex.")
    print("Press the 'A/a' to save the selected area.")

    plt.axis('off')
    plt.show()

    selector.disconnect()




