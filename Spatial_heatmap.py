import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.cm as cm
import glob
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image
import scipy.stats
from matplotlib import rcParams

rcParams['font.family']='Arial'
rcParams['pdf.fonttype']=42 


def spatial_color_map1(vmin,vmax):
    cmap=mpl.cm.viridis
    norm=mpl.colors.Normalize(vmin=vmin,
                              vmax=vmax)
    m=cm.ScalarMappable(norm=norm,cmap=cmap)
    return m,cmap,norm


#create position grid
def create_position_array():    
    dic={}
    for i in range(78):
        for j in range(64):
            if i%2==0:
                dic[(i,2*j)]=[j,i]
            else:
                dic[(i,2*j+1)]=[j+0.5,i]
    return dic

def get_tissue_potition_list(path):
	tissue_positions_list_df=pd.read_csv(path,sep=',',
                                         names=['under_or_out_of_tissue',
                                                'array_position_x','array_position_y',
                                                'image_pixel_position_x','image_pixel_position_y'])
	return tissue_positions_list_df

#spatial heatmap for gene expression levels and scores
#scores_series: Series, a gene expression/signature score vector 
#tissue_positions_list_df: tissue positions list from spaceranger, using get_tissue_potition_list to taking the file
#vmin, vmax:minimum and maximum value of colormap range
#name: name of gene or signature gene set
#path: output file path, e.g. ../test
def spatial_heatmap(scores_series,tissue_positions_list_df,borders,title='',name='',vmin=-0.5,vmax=0.5,path=''):
    m,cmap,norm=spatial_color_map1(vmin,vmax)
    dic=create_position_array()  
    plt.figure(figsize=(11,10))
    gs=mpl.gridspec.GridSpec(5,15)
    gs1=gs[:,:14]
    gs2=gs[1,14]
    gs3=gs[2,14]

    ax=plt.subplot(gs1)
    xs=[]
    ys=[]
    colors=[]
    for i in tissue_positions_list_df.index:
        x,y=tissue_positions_list_df.loc[i,['array_position_x','array_position_y']]
        xs.append(dic[(x,y)][0])
        ys.append(dic[(x,y)][1])
        if tissue_positions_list_df.loc[i,'under_or_out_of_tissue']==1:
            if i in scores_series.index:
                colors.append(m.to_rgba(scores_series[i]))
            else:
                colors.append(m.to_rgba(vmin))
        else:
            colors.append(m.to_rgba(vmin))
    
    plt.scatter(xs,ys,marker='h',s=115,c=colors,linewidths=0)
    
    xs=[]
    ys=[]
    colors=[]
    for i in tissue_positions_list_df.index:
        if i in borders:
            x,y=tissue_positions_list_df.loc[i,['array_position_x','array_position_y']]
            xs.append(dic[(x,y)][0])
            ys.append(dic[(x,y)][1])
            colors.append('white')
    
    plt.scatter(xs,ys,marker='o',s=30,c=colors,linewidths=0)    

    plt.xlim(-1,65)
    plt.ylim(78,-1)
    plt.title(title,size=20)
    plt.axis('off')

    ax=plt.subplot(gs2)
    cb1=mpl.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm)
    cb1.set_label(name,size=18)
    cb1.ax.tick_params(labelsize=15)
    
    
    ax=plt.subplot(gs3)
    patches=[plt.scatter([],[],marker='o',color='white',
                         edgecolors='black',label='Border',s=10)]

    plt.legend(bbox_to_anchor=(3.5,1),handles=patches,markerscale=3,fontsize=18)
    plt.axis('off')
    
    if path:
        plt.savefig('%s.pdf'%(path),bbox_inches='tight')
        plt.savefig('%s.tiff'%(path),bbox_inches='tight')

    plt.show()

