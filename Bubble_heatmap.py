from matplotlib import rcParams
import pandas as pd
import numpy as np
import scipy.stats
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.pyplot as plt

rcParams['font.family']='Arial'
rcParams['pdf.fonttype']=42 

import pandas as pd
import scipy.stats
import os

def create_expression_profile(gene_expr_path,metadata_path,celltypes=None,marker_genes=None,group_by='celltype'):
    if os.path.isfile(gene_expr_path):
        gene_expr=pd.read_csv(gene_expr_path,sep='\t')
    else:
        print('Gene expression file dose not exist!')
        raise FileNotFoundError
        
    if os.path.isfile(metadata_path):
        metadata_df=pd.read_csv(metadata_path,sep='\t')
    else:
        print('Metadata file dose not exist!')
        raise FileNotFoundError
        
    if len(set(gene_expr.columns)&set(metadata_df.index))!=max(len(gene_expr.columns),len(metadata_df.index)):
        print('Ensure columns in Gene expression file and row in Metadata file are consistent!')
        raise
        
    if not celltypes:
        celltypes=list(set(metadata_df[group_by]))
    else:
        check=set(metadata_df[group_by])
        for celltype in celltypes:
            if celltype not in check:
                print('%s dose not exist in column %s of Metadata file'%(celltype,celltypes))
                raise
                
    if not marker_genes:
        marker_genes=list(gene_expr.index)
    else:
        check=set(gene_expr.index)
        for gene in marker_genes:
            if gene not in marker_genes:
                print('%s dose not record in Gene expression file'%(gene))
                raise        
        
    means=[]
    proportions=[]
    for gene in marker_genes:
        temp1=[]
        temp2=[]
        for celltype in celltypes:
            temp1.append(gene_expr.loc[gene,
                         metadata_df[metadata_df[group_by]==celltype].index].mean())
            temp2.append((gene_expr.loc[gene,
                         metadata_df[metadata_df[group_by]==celltype].index]!=0).mean()) 
        means.append(list(scipy.stats.zscore(temp1)))
        proportions.append(temp2) 
        
    means_df=pd.DataFrame(means,index=marker_genes,columns=celltypes)
    proportions_df=pd.DataFrame(proportions,index=marker_genes,columns=celltypes)
    return means_df,proportions_df

def bubble_heatmap(means_df,proportions_df,celltypes=None,marker_genes=None,vmin=-2,vmax=2,scale=300,cmap=mpl.cm.viridis):
    if celltypes:
        if len(set(celltypes)&set(means_df.columns))!=len(celltypes):
            raise
    else:
        celltypes=means_df.columns

    if marker_genes:
        if len(set(marker_genes)&set(means_df.index))!=len(marker_genes):
            raise        
    else:
        marker_genes=means_df.index
    
    plt.figure(figsize=(means_df.shape[0]/2,means_df.shape[1]/2))

    gs=mpl.gridspec.GridSpec(5,30,wspace=0.1)

    ax=plt.subplot(gs[:,:29])
    xs=[]
    ys=[]
    ss=[]
    cs=[]
    for i,g in enumerate(marker_genes):
        for j,c in enumerate(celltypes):
            xs.append(i)
            ys.append(j)
            ss.append(proportions_df.loc[g,c]*scale)
            cs.append(means_df.loc[g,c])
        
    plt.scatter(xs,ys,s=ss,c=cs,cmap=cmap,
                        vmin=vmin,vmax=vmax,edgecolors='black',
                        linewidths=0.1)

    plt.xticks([i for i in range(len(marker_genes))],marker_genes,size=20,rotation=90)
    plt.yticks([i for i in range(len(celltypes))],celltypes,size=20,rotation=0)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(which='major',length=8,width=2)
    ax.tick_params(which='minor',length=4,width=2)
    for axis in ['bottom','left']:
        ax.spines[axis].set_linewidth(2)

    ax=plt.subplot(gs[1,29])
    norm=mpl.colors.Normalize(vmin=vmin,
                              vmax=vmax)
    cb1=mpl.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm)
    cb1.set_label('Scaled expression\nlevel',size=12)
    cb1.ax.tick_params(labelsize=15)
    cb1.set_ticks([vmin,vmax])
    cb1.ax.set_yticklabels(['Low','High'],size=12)

    ax=plt.subplot(gs[3,29])
    patches=[plt.scatter([],[],marker='o',color='white',
                             edgecolors='black',label='25',s=scale*0.25),
            plt.scatter([],[],marker='o',color='white',
                             edgecolors='black',label='50',s=scale*0.50),
            plt.scatter([],[],marker='o',color='white',
                             edgecolors='black',label='75',s=scale*0.75),
            plt.scatter([],[],marker='o',color='white',
                             edgecolors='black',label='100',s=scale*1)]

    plt.legend(bbox_to_anchor=(6.5,0.5),handles=patches,markerscale=1,fontsize=12,
               title='%cells expressed',title_fontsize=12)
    plt.axis('off')

    plt.show()
    
bubble_heatmap(means_df,proportions_df,vmin=-2,vmax=2,scale=300,cmap=mpl.cm.viridis)