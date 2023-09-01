library(SPOTlight)
library(Seurat)
library(SeuratData)
library(dplyr)

set.seed(1234)

print('taking in scRNA-seq raw read count...')
#immune cells
TorNK.norm.counts<-read.table('../LUAD_scRNA_norm_TorNK_20712sc_20710gn.xls',sep='\t',header=T,row.names=1)
B.norm.counts<-read.table('../LUAD_scRNA_norm_B_1950sc_15446gn.xls',sep='\t',header=T,row.names=1)
Myeloid.norm.counts<-read.table('../LUAD_scRNA_norm_Myeloid_13416sc_20928gn.xls',sep='\t',header=T,row.names=1)
Mast.norm.counts<-read.table('../LUAD_scRNA_norm_Mast_841sc_12111gn.xls',sep='\t',header=T,row.names=1)

print('creating seurat objects...')
#seurat objects
TorNK.cells<-CreateSeuratObject(counts=TorNK.norm.counts,min.cells=3,min.genes=200,project="TorNK.cells")
B.cells<-CreateSeuratObject(counts=B.norm.counts,min.cells=3,min.genes=200,project="B.cells")
Myeloid.cells<-CreateSeuratObject(counts=Myeloid.norm.counts,min.cells=3,min.genes=200,project="Myeloid.cells")
Mast.cells<-CreateSeuratObject(counts=Mast.norm.counts,min.cells=3,min.genes=200,project="Mast.cells")

print('merging seurat objects...')
#merge objects
merge.obj<-merge(TorNK.cells,c(B.cells,Myeloid.cells,Mast.cells))
merge.obj@assays$RNA@data<-log2(as.matrix(merge.obj@assays$RNA@counts)+1)

print('taking in immune cell annotation information...')
#immune cell annotation data
annotation.df<-read.table('../all_immune_cells_annotations.txt',sep='\t',header=T,stringsAsFactors=F)
rownames(annotation.df)<-annotation.df$cell.name

print('subseting...')
#subseting
merge.obj<-subset(merge.obj,cells=annotation.df$cell.name)
merge.obj@meta.data$celltype<-annotation.df[rownames(merge.obj@meta.data),'sub.cell.type']
Seurat::Idents(object=merge.obj)<-merge.obj@meta.data$celltype

markers<-FindAllMarkers(merge.obj,only.pos=TRUE,min.pct=0.1,logfc.threshold=0.25)

subset.obj<-subset(merge.obj,cells=annotation.df2$cell.name)

P01<-Load10X_Spatial('../P01/',
                     filename='filtered_feature_bc_matrix.h5',
                     assay='Spatial')

spotlight_ls<-spotlight_deconvolution(se_sc=subset.obj,
                                      counts_spatial=P01@assays$Spatial@counts,
                                      clust_vr="celltype",
                                      cluster_markers=markers[c(markers$p_val_adj<0.01),],
                                      cl_n=100,
                                      hvg=3000,
                                      ntop=NULL,
                                      transf="uv",
                                      method="nsNMF",
                                      min_cont=0.01)

decon_mtrx<-spotlight_ls[[2]]
P01@meta.data<-cbind(P01@meta.data,decon_mtrx)

write.table(P01@meta.data,'../P01_immune_cell_type_proportions.txt',sep='\t',quote=F)

nmf_mod<-spotlight_ls[[1]]
h<-NMF::coef(nmf_mod[[1]])
rownames(h)<-paste("Topic",1:nrow(h),sep = "_")
write.table(h,'../P01_coefs.txt',sep='\t',quote=F)

w<-NMF::basis(nmf_mod[[1]])
colnames(w)<-paste("Topic",1:ncol(w),sep = "_")
write.table(w,'../P01_basis.txt',sep='\t',quote=F)

