library(SPOTlight)
library(Seurat)
library(SeuratData)
library(dplyr)

#---------------------------------------------------------------------------------------------------------------------
#deconvoluion using SPOTlight
#---------------------------------------------------------------------------------------------------------------------

set.seed(1234)

print('taking in scRNA-seq norm read count...')
#immune cells
TorNK.norm.counts<-read.table('../LUAD_scRNA_norm_TorNK_20712sc_20710gn.xls',sep='\t',header=T,row.names=1)

B.norm.counts<-read.table('../LUAD_scRNA_norm_B_1950sc_15446gn.xls',sep='\t',header=T,row.names=1)

Myeloid.norm.counts<-read.table('../LUAD_scRNA_norm_Myeloid_13416sc_20928gn.xls',sep='\t',header=T,row.names=1)

Mast.norm.counts<-read.table('../LUAD_scRNA_norm_Mast_841sc_12111gn.xls',sep='\t',header=T,row.names=1)

#stromal cells
Ciliated.epithelia.norm.counts<-read.table('../LUAD_scRNA_norm_Ciliated.epithelia_1297sc_16321gn.xls',sep='\t',header=T,row.names=1)

Endothelia.norm.counts<-read.table('../LUAD_scRNA_norm_Endothelia_1690sc_15477gn.xls',sep='\t',header=T,row.names=1)

Fibroblast.norm.counts<-read.table('../LUAD_scRNA_norm_Fibroblast_776sc_13551gn.xls',sep='\t',header=T,row.names=1)

#epithelial/malignant cells
Epithelia.norm.counts<-read.table('..LUAD_scRNA_norm_Epithelia_3618sc_19067gn.xls',sep='\t',header=T,row.names=1)

print('creating seurat objects...')
#seurat objects
TorNK.cells<-CreateSeuratObject(counts=TorNK.norm.counts,min.cells=3,min.genes=200,project="TorNK.cells")
B.cells<-CreateSeuratObject(counts=B.norm.counts,min.cells=3,min.genes=200,project="B.cells")
Myeloid.cells<-CreateSeuratObject(counts=Myeloid.norm.counts,min.cells=3,min.genes=200,project="Myeloid.cells")
Mast.cells<-CreateSeuratObject(counts=Mast.norm.counts,min.cells=3,min.genes=200,project="Mast.cells")
Ciliated.epithelia.cells<-CreateSeuratObject(counts=Ciliated.epithelia.norm.counts,min.cells=3,min.genes=200,project="Ciliated.epithelia.cells")
Endothelia.cells<-CreateSeuratObject(counts=Endothelia.norm.counts,min.cells=3,min.genes=200,project="Endothelia.cells")
Fibroblast<-CreateSeuratObject(counts=Fibroblast.norm.counts,min.cells=3,min.genes=200,project="Fibroblast")
Epithelia.cells<-CreateSeuratObject(counts=Epithelia.norm.counts,min.cells=3,min.genes=200,project="Epithelia.cells")

annotation.df<-read.table('../epithelial_cells_annotations.txt',sep='\t',header=T,stringsAsFactors=F)

Epithelia.cells.subset.obj<-subset(Epithelia.cells,cells=row.names(annotation.df))

Epithelia.cells.subset.obj@meta.data$celltype<-annotation.df[rownames(Epithelia.cells.subset.obj@meta.data),'annotation_type']
Seurat::Idents(object=Epithelia.cells.subset.obj)<-Epithelia.cells.subset.obj@meta.data$celltype

print('merging seurat objects...')
#merge objects
merge.obj<-merge(TorNK.cells,c(B.cells,Myeloid.cells,Mast.cells,Ciliated.epithelia.cells,Endothelia.cells,Fibroblast,Epithelia.cells.subset.obj))

merge.obj@meta.data[['tissue_type']]<-as.character(lapply(strsplit(colnames(merge.obj),'_'),function(x){x[4]}))
merge.obj@meta.data[['subtype']]<-as.character(lapply(strsplit(colnames(merge.obj),'_'),function(x){return(strsplit(x[8],'\\.')[[1]][1])}))
merge.obj@meta.data[['celltype']]<-c(rep("TorNK.cells",dim(TorNK.cells)[[2]]),rep("B.cells",dim(B.cells)[[2]]),rep("Myeloid.cells",dim(Myeloid.cells)[[2]]),rep("Mast.cells",dim(Mast.cells)[[2]]),rep("Ciliated.epithelia.cells",dim(Ciliated.epithelia.cells)[[2]]),rep("Endothelia.cells",dim(Endothelia.cells)[[2]]),rep("Fibroblast",dim(Fibroblast)[[2]]),Epithelia.cells.subset.obj@meta.data$celltype)


print('identifing marker genes...')
merge.obj@assays$RNA@data<-log2(as.matrix(merge.obj@assays$RNA@counts)+1)

Seurat::Idents(object=merge.obj)<-merge.obj@meta.data$celltype
celltype.markers<-FindAllMarkers(object=merge.obj,assay="RNA",only.pos=TRUE)

print('UMAP visualization of major cell types...')
merge.obj<-FindVariableFeatures(merge.obj,selection.method="vst",nfeatures=2000)
merge.obj<-ScaleData(merge.obj)
merge.obj<-RunPCA(merge.obj)
merge.obj<-RunUMAP(merge.obj,dims=1:10)

DimPlot(merge.obj,group.by='celltype',label=TRUE)


print('Deconvoluting a spot into a combination of major celltyes...')
P01<-Load10X_Spatial('../P01/',
                      filename='filtered_feature_bc_matrix.h5',
                      assay='Spatial')

spotlight_ls<-spotlight_deconvolution(se_sc=merge.obj,
                                      counts_spatial=P01@assays$Spatial@counts,
                                      clust_vr="celltype",
                                      cluster_markers=celltype.markers[c(celltype.markers$p_val_adj<0.01),],
                                      cl_n=100,
                                      hvg=3000,
                                      ntop=NULL,
                                      transf="uv",
                                      method="nsNMF",
                                      min_cont=0.01)

decon_mtrx<-spotlight_ls[[2]]
P01@meta.data<-cbind(P01@meta.data,decon_mtrx)

write.table(P01@meta.data,'P01_major_cell_type_proportions.txt',sep='\t',quote=F)

nmf_mod<-spotlight_ls[[1]]
h<-NMF::coef(nmf_mod[[1]])
rownames(h)<-paste("Topic",1:nrow(h),sep = "_")
write.table(h,'P01_coefs.txt',sep='\t',quote=F)

w<-NMF::basis(nmf_mod[[1]])
colnames(w)<-paste("Topic",1:ncol(w),sep = "_")
write.table(w,'P01_basis.txt',sep='\t',quote=F)

#---------------------------------------------------------------------------------------------------------------------
#identifying cellular organization patterns using standrad single cell RNA-seq analysis pipeline
#---------------------------------------------------------------------------------------------------------------------

compositions<-read.table('../major_cell_types_composition.txt',sep='\t',header=T,stringsAsFactors=F,check.names=F)
major.cell.tpyes.compositions<-CreateSeuratObject(counts=compositions,
                                                  min.cells=1,min.genes=1,project="Major.cell.tpyes")
major.cell.tpyes.compositions@assays$RNA@data<-major.cell.tpyes.compositions@assays$RNA@counts

all.features<-rownames(major.cell.tpyes.compositions)
major.cell.tpyes.compositions<-ScaleData(major.cell.tpyes.compositions,features=all.features)
major.cell.tpyes.compositions<-RunPCA(major.cell.tpyes.compositions,features=all.features)
major.cell.tpyes.compositions<-FindNeighbors(major.cell.tpyes.compositions,dims=1:5)
major.cell.tpyes.compositions<-FindClusters(major.cell.tpyes.compositions,resolution=0.5)
major.cell.tpyes.compositions<-RunUMAP(major.cell.tpyes.compositions,dims=1:5)

DimPlot(major.cell.tpyes.compositions,label=TRUE)
FeaturePlot(major.cell.tpyes.compositions,features='B.cells',min.cutoff=0.02,max.cutoff=0.2)

write.table(major.cell.tpyes.compositions@meta.data,'../meta_data.txt',sep='\t',quote=F)





