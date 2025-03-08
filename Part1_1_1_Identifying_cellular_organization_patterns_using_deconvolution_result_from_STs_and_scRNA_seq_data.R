# deconvolution using CRAD
# identifying cellular organization patterns
# author: Xinyue Zhang & Haojie Chen
# date: 2025-03-08
# email: zhangxinyue2021@sinh.ac.cn & chenhaojie2017@sinh.ac.cn

library(CARD)
library(Seurat)
library(SeuratData)
library(dplyr)

set.seed(1234)

# +++++++++++++++++++++++++++++++++++++++++++++
# load data
print('taking in scRNA-seq norm read count...')
# immune cells
TorNK.norm.counts<-read.table('./data/scRNAseq/LUAD_scRNA_norm_TorNK_20712sc_20710gn.xls',sep='\t',header=T,row.names=1)
B.norm.counts<-read.table('./data/scRNAseq/LUAD_scRNA_norm_B_1950sc_15446gn.xls',sep='\t',header=T,row.names=1)
Myeloid.norm.counts<-read.table('./data/scRNAseq/LUAD_scRNA_norm_Myeloid_13416sc_20928gn.xls',sep='\t',header=T,row.names=1)
Mast.norm.counts<-read.table('./data/scRNAseq/LUAD_scRNA_norm_Mast_841sc_12111gn.xls',sep='\t',header=T,row.names=1)

# stromal cells
Ciliated.epithelia.norm.counts<-read.table('./data/scRNAseq/LUAD_scRNA_norm_Ciliated.epithelia_1297sc_16321gn.xls',sep='\t',header=T,row.names=1)
Endothelia.norm.counts<-read.table('./data/scRNAseq/LUAD_scRNA_norm_Endothelia_1690sc_15477gn.xls',sep='\t',header=T,row.names=1)
Fibroblast.norm.counts<-read.table('./data/scRNAseq/LUAD_scRNA_norm_Fibroblast_776sc_13551gn.xls',sep='\t',header=T,row.names=1)

# malignant cells
Epithelia.norm.counts<-read.table('./data/scRNAseq/LUAD_scRNA_norm_Epithelia_3618sc_19067gn.xls',sep='\t',header=T,row.names=1)
Epithelia.cells<-CreateSeuratObject(counts=Epithelia.norm.counts,min.cells=3,min.features=200,project="Epithelia.cells")
annotation.df<-read.table('./data/scRNAseq/epithelial_cells_annotations.txt',sep='\t',header=T,stringsAsFactors=F)

Epithelia.cells.subset.obj<-subset(Epithelia.cells,cells=row.names(annotation.df))
Epithelia.cells.subset.obj@meta.data$celltype<-annotation.df[rownames(Epithelia.cells.subset.obj@meta.data),'annotation_type']
Seurat::Idents(object=Epithelia.cells.subset.obj)<-Epithelia.cells.subset.obj@meta.data$celltype


# ++++++++++++++++++++++++++++++++++++++++++++
# create and merge seurat objects
print('creating seurat objects...')
# seurat objects
TorNK.cells<-CreateSeuratObject(counts=TorNK.norm.counts,min.cells=3,min.features=200,project="TorNK.cells")
B.cells<-CreateSeuratObject(counts=B.norm.counts,min.cells=3,min.features=200,project="B.cells")
Myeloid.cells<-CreateSeuratObject(counts=Myeloid.norm.counts,min.cells=3,min.features=200,project="Myeloid.cells")
Mast.cells<-CreateSeuratObject(counts=Mast.norm.counts,min.cells=3,min.features=200,project="Mast.cells")
Ciliated.epithelia.cells<-CreateSeuratObject(counts=Ciliated.epithelia.norm.counts,min.cells=3,min.features=200,project="Ciliated.epithelia.cells")
Endothelia.cells<-CreateSeuratObject(counts=Endothelia.norm.counts,min.cells=3,min.features=200,project="Endothelia.cells")
Fibroblast<-CreateSeuratObject(counts=Fibroblast.norm.counts,min.cells=3,min.features=200,project="Fibroblast")


print('merging seurat objects...')
# merge objects
merge.obj<-merge(TorNK.cells,c(B.cells,Myeloid.cells,Mast.cells,Ciliated.epithelia.cells,Endothelia.cells,Fibroblast,Epithelia.cells.subset.obj))
merge.obj<-JoinLayers(merge.obj)

merge.obj@meta.data[['tissue_type']]<-as.character(lapply(strsplit(colnames(merge.obj),'_'),function(x){x[4]}))
merge.obj@meta.data[['subtype']]<-as.character(lapply(strsplit(colnames(merge.obj),'_'),function(x){return(strsplit(x[8],'\\.')[[1]][1])}))
merge.obj@meta.data[['celltype']]<-c(rep("TorNK.cells",dim(TorNK.cells)[[2]]),rep("B.cells",dim(B.cells)[[2]]),
                                     rep("Myeloid.cells",dim(Myeloid.cells)[[2]]),rep("Mast.cells",dim(Mast.cells)[[2]]),
                                     rep("Ciliated.epithelia.cells",dim(Ciliated.epithelia.cells)[[2]]),
                                     rep("Endothelia.cells",dim(Endothelia.cells)[[2]]),
                                     rep("Fibroblast",dim(Fibroblast)[[2]]),Epithelia.cells.subset.obj@meta.data$celltype)
merge.obj@meta.data[['donor']]<-as.character(lapply(strsplit(colnames(merge.obj),'_'),function(x){return(paste(x[3],x[4],sep='_'))}))
unique(merge.obj@meta.data$donor)

sc_count<-merge.obj@assays$RNA$counts
sc_meta<-merge.obj@meta.data



# ++++++++++++++++++++++++++++++++++++++++++++++
# deconvolution using CARD(example: P01)
print('Deconvoluting a spot into a combination of major celltypes...')
P01 <- Load10X_Spatial('./data/ST/P01/',
                         filename='filtered_feature_bc_matrix.h5',
                         assay='Spatial')

SpatialFeaturePlot(P01, features = "nCount_Spatial")


spatial_count<-P01@assays$Spatial$counts[,P01$nCount_Spatial>1000]
spatial_location<-P01@images$slice1@boundaries$centroids@coords[P01$nCount_Spatial>1000,]
rownames(spatial_location)<-P01@images$slice1@boundaries$centroids@cells[P01$nCount_Spatial>1000]
colnames(spatial_location)<-c('y','x')
spatial_location<-as.data.frame(spatial_location)

CARD_obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = "celltype",
  ct.select = unique(sc_meta$celltype),
  sample.varname = NULL,
  minCountGene = 100,
  minCountSpot = 5) 

CARD_obj <-CARD_deconvolution(CARD_object = CARD_obj)

CARD.visualize.prop(
  proportion = CARD_obj@Proportion_CARD,        
  spatial_location = CARD_obj@spatial_location, 
  ct.visualize = c('TorNK.cells'),                 
  colors = c("lightblue","lightyellow","red"), 
  NumCols = 4,                                
  pointSize = 1.0)

write.table(CARD_obj@Proportion_CARD,'./data/ST/CARD_deconvolution/major_cell_types/P01_major_celltypes_composition.txt',sep='\t',quote=F)





# ++++++++++++++++++++++++++++++++++++++++++++
# identifying cellular organization patterns
print('identifying cellular organization patterns using standrad single cell RNA-seq analysis pipeline')
compositions<-read.table('./data/ST/CARD_deconvolution/major_cell_types/merged_major_celltypes_composition.txt',sep='\t',header=T,stringsAsFactors=F,check.names=F)
major.cell.types.compositions<-CreateSeuratObject(counts=compositions,
                                                  min.cells=1,min.genes=1,project="Major.cell.tpyes")
major.cell.types.compositions@assays$RNA@data<-major.cell.types.compositions@assays$RNA@counts

all.features<-rownames(major.cell.types.compositions)
major.cell.types.compositions<-ScaleData(major.cell.types.compositions,features=all.features)
major.cell.types.compositions<-RunPCA(major.cell.types.compositions,features=all.features)
major.cell.types.compositions<-FindNeighbors(major.cell.types.compositions,dims=1:7)
major.cell.types.compositions<-FindClusters(major.cell.types.compositions,resolution=0.25)
major.cell.types.compositions<-RunUMAP(major.cell.types.compositions,dims=1:7)

DimPlot(major.cell.types.compositions, label=T, group.by = 'seurat_clusters')
FeaturePlot(major.cell.types.compositions, features='B.cells', min.cutoff=0.02, max.cutoff=0.2)

write.table(major.cell.types.compositions@meta.data,
            './data/ST/CARD_deconvolution/major_cell_types/COPs_metadata.txt',
            sep='\t',quote=F)

