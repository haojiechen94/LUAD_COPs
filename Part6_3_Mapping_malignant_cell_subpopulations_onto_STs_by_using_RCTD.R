# mapping malignant cell subpopulations onto STs by using RCTD to integrate the scRNA-seq and ST data for each sample
# 1. deconvolution by RCTD
# author: Xinyue Zhang & Haojie Chen
# date: 2025-03-08
# email: zhangxinyue2021@sinh.ac.cn & chenhaojie2017@sinh.ac.cn


library(scater)
library(stringr)
library(dplyr)
library(patchwork)
library(Seurat)
library(spacexr)


set.seed(1234)
setwd('/data/scRNAseq/')


# +++++++++++++++++++++++++++++++++++++++++++++++++++++
# deconvolution by RCTD

# load seurat obj of scRNA data(P12 for example)
seu <- readRDS('./malignant_ana/P12/P12_Epi_seu_obj.rds')

# create scRNA ref
select_cells <- rownames(seu@meta.data[which(seu$anno_before!='Club_cell'),])   # remove the celltype with few cells
seu_select <- subset(seu, cells=select_cells)  

sc_count <- seu_select@assays$RNA@counts
sc_celltype <- seu_select$anno
names(sc_celltype) <- rownames(seu_select@meta.data)
sc_celltype <- as.factor(sc_celltype)

# nUMI
nUMI <- seu_select$nCount_RNA
names(nUMI) <- rownames(seu_select@meta.data)

# ref obj
sc_ref <- Reference(sc_count, sc_celltype, nUMI)


# load ST data
P12_ST <- Load10X_Spatial('/data/ST/P12/',
                         filename='filtered_feature_bc_matrix.h5',
                         assay='Spatial')
# remove the spots with umi_count<=1000
spatial_count <- P12_ST@assays$Spatial$counts[,P12_ST$nCount_Spatial>1000]
spatial_location <- P12_ST@images$slice1@coordinates[P12_ST$nCount_Spatial>1000,]
spatial_location <- spatial_location[,2:3]
colnames(spatial_location) <- c("x","y")

# create spatial obj
coords <- spatial_location
nUMI <- colSums(spatial_count)
spatial.obj <- SpatialRNA(coords, spatial_count, nUMI)

# create RCTD obj
myRCTD <- create.RCTD(spatial.obj, sc_ref, max_cores = 10)

# run RCTD
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

# normalize
weight_matrix <- myRCTD@results$weights
nor_weight <- normalize_weights(weight_matrix)
nor_weight <- as.matrix(nor_weight)

malignant_subtypes <- grep(pattern='Malignant', unique(seu_select$anno), value=T)
mal_nor_weight <- nor_weight[, malignant_subtypes]
mal_nor_weight <- mal_nor_weight/rowSums(mal_nor_weight)

# save 
write.table(nor_weight, 
	'/data/inte_scRNA_ST/P12_malignant_subtypes_composition_RCTD.txt',
	sep='\t',quote=F)
