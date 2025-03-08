# for each IAC patients, modeling heterogeneity using separate for epithelial cells to identify malignant subtypes
# author: Xinyue Zhang & Haojie Chen & Yana Li
# date: 2025-03-08
# email: zhangxinyue2021@sinh.ac.cn & chenhaojie2017@sinh.ac.cn & liyana2017@sinh.ac.cn

library(separate)
library(SingleCellExperiment)
library(scater)
library(stringr)
library(readxl)
library(openxlsx)
library(dplyr)
library(patchwork)
library(clustree)
library(ggrepel)
library(coda)
library(Seurat)


set.seed(1234)
setwd('/data/scRNAseq/')

# Epi_exp
Epithelia.UMI.counts<-read.table('./LUAD_scRNA_UMI_Epithelia_3618sc_19067gn.xls',sep='\t',header=T,row.names=1)
Epithelia.cells<-CreateSeuratObject(counts=Epithelia.UMI.counts, min.cells=3,min.features=200,project="Epithelia.cells")
Epi.anno.df<-read.table('./Epithelial_cells_annotations.txt',sep='\t',header=T,stringsAsFactors=F)
Epithelia.cells.subset.obj<-subset(Epithelia.cells,cells=row.names(Epi.anno.df))

# P12 for example 
donor = 'P12'
dir.create('./malignant_ana/P12/')

donor.Epi.anno <- Epi.anno.df[which(Epi.anno.df$donor==donor),]
donor_exp <- as.matrix(Epithelia.cells.subset.obj@assays$RNA@counts[,rownames(donor.Epi.anno)])

# ++++++++++++++++++++++++++++++++++++++++
# create separate obj
sce <- SingleCellExperiment(assays = list(counts = as.matrix(donor_exp)), colData = donor.Epi.anno)
sce <- filterGenes(sce, count_thres = 2, min_ncell = 5)
sce <- robustTransform(sce, method = rankProbit)
# gene GMM model
sce <- modelGenes(sce,
                  genes = rownames(sce),
                  if_assignments = TRUE,
                  n_cores = 8,    # parallel::detectCores() - 2,
                  verbose = TRUE
)

# select the features with adj.p<0.05
sce <- findFeatures(sce, n_features = 800, logit_flag = TRUE, parametricFit_flag = FALSE, test_family = "Normal") 
features.df <- metadata(sce)$separate@featuresInfo$findFeatures
features <- rownames(features.df)[features.df$if.feature]

# posterior probability
sce <- getChProbability(sce, genes = rownames(sce)) 
sce_prob_matrix <- metadata(sce)$separate@GMMInfo$ChProbability



# +++++++++++++++++++++++++++++++++++++++
# use posterior probability as scaled data to cluster the malignant cells
seu <- CreateSeuratObject(counts=donor_exp[rownames(sce_prob_matrix), colnames(sce_prob_matrix)], 
                          project = donor, 
                          meta.data = seu.meta)   
seu <- NormalizeData(seu, normalization.method = "LogNormalize")
seu@assays$RNA@scale.data <- as.matrix(metadata(sce)$separate@GMMInfo$ChProbability)
seu <- RunPCA(seu, verbose = F, features = features) # scale data
ElbowPlot(seu, ndims = 40)  

# cluster
seu <- FindNeighbors(seu, dims = 1:15)
seu <- FindClusters(object = seu, resolution = c(seq(0,1.5,0.1)))
cluster_tree <- clustree(seu@meta.data,prefix="RNA_snn_res.")+coord_flip()
seu <- FindClusters(seu, resolution = 0.7)
seu <- RunTSNE(seu, dims = 1:15)
seu <- RunUMAP(seu, dims = 1:15)
tsne_epi_p <- DimPlot(seu, reduction = 'tsne', group.by='Epi_anno', label=T, 
                      cols = c('AT1'='#C5E4A6', 'AT2'='#3FAE45','Club cells'='#A24D56',
                               'Malignant_C1'='#EF7E33', 'Malignant_C2'='#B292F4', 'Malignant_C3'='#1A60D9'))
tsne_p <- DimPlot(seu, reduction = "tsne",group.by = 'seurat_clusters',label = TRUE)

# relabel 
seu@meta.data <- seu@meta.data %>% mutate(anno = case_when(
  seurat_clusters==3 ~ 'AT1',
  seurat_clusters==5 ~ 'AT1',
  seurat_clusters==4 ~ 'AT2',
  seurat_clusters==6 ~ 'Club_cell',
  seurat_clusters==2 ~ 'Malignant_C1',
  seurat_clusters==1 ~ 'Malignant_C2',
  seurat_clusters==0 ~ 'Malignant_C3',
))





# +++++++++++++++++++++++++++++++++++++++++++++++++++++
# exp of immune inhibitory gene
seu %>% FeaturePlot_scCustom(features = c('CD74','NRP1','LGALS9'), reduction = 'tsne',
                             colors_use = viridis_plasma_light_high, 
                             # na_cutoff = NA, 
                             layer = 'data')




saveRDS(seu, './malignant_ana/P12/P12_Epi_seu_obj.rds')