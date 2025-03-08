# inferring CNVs from scRNA data using inferCNV
# author: Xinyue Zhang & Haojie Chen
# date: 2025-03-08
# email: zhangxinyue2021@sinh.ac.cn & chenhaojie2017@sinh.ac.cn

# 
library(scater)
library(stringr)
library(readxl)
library(openxlsx)
library(dplyr)
library(patchwork)
library(ggrepel)
library(coda)
library(Seurat)
library(infercnv)
library(scCustomize)
library(ggthemes)



set.seed(1234)
setwd('/data/scRNAseq/')

# load data
seu <- readRDS('./malignant_ana/P12/P12_Epi_seu_obj.rds')
epi_cell_anno <- data.frame(cell.name = rownames(seu@meta.data),
                            cell.type = seu$anno)

gene_loc <- read.table('/data/gene_pos.txt', sep = '\t', header = F)
rownames(gene_loc) <- gene_loc$V1
gene_loc <- gene_loc[-1]


# gene expression matrix
Epithelia.UMI.counts<-read.table('./LUAD_scRNA_UMI_Epithelia_3618sc_19067gn.xls',sep='\t',header=T,row.names=1)
Epithelia.cells<-CreateSeuratObject(counts=Epithelia.UMI.counts, min.cells=3,min.features=200,project="Epithelia.cells")
Epi.anno.df<-read.table('./Epithelial_cells_annotations.txt',sep='\t',header=T,stringsAsFactors=F)
Epithelia.cells.subset.obj<-subset(Epithelia.cells,cells=row.names(Epi.anno.df))
donor.Epi.anno <- Epi.anno.df[which(Epi.anno.df$donor==donor),]
donor_exp <- as.matrix(Epithelia.cells.subset.obj@assays$RNA@counts[,rownames(donor.Epi.anno)])

all_raw_exp <- donor_exp[,rownames(all_cell_anno)]



# refercenced cell type: T and NK cells
immune_cell_anno <- read.table('./All_immune_cells_annotations_v2.txt', sep = '\t', header = T)
donor_immune_cell_anno <- immune_cell_anno[grep(pattern = donor, x=immune_cell_anno$cell.name),]
## T cell
T_subtype <- c('CD4+ Tcm', 'CD4+ Treg', 'Exhausted CD8+ T cell', 'GIMAP7+ CD8+ effector T cell', 
               'KLRB1+ CD8+ T cell','CD4+ Naive T cell', 'CXCL13+ CD4+ T cell','GZMK+ CD8+ effector T cell',
               'HSP1A1+ CD8+ effector T cell','XCL+ CD8+ T cell')
T_cell_anno <- donor_immune_cell_anno[which(donor_immune_cell_anno$sub.cell.type %in% T_subtype),]
T_cell_anno$sub.cell.type <- 'T cells'
## NK cell
NK_subtype <- c('GZMK+ NK', 'LTB+ NK', 'TNK','FCGR3A+ AREG+ NK','GIMAP7+ NK')
NK_cell_anno <- donor_immune_cell_anno[which(donor_immune_cell_anno$sub.cell.type %in% NK_subtype),]
NK_cell_anno$sub.cell.type <- 'NK cells'

ref_cell_anno <- rbind(T_cell_anno, NK_cell_anno)
colnames(ref_cell_anno) <- c('cell.name','cell.type')
rownames(ref_cell_anno) <- ref_cell_anno$cell.name

all_cell_anno <- rbind(epi_cell_anno, ref_cell_anno)
all_cell_anno <- all_cell_anno[-1]





# +++++++++++++++++++++++++++++++++++++++++++++++++++++
# infer CNV from scRNA data using inferCNV
dir.create('./malignant_ana/P12/inferCNV/')
# create inferCNV object
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = all_raw_exp,
                                     annotations_file = all_cell_anno,
                                     delim="\t",
                                     gene_order_file = gene_loc,
                                     min_max_counts_per_cell = c(100, +Inf),
                                     ref_group_names = c("NK cells", "T cells"))
# run inferCNV
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1,    
                              min_cells_per_gene = 3,  
                              window_length = 101,    
                              out_dir = './malignant_ana/P12/inferCNV/', 
                              cluster_by_groups = T,   
                              # k_obs_groups = 2,        
                              HMM = FALSE,
                              denoise = TRUE,
                              num_threads = 10,
                              scale_data = T,
                              analysis_mode ='samples'  
)




# +++++++++++++++++++++++++++++++++++++++++++++
# plot
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
library(ggpubr)
infercnv_obj.result <-  readRDS('./malignant_ana/P12/inferCNV/run.final.infercnv_obj')


# plot CNV pattern
## down sampling 300 ref cells
ref_cell <- ref_cell_anno$cell.name
downsample_300refcell <- sample(ref_cell, size=300, replace = F)
rm_refcell <- setdiff(ref_cell, downsample_300refcell)

expr <- infercnv_obj.result@expr.data
expr.choose <- expr[, rownames(all_cell_anno)]
sub_geneloc <- gene_loc[rownames(expr.choose),]

## cluster the CNV_score by k-means
kmeans.result <- kmeans(t(expr),9)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
## add kmeans_anno
all_anno <- merge(all_cell_anno, kmeans_df, by='row.names', all.x=T)
rownames(all_anno) <- all_anno$Row.names
all_anno$Row.names <- NULL
all_anno$cell.type <- NULL
all_anno$anno_all <- as.factor(all_anno$anno_all)
all_anno$anno_before <- as.factor(all_anno$anno_before)
all_anno$kmeans_class <- as.factor(all_anno$kmeans_class)
## order
all_anno <- arrange(all_anno, kmeans_class, anno)



top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v=RColorBrewer::brewer.pal(11, "Paired")[1:9]
names(color_v)=as.character(1:9)

left_anno <- rowAnnotation(df = all_anno, col=list(anno_all=c('AT1'='#C5E4A6', 'AT2'='#3FAE45','Club_cell'='#A24D56',
                                                              'Malignant_C1'='#EF7E33', 'Malignant_C2'='#B292F4', 'Malignant_C3'='#1A60D9',
                                                              'T cells' = '#EDD6D7', 'NK cells'='#F2D086'),
                                                   anno_before= c('AT1'='#C5E4A6', 'AT2'='#3FAE45','Club_cell'='#A24D56',
                                                                  'Mal_c1'='#EF7E33', 'Mal_c2'='#B292F4', 'Mal_c3'='#1A60D9',
                                                                  'Mal_c4'='#FFE101', 'T cells' = '#EDD6D7', 'NK cells'='#F2D086'),
                                                    kmeans_class=color_v),
                           border = TRUE,
                           simple_anno_size = unit(4, "mm"),
                           gap = unit(1,'mm'),
                           annotation_legend_param = list(anno_all = list(title_gp = gpar(fontsize = 12),
                                                                                 labels_gp = gpar(fontsize = 10)),
                                                          anno_before = list(title_gp = gpar(fontsize = 12),
                                                                                 labels_gp = gpar(fontsize = 10)),
                                                          kmeans_class = list(title_gp = gpar(fontsize = 12),
                                                                              labels_gp = gpar(fontsize = 10))
                           )
                           
)


tiff('./malignant_ana/P12/CNV_pattern.tiff', width = 40, height = 25, unit='cm', res = 800)
ht = Heatmap(t(expr.choose)[rownames(all_anno),], 
             col = colorRamp2(c(0.4,1,1.6), c("#377EB8","#F0F0F0","#E41A1C")), 
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneloc$V2, paste("chr",1:22,sep = "")), 
             column_gap = unit(2, "mm"),
             heatmap_legend_param = list(title = "Modified expression",
                                         direction = "vertical",
                                         title_position = "leftcenter-rot",
                                         at=c(0.4,1,1.6),
                                         legend_height = unit(3, "cm")),
             
             top_annotation = top_anno,left_annotation = left_anno, #添加注释
             row_title = NULL,column_title = NULL)
draw(ht, heatmap_legend_side = "right")
dev.off()




# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot and compare cnv level of each cell type or subtype
expr2 <- (expr-1)^2
CNV_score <- as.data.frame(colMeans(expr2))
colnames(CNV_score) <- 'cell.CNV_score'
CNV_score$cellbc <- rownames(CNV_score)

CNV_score <- merge(all_anno, CNV_score, by='row.names', all.x = T)
rownames(CNV_score) <- CNV_score$Row.names
CNV_score$Row.names <- NULL
CNV_score$anno_all <- as.factor(CNV_score$anno_all)
CNV_score$anno_before <- as.factor(CNV_score$anno_before)
CNV_score$kmeans_class <- as.factor(CNV_score$kmeans_class)
CNV_score$cellbc <- NULL


CNV_score <- CNV_score[which(CNV_score$anno_all!='NK cells' & CNV_score$anno_all!='T cells'),]
CNV_score %>% ggplot(aes(anno_all, cell.CNV_score))+
  geom_violin(aes(fill=anno_all), scale = 'width')+
  # geom_boxplot(width=0.05,position=position_dodge(0.9))+
  geom_signif(comparisons = list(c("Malignant_C1","AT2"),
                                 c("Malignant_C1","Malignant_C2"),c("Malignant_C1","Malignant_C3"),
                                 c("Malignant_C2","Malignant_C3")
                                 ),
              map_signif_level=F,
              textsize=3,test=wilcox.test,step_increase=0.08)+
  scale_fill_manual(values=c('AT1'='#C5E4A6', 'AT2'='#3FAE45','Club cells'='#A24D56',
                             'Malignant_C1'='#EF7E33', 'Malignant_C2'='#B292F4', 'Malignant_C3'='#1A60D9'))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 0.8))
