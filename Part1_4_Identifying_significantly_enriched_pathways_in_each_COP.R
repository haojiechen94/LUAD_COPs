# identifying significantly enriched pathways of each COP using ssGSEA+limma
# author: Xinyue Zhang & Haojie Chen
# date: 2025-03-08
# email: zhangxinyue2021@sinh.ac.cn & chenhaojie2017@sinh.ac.cn

library(GSEABase)
library(Seurat)
library(GSVA)
library(GSVAdata)
library(dplyr)
library(stringr)
library(gson)
library(msigdbr)
library(tidyverse)
library(BiocParallel)
library(patchwork)
library(ggpubr)

library(limma)
library(tibble)
library(stringr)

set.seed(0909)


# ++++++++++++++++++++++++++++++++++++++++++++++++
# load expr mtx and genesets

print('load the normalize expr of spots from all samples')
# load the normalize expr of spots from all patients
st.merge.obj <- readRDS('./data/ST/merge_all_spots_obj.rds')
st.combined.nor_data <- st.merge.obj@assays$SCT@data %>% as.matrix()
dim(st.combined.nor_data)

# C2 genesets
all_c2_genesets <- read.gmt('./data/c2.all.v2024.1.Hs.symbols.gmt')
all_c2_geneset_df <- all_c2_genesets %>%
  group_by(term) %>%
  summarise(gene_list = list(gene)) %>% as.data.frame()
rownames(all_c2_geneset_df) <- all_c2_geneset_df$term
c2_geneset_list <- as.list(all_c2_geneset_df$gene_list)
names(c2_geneset_list) <- rownames(all_c2_geneset_df)
# select KEGG, REACTOME, BIOCARTA
filtered_terms.2 <- as.character(all_c2_geneset_df[grepl(pattern = 'REACTOME|BIOCARTA|KEGG', 
                                                         all_c2_geneset_df$term),]$term)
filtered.c2_geneset_list <- as.list(all_c2_geneset_df[filtered_terms.2,]$gene_list)
names(filtered.c2_geneset_list) <- rownames(all_c2_geneset_df[filtered_terms.2,])

# hallmark geneset
hallmark_term_gene <- read.gmt('./h.all.v2023.2.Hs.symbols.gmt')
hallmark_geneset_df <- hallmark_term_gene %>%
  group_by(term) %>%
  summarise(gene_list = list(gene))
hallmark_geneset_df$term <- as.character(hallmark_geneset_df$term)

hallmark_geneset_list <- as.list(hallmark_geneset_df$gene_list)
names(hallmark_geneset_list) <- hallmark_geneset_df$term

# combine two geneset
combines.geneset <- c(hallmark_geneset_list, filtered.c2_geneset_list)



# +++++++++++++++++++++++++++++++++++++++++++++
# ssGSEA for each sample(patient)
print('ssGSEA for each patient')
st.combined.params <- ssgseaParam(exprData = st.combined.nor_data[, p.spots], 
                                    geneSets = combines.geneset,
                                    minSize = 2, maxSize = Inf,
                                    normalize = TRUE)
# run gsva
st.combined.gsva.result <- gsva(st.combined.params,
                                  verbose=T,
                                  BPPARAM = MulticoreParam(workers = 40, progressbar = T))
st.combined.gsva.score <- as.data.frame(st.combined.gsva.result)
write.table(st.combined.gsva.result, './data/ST/ssGSEA_COP/P01_st_ssGSEA_filteredC2hallmark.txt',  sep='\t', row.names = T, col.names = T, quote = F)





# +++++++++++++++++++++++++++++++++++++++++
# identifying significantly enriched pathways of each COP using limma

# load ssGSEA result
merged.ssgsea <- read.table('./data/ST/ssGSEA_COP/merged_st_ssGSEA_filteredC2hallmark.txt', 
                            sep='\t', header = T, check.names =F)
colnames(merged.ssgsea) <- gsub('\\.','-',colnames(merged.ssgsea))

# load COPs meta data
COPs_meta <- read.csv('./data/ST/CARD_deconvolution/major_cell_types/COPs_metadata.txt',
                      sep='\t')
table(COPs_meta$COPs)

# limma
rownames(COPs_meta) <- str_split_fixed(rownames(COPs_meta), '_', 2)[,2]
spots <- intersect(rownames(COPs_meta), colnames(merged.ssgsea))
COPs_meta <- COPs_meta[spots,]
gsva.score <- merged.ssgsea[, spots]

st.combined.top_col <- COPs_meta[,c('COPs'), drop=F]
st.combined.top_col <- arrange(st.combined.top_col, COPs)
cops <- names(table(COPs_meta$COPs))

# initialize
diff_all <- data.frame()

for (focus_cop in cops){
  other_cops <- 'others'
  # contrast pairs
  pair <- paste(focus_cop, other_cops, sep='-')
  print(pair)
  # group col
  group_col <- COPs_meta[,c('COPs'), drop=F]
  group_col$COP <- ifelse(group_col$COPs==focus_cop, focus_cop, other_cops)
  group_col$COPs <- NULL
  colnames(group_col) <- c('COPs')
  group_col <- arrange(group_col, desc(COPs))  # others and focus_cop
  # rownames(group_col) <- str_split_fixed(rownames(group_col),'_',n=2)[,2]
  
  # design
  DE_design <- model.matrix(~0 + group_col$COPs)
  colnames(DE_design) <- str_split_fixed(colnames(DE_design), 'COPs', n=2)[,2]
  rownames(DE_design) <- rownames(group_col)

  # contrast
  contrasts.matrix <- makeContrasts(contrasts = pair, levels = DE_design)
  fit <- lmFit(gsva.score[,rownames(group_col)], DE_design)
  fit2 <- contrasts.fit(fit, contrasts = contrasts.matrix)
  fit2 <- eBayes(fit2)
  
  # result
  diff <- topTable(fit2, coef = 1, n=Inf, adjust.method = "BH", sort.by = "logFC")  
  diff <- diff[which(diff$adj.P.Val<0.05),]
  diff <- rownames_to_column(diff, 'pathways')
  diff <- diff %>% mutate(group = case_when(
    diff$logFC>0 ~ 'Up',
    diff$logFC<0 ~ 'Down'
  ))
  diff <- arrange(diff, logFC)
  diff$pathways <- factor(diff$pathways, levels = diff$pathways)
  diff$type <- pair
  print(nrow(diff))
  diff_all <- rbind(diff_all, diff)
}

write.table(diff_all, './data/ST/ssGSEA_COP/COPenriched_sig_pathways.txt', sep='\t', quote = F,
            row.names = F, col.names = T)
