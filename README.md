# Spatial coupling between malignant cell evolution and local immune microenvironment remodeling during LUAD pathological progression
## Joint analysis of single-cell and spatial transcriptomic data unveils the spatial coupling of malignant cell evolution and local microenvironment remodeling in early progression of lung adenocarcinoma

Spatial cellular organization patterns (COPs) within the tumor microenvironment (TME) evolve during tumor progression and influence patients' therapeutic responses. However, the cellular composition and functional potential of these multicellular structures during the progression of lung adenocarcinoma (LUAD) remain poorly understood. More importantly, the relationship between the evolution of malignant cells and the remodeling of the local immune microenvironment has not been thoroughly explored. In this study, we integrate spatial transcriptomics (ST) and single-cell RNA sequencing (scRNA-seq) data to characterize the local tumor and immunological landscape of tissue samples from eight patients with early-stage LUAD at various pathological stages. We identified eight major COPs that show distinct associations with pathological stages, local immune states, and clinical outcomes, including patient survival and therapy response. Notably, the infiltration of dysfunctional immune cells around malignant cells increases with pathological progression. Additionally, we explored cell-to-cell interactions between malignant cells and TME cells, revealing their role in remodeling the pro-tumor immune environment in invasive adenocarcinoma (IAC). Furthermore, integrating scRNA-seq, ST, and whole-exome sequencing (WES) data highlights the coupling between transcriptional alterations in malignant cells and the remodeling of the local microenvironment in IAC. This process leads to the establishment of a locally immune-suppressive environment, and our inferred copy number variations (CNVs) suggest that these genomic alterations drive the transcriptional evolution of malignant cells, contributing to the formation of a suppressive immune microenvironment in IAC. Overall, these findings enhance our understanding of the local microenvironment characteristics that contribute to LUAD progression and may help identify diagnostic markers and potential drug targets to inhibit invasive progression. 

![Analysis pipeline](https://github.com/haojiechen94/LUAD_COPs/blob/master/images/analysis_pipeline.png)
Deconvolution is applied to quantify the cellular composition in each spot and then unsupervised clustering is performed to integrate spots from different slides and identify cellular organization patterns (COPs), COPs are visualized in two dimensional UMAP. Finally, COPs are mapping into the pathologically annotated areas in HE images, by this mean to link the COPs to spatial structures. TAT, tumor adjacent tissues, TLSs, tertiary lymphoid structures, LAs, lymphoid aggregates.

![Neighboring_immune_cell_states](https://github.com/haojiechen94/LUAD_COPs/blob/master/images/quantifying_neighboring_immune_cells_and%20evaluating_local_immune_cell_states.png)
To quantify the infiltrating level of cell types (immune cell types) in each COPs, for a specific cell type, we first calculated the average proportion value as the observed average. Then we permuted the labels of spots and randomly selected a set of spots that matched the number of spots in the COP, computed the average proportion value as the permuted average, repeated this procedure 1,000 times. The differences between observed average and permuted average were computed and the mean value of the differences scaled by the its standard deviation was taken as the relative infiltrating score for each COP. Specifically, for a certain cell type, let $μ_{observed}$ denotes the average proportion value of this cell type and let $μ_{permuted}^{i}$ denotes the average proportion value in $i^{th}$ permutation. The relative infiltrating score is defined as follow:

$$δ_{i}=μ_{observed}-μ_{permuted}^{i}$$

$$relative\ infiltrating\ score=\frac{\frac{\sum_{i}δ_{i}}{1000}}{\sqrt{(δ_{i}-\frac{\sum_{i}δ_{i}}{1000})^2/1000}}$$


Positive scores indicate relative enrichment but negative scores represent relative depletion





Processed Data have been deposited in Zenodo (https://zenodo.org/record/8417887).

We are providing codes to reproduce our analysis results and tools for spatial transcriptomic analysis.

1. Part1_1_1_Identifying_cellular_organization_patterns_using_deconvolution_result_from_STs_and_scRNA_seq_data.R
2. Part1_1_2_Visualizing_deconvolution_result_and_COPs_on_the_top_of_tissue_histology.py
3. Part1_2_Linking_COPs_to_pathologically_annotated_areas_using_Fisher_exact_test.py
4. Part1_3_Quantifying_the_local_immune_cell_infiltration_level_in_normal_epithelal_or_malignant_cell_enriched_spot_group_or_COP.py
5. Part1_4_Identifying_significantly_enriched_pathways_in_each_COP.R
6. Part6_1_Categorizing_epithelial_cells_into subpopulations_using_the_in_house_analytical_pipelines.R
7. Part6_2_Inferring_CNVs_from _malignant_subtypes_scRNA_seq_data_using_inferCNV.R
8. Part6_3_Mapping_malignant_cell_subpopulations_onto_STs_by_using_RCTD.R
9. Part6_4_Analyzing_the_colocalization_relationship_between_malignant_cell_subpopulations_and_immune_cell_types.py
10. Part7_Identifying_CNVs_from_tumor_normal_matched_WES_data.sh

Citation: Haojie Chen et al. 2023. Spatiotemporal Profiling Unveiling the Cellular Organization Patterns and Local Protumoral Immune Microenvironment Remodeling in Early Lung Adenocarcinoma Progression. bioRxiv. https://doi.org/10.1101/2023.10.29.564580

<p align="center">
  <a href="#">
     <img src="https://api.visitorbadge.io/api/visitors?path=https://github.com/haojiechen94/LUAD_COPs" />
   </a>
</p>



