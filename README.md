# LUAD_COPs
## Spatiotemporal Profiling of Lung Adenocarcinoma Unveils the Microscopic Cellular Organization Patterns and Protumoral Immune Microenvironment Remodeling in Early Progression

Spatial cellular organization patterns (COPs) in tumor microenvironment influence the tumor progression and therapeutic response, however, little is known about the cellular composition and functional potential of these multicellular structures during lung adenocarcinoma progression. Here, we integrate spatial transcriptomics with single cell RNA sequencing to characterize the local tumor and immunological landscape of samples from 8 patients with early-stage lung adenocarcinoma at different pathological stages. We identified ten COPs that show distinct associations with local immune states and clinical outcomes, including survival and therapy response. The local infiltration levels of regulatory and dysfunctional immune cells are increased with pathological progression. Cell-to-cell interactions between malignant cells and tumor microenvironment (TME) cells were involved in protumor immune state remodeling. Finally, we detected a group of malignant cells that were specifically located at the tumor boundary, representing a more aggressive state, were involved in the invasion of invasive adenocarcinoma (IAC). Altogether, these results can improve our understanding of the local microenvironment characteristics that underlie LUAD progression and may facilitate the identification of drug targets to prevent invasive progression and biomarkers for diagnosis.

![Analysis pipeline](https://github.com/haojiechen94/LUAD_COPs/blob/master/images/analysis_pipeline.png)
Deconvolution is applied to quantify the cellular composition in each spot and then unsupervised clustering is performed to integrate spots from different slides and identify cellular organization patterns (COPs), COPs are visualized in two dimensional UMAP. Finally, COPs are mapping into the pathologically annotated areas in HE images, by this mean to link the COPs to spatial structures. TAT, tumor adjacent tissues, TLSs, tertiary lymphoid structures, LAs, lymphoid aggregates.

![Neighboring_immune_cell_states](https://github.com/haojiechen94/LUAD_COPs/blob/master/images/quantifying_neighboring_immune_cells_and%20evaluating_local_immune_cell_states.png)
To quantify the infiltrating level of cell types (immune cell types) in each COPs, for a specific cell type, we first calculated the average proportion value as the observed average. Then we permuted the labels of spots and randomly selected a set of spots that matched the number of spots in the COP, computed the average proportion value as the permuted average, repeated this procedure 1,000 times. The differences between observed average and permuted average were computed and the mean value of the differences scaled by the its standard deviation was taken as the relative infiltrating score for each COP. Specifically, for a certain cell type, let $μ_{observed}$ denotes the average proportion value of this cell type and let $μ_{permuted}^{i}$ denotes the average proportion value in $i^{th}$ permutation. The relative infiltrating score is defined as follow:

$δ_{i}=μ_{observed}-μ_{permuted}^{i}$

relative infiltrating score=$\sum_{i}$


Positive scores indicate relative enrichment but negative scores represent relative depletion





Processed Data have been deposited in Zenodo (https://zenodo.org/record/8417887).

We are providing codes to reproduce our analysis results and tools for spatial transcriptomic analysis.

1. Identifyiing_cellular_organization_patterns_by_integratiing_STs_and_scRNA_seq.R
2. Infering_the_composition_of_STs_spots_using_scRNA_seq_characterized_immune_cells.R
3. Computing_relative_celltype_enrichment_or_depletion_scores_in_each_COP.py
4. Spatial_heatmap.py
5. select_spots_on_HE.py
6. plot_spots_on_HE.py
7. Bubble_heatmap.py

Citation: Haojie Chen et al. 2023. Spatiotemporal Profiling Unveiling the Cellular Organization Patterns and Local Protumoral Immune Microenvironment Remodeling in Early Lung Adenocarcinoma Progression. bioRxiv. https://doi.org/10.1101/2023.10.29.564580
