# LUAD_COP
## Spatiotemporal Profiling Unveiling the Cellular Organization Patterns and Local Protumoral Immune Microenvironment Remodeling in Early Lung Adenocarcinoma Progression
![Analysis pipeline](https://github.com/haojiechen94/LUAD_COPs/blob/master/images/analysis_pipeline.png)
Deconvolution is applied to quantify the cellular composition in each spot and then unsupervised clustering is performed to integrate spots from different slides and identify cellular organization patterns (COPs), COPs are visualized in two dimensional UMAP. Finally, COPs are mapping into the pathologically annotated areas in HE images, by this mean to link the COPs to spatial structures. TAT, tumor adjacent tissues, TLSs, tertiary lymphoid structures, LAs, lymphoid aggregates.


Processed Data have been deposited in Zenodo (https://zenodo.org/record/8417887).

We are providing codes to reproduce our analysis results and tools for spatial transcriptomic analysis.

1. Identifyiing_cellular_organization_patterns_by_integratiing_STs_and_scRNA_seq.R
2. Infering_the_composition_of_STs_spots_using_scRNA_seq_characterized_immune_cells.R
3. Computing_relative_celltype_enrichment_or_depletion_scores_in_each_COP.py
4. Spatial_heatmap.py
5. select_spots_on_HE.py
6. plot_spots_on_HE.py
7. Bubble_heatmap.py


