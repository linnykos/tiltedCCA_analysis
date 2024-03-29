# Purpose

This repository contains all the preprocessing and analyses for the paper "Quantifying common and distinct information in single-cell multimodal data with Tilted-CCA", in order to illustrate Tilted-CCA. See the companion GitHub package https://github.com/linnykos/tiltedCCA, which contains the primary R functions for the standalone method.

This code was developed and tested primarily on R 4.1.2. on a Macbook (macOS 11.6.8 Big Sur) equipped with an i7 processor.

# Dependencies to run all code

The data analysis requires several additional packages in order to run all the analyses, in addition to those needed to install `tiltedCCA`. These include `ArchR`, `igraph`, `MAST`, `MOFA2`, `SeuratData`, `Signac`, and `slingshot`. See the 
last section of this README to see where (i.e., CRAN, Bioconductor, or GitHub) to download all such packages.

# Small simulated dataset to demo the software

See https://github.com/linnykos/tiltedCCA_analysis/tree/master/simulation for the small demo on how to use Tilted-CCA.

# Data

All data used in our paper publicly available. 
* The human bone marrow CITE-seq data[1] is available in the R package `SeuratData` under `bmcite`.
* The human bone marrow Abseq data[2] (both, either equipped with the whole transcriptome or of only 461 genes) is availabe at https://figshare.com/projects/Single-cell_proteo-genomic_reference_maps_of_the_human_hematopoietic_system/94469, under `WTA_projected.rds` or `Healthy.rds` respectively. 
* The human PBMC CITE-seq data[3] is available at https://atlas.fredhutch.org/nygc/multimodal-pbmc/.
* The human PBMC 10x Multiome data is available at https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k. 
* The human brain development 10x Multiome data[4] is available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162170 and https://github.com/GreenleafLab/brainchromatin.
* The mouse brain development 10x Multiome data is available at https://www.10xgenomics.com/resources/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-2-0-0.
* The mouse embryo development 10x Multiome data[5]
is available at https://github.com/rargelaguet/mouse_organogenesis_10x_multiome_publication. 

The majority of the computational tools used in this manuscript consisted of Seurat (v4.1.1), Signac (v1.7.0), and Slingshot (v2.3.1). We derived the cell-cycle genes for humans via `Seurat::cc.genes`, while we used the lists of human housekeeping genes[6] or mouse cell-cycling genes[7] from particular published work.

# Reproducing the results

All the code to reproduce the analyses in the paper are located in the `main` folder. We first enumerate all the code to preprocess the data, and then we enumerate the code to reproduce the figures. When reproducing the code, the specific filepaths and filenames will differ from your setup.

## Preprocessing the data

* Human bone marrow CITE-seq: See https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_preprocess.R
* Human bone marrow Abseq (Full transcriptome): See https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_97antibody_preprocess.R
* Human bone marrow Abseq (461 genes): See https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/abseq_bm_97antibodyRef_preprocess.R
* Human PBMC CITE-seq: See https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/pbmc_224antibody_preprocess.R
* Human PBMC 10x Multiome: See https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/pbmc_10x_preprocess.R
* Human brain development 10x Multiome: See https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_preprocess.R and https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_preprocess_customGeneActivity.R
* Mouse brain development 10x Multiome: See https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/mouseembryo_preprocess.R
* Mouse embryo development 10x Multiome: See https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/reik_preprocess.R

## Reproducing the figures

* Figure 1a,f: To run Tilted-CCA on the human bone marrow CITE-seq data, see https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_tiltedcca.R. The UMAPs for the RNA, protein and Tilted-CCA are plotted via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_tiltedcca_postprocess.R. 

* Figure 1c: The Consensus PCA results and plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_consensusPca.R

* Figure 1d: After running Tilted-CCA on the human bone marrow CITE-seq data, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_2dplots.R 

* Figure 2a: The results and plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_tiltedcca_extremeExamples.R

* Figure 3a: After running Tilted-CCA on the human bone marrow CITE-seq data, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_tiltedcca_postprocess.R.  

* Figure 3b: After running Tilted-CCA on the human bone marrow CITE-seq data, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_cell_enrichment.R. 

* Figure 3c: After running Tilted-CCA on the human bone marrow CITE-seq data, the differential-expression for genes is computed in https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_differential.R. Then, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_differential_postprocess.R

* Figure 3e (Human PBMC CITE-seq data): To run Tilted-CCA see https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/pbmc_224antibody_tiltedcca-2000genes.R. The differential-expression for genes and antibodies is computed in https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/pbmc_224antibody_differential.R. Then, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/pbmc_224antibody_differential_postprocess.R

* Figure 3e (Human PBMC 10x Multiome): To run Tilted-CCA see https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/pbmc_10x_tiltedcca-2000genes.R. The differential-expression for genes is computed in https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/pbmc_10x_differential.R. Then, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/pbmc_10x_differential-2000genes_postprocess.R

* Figure 3e (Human brain development 10x Multiome): To run Tilted-CCA (RNA+ATAC) see https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_tiltedcca_RNA-ATAC.R. The differential-expression for genes is computed in https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_differential.R. Then, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_differential_postprocess.R

* Figure 4a: The plots were all created during preprocessing (see https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/abseq_bm_97antibodyRef_preprocess.R).

* Figure 4b: To run Tilted-CCA on the human bone marrow Abseq data (461 genes), see https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/abseq_bm_97antibodyRef_tcca.R. Each antibody's differential-expression of the distinct component is computed via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/abseq_bm_97antibodyRef_differential-distinct.R. Then, both results are used when selecting the 10 antibodies and plots are made via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/abseq_bm_97antibodyRef_antibodySelection.R

* Figure 4c: After selecting the antibodies for the human bone marrow Abseq data (461 genes), the heatmap is plotted via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/abseq_bm_97antibodyRef_antibodySelection_heatmap.R

* Figure 4d: After selecting the antibodies for the human bone marrow Abseq data (461 genes), the heatmap is plotted via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/abseq_bm_97antibodyRef_antibodySelection_graph.R

* Figure 4e: The alternative antibody panels for the human bone marrow Abseq data (461 genes) are computed via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/abseq_bm_97antibodyRef_antibodySelection_alternatives.R. These results are combined with the selected antibodies via Tilted-CCA and plotted via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/abseq_bm_97antibodyRef_antibodySelection_enrichment.R

* Figure 5a,f,g: To run Tilted-CCA (RNA+ATAC) on the human brain development 10x Multiome data, see https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_tiltedcca_RNA-ATAC.R. The synchrony scores are computed and plotted via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_steadystate.R

* Figure 5b,c: To run Tilted-CCA on the mouse embryo development 10x Multiome data, see https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/reik_tiltedcca.R. The synchrony scores are computed and plotted via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/reik_steadystate.R

* Figure 5d,e: To run Tilted-CCA (RNA+ATAC) on the mouse brain development 10x Multiome data, see https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/mouseembryo_tiltedcca_RNA-geneActivity.R. The synchrony scores are computed and plotted via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/mouseembryo_steadystate.R

* Figure 5j: To run Tilted-CCA (RNA+ChromatinActivity) on the mouse brain development 10x Multiome data, see https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/mouseembryo_tiltedcca_RNA-geneActivity.R. The developmentally-informative genes are selected via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/mouseembryo_developmentalGenes.R, and the heatmap is plotted via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/mouseembryo_developmentalGenes_heatmap.R

* Figure 5k: To run Tilted-CCA (RNA+ChromatinActivity) on the human brain development 10x Multiome data, see https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_tiltedcca_RNA-geneActivity.R. The developmentally-informative genes are selected via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_developmentalGenes.R, and the heatmap is plotted via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_developmentalGenes_heatmap.R

* Figure 6: After running Tilted-CCA on the human brain development 10x Multiome data, the phase portraits are computed and plotted via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_leafplots.R. 

* Extended Figure 1a: Consensus PCA is applied on the human bone marrow CITE-seq data via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_consensusPca.R

* Extended Figure 1b: JIVE is applied on the human bone marrow CITE-seq data via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_jive.R

* Extended Figure 1c: scAI is applied on the human bone marrow CITE-seq data via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_scai.R

* Extended Figure 1d: MOFA+ is applied on the human bone marrow CITE-seq data via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_mofa.R

* Extended Figure 1e,f: WNN is applied on the human bone marrow CITE-seq data during preprocessing, see https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_preprocess.R

* Extended Figure 1g: DCCA is applied on the human bone marrow CITE-seq data via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_dcca.R

* Extended Figure 2a: After running Tilted-CCA on the human bone marrow CITE-seq data, the cell enrichment plots are constructed via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_cell_enrichment.R

* Extended Figure 2b: After running Consensus PCA and Tilted-CCA on the human bone marrow CITE-seq data, the cell enrichment plots are constructed via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_enrichment_consensus-tcca.R

* Extended Figure 2c: These alignment-differentiability plots are made in the same script that produces Figure 3c

* Extended Figure 2d: The violin plots of the ATAC libraries is plotted via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/cell-gene_library_violinplots.R

* Extended Figure 2e (Human PBMC 10x Multiome): The plot is made in the same script that produces Figure 3e (Human PBMC 10x Multiome)

* Extended Figure 2e (Mouse brain development 10x Multiome): After running Tilted-CCA on the mouse brain development 10x Multiome data, the differential-expression for genes is computed in https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/mouseembryo_differential.R. Then, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/mouseembryo_differential_postprocess.R

* Extended Figure 2e (Mouse embryo development 10x Multiome): After running Tilted-CCA on the mouse embryo development 10x Multiome data, the differential-expression for genes is computed in https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/reik_differential.R. Then, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/reik_differential_postprocess.R

* Extended Figure 3a: After running Tilted-CCA on the human bone marrow Abseq data (whole transcriptome), the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_97antibody_tiltedcca_postprocess.R

* Extended Figure 3b: After running Tilted-CCA on the human PBMC CITE-seq data, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/pbmc_224antibody_tiltedcca_postprocess.R

* Extended Figure 3c: After running Tilted-CCA on the human PBMC 10x Multiome data, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/pbmc_10x_tiltedcca_postprocess.R

* Extended Figure 4a: After running Tilted-CCA on the human bone marrow Abseq data (461 genes), the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/abseq_bm_97antibodyRef_tcca_postprocess.R

* Extended Figure 4b: After running Tilted-CCA on the human bone marrow Abseq data (461 genes) and selecting the antibodies, the feature plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/abseq_bm_97antibodyRef_antibodySelection_postprocess.R

* Extended Figure 4c: After running Tilted-CCA on the human bone marrow Abseq data (461 genes), the differential-expression for antibodies is computed in https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/abseq_bm_97antibodyRef_differential.R. Then, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/abseq_bm_97antibodyRef_differential_postprocess.R

* Extended Figure 4d,e: After designing the alternative antibody panels, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/abseq_bm_97antibodyRef_antibodySelection_alternatives_postprocess.R

* Extended Figure 5a: After running Tilted-CCA on the mouse embryo development 10x Multiome data, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/reik_tiltedcca_postprocess.R

* Extended Figure 5b: After running Tilted-CCA (RNA+ATAC) on the mouse brain development 10x Multiome data, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/mouseembryo_tiltedcca_RNA-ATAC_postprocess.R

* Extended Figure 5c: After running Tilted-CCA (RNA+ChromatinActivity) on the mouse brain development 10x Multiome data, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/mouseembryo_tiltedcca_RNA-geneActivity_postprocess.R

* Extended Figure 5d: After running Tilted-CCA (RNA+ATAC) on the human brain development 10x Multiome data, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_tiltedcca_RNA-ATAC_postprocess.R

* Extended Figure 5e: After running Tilted-CCA (RNA+ChromatinActivity) on the human brain development 10x Multiome data, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_tiltedcca_RNA-geneActivity_postprocess.R

* Extended Figure 6a (Mouse brain development, 10x Multiome): After running Tilted-CCA (RNA+Chromatin Activity), the plot is created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/mouseembryo_steadystate_geneActivity.R

* Extended Figure 6a (Human brain development, 10x Multiome): After running Tilted-CCA (RNA+Chromatin Activity), the plot is created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_steadystate_RNA-geneActivity.R

* Extended Figure 6b (Mouse brain development, 10x Multiome): After selected the genes and ordering them to construct Figure 5j, the heatmap is created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/mouseembryo_developmentalGenes_heatmap-geneActivity.R

* Extended Figure 6b (Human brain development, 10x Multiome): After selected the genes and ordering them to construct Figure 5k, the heatmap is created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_developmentalGenes_heatmap-geneActivity.R

* Extended Figure 6c: After selected the genes to construct Figure 5k, the phase portraits are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/greenleaf_leafplots_additional.R

* Supplementary Figure 1: After running scAI on the human bone marrow CITE-seq data, the plot is created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_25antibody_scai_postprocess.R

* Supplementary Figure 2a,c (Human bone marrow CITE-seq): This plot is created in the same script that produces Figure 3c

* Supplementary Figure 2a,c (Human bone marrow Abseq, whole transcriptome): After running Tilted-CCA,  the differential-expression for genes and antibodies is computed in https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_97antibody_differential.R. Then, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_97antibody_differential_postprocess.R

* Supplementary Figure 2b (Human bone marrow Abseq, whole transcriptome): To run Tilted-CCA for only 25 antibodies for the protein modality, see https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_97antibody_tiltedcca-25antibodies.R. Afterwards, combine with the differential expression of each gene as calculated for Supplementary Figure 2a. Then, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/bm_97antibody_differential_postprocess.R

* Supplementary Figure 2b (Human PBMC CITE-seq): To run Tilted-CCA for only 25 antibodies for the protein modality, see https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/pbmc_224antibody_tiltedcca-25antibodies.R. Afterwards, combine with the differential expression of each gene as calculated for Figure 3e. Then, the plots are created via https://github.com/linnykos/tiltedCCA_analysis/blob/master/main/pbmc_224antibody_differential-25antibodies_postprocess.R

* Supplementary Figure 2c (Human PBMC CITE-seq): This plot is created in the same script that produces Figure 3e


# Setup

The following shows the suggested package versions that the developer (GitHub username: linnykos) used when developing the Tilted-CCA package.

```R
> devtools::session_info()
─ Session info ─────────────────────────────────────────────────────
 setting  value
 version  R version 4.1.2 (2021-11-01)
 os       Red Hat Enterprise Linux
 system   x86_64, linux-gnu
 ui       X11
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       America/New_York
 date     2022-10-04
 pandoc   1.12.3.1 @ /usr/bin/pandoc

─ Packages ─────────────────────────────────────────────────────────
 package              * version    date (UTC) lib source
 abind                  1.4-5      2016-07-21 [1] CRAN (R 4.1.2)
 ArchR                * 1.0.2      2022-07-14 [1] Github (GreenleafLab/ArchR@79953a9)
 assertthat             0.2.1      2019-03-21 [1] CRAN (R 4.1.2)
 basilisk               1.6.0      2021-10-26 [1] Bioconductor
 basilisk.utils         1.6.0      2021-10-26 [1] Bioconductor
 Biobase              * 2.54.0     2021-10-26 [1] Bioconductor
 BiocGenerics         * 0.40.0     2021-10-26 [1] Bioconductor
 BiocParallel           1.28.3     2021-12-09 [1] Bioconductor
 Biostrings             2.62.0     2021-10-26 [1] Bioconductor
 bitops                 1.0-7      2021-04-24 [1] CRAN (R 4.1.2)
 bmcite.SeuratData    * 0.3.0      2022-04-12 [1] local
 brio                   1.1.3      2021-11-30 [1] CRAN (R 4.1.2)
 cachem                 1.0.6      2021-08-19 [1] CRAN (R 4.1.2)
 callr                  3.7.1      2022-07-13 [1] CRAN (R 4.1.2)
 cli                    3.3.0      2022-04-25 [1] CRAN (R 4.1.2)
 cluster                2.1.2      2021-04-17 [2] CRAN (R 4.1.2)
 codetools              0.2-18     2020-11-04 [2] CRAN (R 4.1.2)
 colorspace             2.0-3      2022-02-21 [1] CRAN (R 4.1.2)
 corrplot               0.92       2021-11-18 [1] CRAN (R 4.1.2)
 cowplot                1.1.1      2020-12-30 [1] CRAN (R 4.1.2)
 crayon                 1.5.1      2022-03-26 [1] CRAN (R 4.1.2)
 data.table           * 1.14.2     2021-09-27 [1] CRAN (R 4.1.2)
 DBI                    1.1.3      2022-06-18 [1] CRAN (R 4.1.2)
 dbscan               * 1.1-10     2022-01-15 [1] CRAN (R 4.1.2)
 DelayedArray           0.20.0     2021-10-26 [1] Bioconductor
 deldir                 1.0-6      2021-10-23 [1] CRAN (R 4.1.2)
 devtools               2.4.4      2022-07-20 [1] CRAN (R 4.1.2)
 digest                 0.6.29     2021-12-01 [1] CRAN (R 4.1.2)
 dir.expiry             1.2.0      2021-10-26 [1] Bioconductor
 dplyr                  1.0.9      2022-04-28 [1] CRAN (R 4.1.2)
 ellipsis               0.3.2      2021-04-29 [1] CRAN (R 4.1.2)
 fansi                  1.0.3      2022-03-24 [1] CRAN (R 4.1.2)
 fastmap                1.1.0      2021-01-25 [1] CRAN (R 4.1.2)
 fastmatch              1.1-3      2021-07-23 [1] CRAN (R 4.1.2)
 filelock               1.0.2      2018-10-05 [1] CRAN (R 4.1.2)
 fitdistrplus           1.1-8      2022-03-10 [1] CRAN (R 4.1.2)
 forcats                0.5.1      2021-01-27 [1] CRAN (R 4.1.2)
 fs                     1.5.2      2021-12-08 [1] CRAN (R 4.1.2)
 future                 1.27.0     2022-07-22 [1] CRAN (R 4.1.2)
 future.apply           1.9.0      2022-04-25 [1] CRAN (R 4.1.2)
 generics               0.1.3      2022-07-05 [1] CRAN (R 4.1.2)
 GenomeInfoDb         * 1.30.1     2022-01-30 [1] Bioconductor
 GenomeInfoDbData       1.2.7      2022-03-09 [1] Bioconductor
 GenomicRanges        * 1.46.1     2021-11-18 [1] Bioconductor
 ggplot2              * 3.3.6      2022-05-03 [1] CRAN (R 4.1.2)
 ggrepel              * 0.9.1      2021-01-15 [1] CRAN (R 4.1.2)
 ggridges               0.5.3      2021-01-08 [1] CRAN (R 4.1.2)
 globals                0.16.0     2022-08-05 [1] CRAN (R 4.1.2)
 glue                   1.6.2      2022-02-24 [1] CRAN (R 4.1.2)
 goftest                1.2-3      2021-10-07 [1] CRAN (R 4.1.2)
 gridExtra            * 2.3        2017-09-09 [1] CRAN (R 4.1.2)
 gtable               * 0.3.0      2019-03-25 [1] CRAN (R 4.1.2)
 gtools               * 3.9.3      2022-07-11 [1] CRAN (R 4.1.2)
 HDF5Array              1.22.1     2021-11-14 [1] Bioconductor
 htmltools              0.5.3      2022-07-18 [1] CRAN (R 4.1.2)
 htmlwidgets            1.5.4      2021-09-08 [1] CRAN (R 4.1.2)
 httpuv                 1.6.5      2022-01-05 [1] CRAN (R 4.1.2)
 httr                   1.4.3      2022-05-04 [1] CRAN (R 4.1.2)
 ica                    1.0-3      2022-07-08 [1] CRAN (R 4.1.2)
 igraph               * 1.3.4      2022-07-19 [1] CRAN (R 4.1.2)
 IRanges              * 2.28.0     2021-10-26 [1] Bioconductor
 irlba                * 2.3.5      2021-12-06 [1] CRAN (R 4.1.2)
 jsonlite               1.8.0      2022-02-22 [1] CRAN (R 4.1.2)
 KernSmooth             2.23-20    2021-05-03 [2] CRAN (R 4.1.2)
 later                  1.3.0      2021-08-18 [1] CRAN (R 4.1.2)
 lattice                0.20-45    2021-09-22 [2] CRAN (R 4.1.2)
 lazyeval               0.2.2      2019-03-15 [1] CRAN (R 4.1.2)
 leiden                 0.4.2      2022-05-09 [1] CRAN (R 4.1.2)
 lifecycle              1.0.1      2021-09-24 [1] CRAN (R 4.1.2)
 listenv                0.8.0      2019-12-05 [1] CRAN (R 4.1.2)
 lmtest                 0.9-40     2022-03-21 [1] CRAN (R 4.1.2)
 magrittr             * 2.0.3      2022-03-30 [1] CRAN (R 4.1.2)
 MASS                 * 7.3-54     2021-05-03 [2] CRAN (R 4.1.2)
 MAST                 * 1.20.0     2021-10-26 [1] Bioconductor
 Matrix               * 1.3-4      2021-06-01 [2] CRAN (R 4.1.2)
 MatrixGenerics       * 1.6.0      2021-10-26 [1] Bioconductor
 matrixStats          * 0.62.0     2022-04-19 [1] CRAN (R 4.1.2)
 memoise                2.0.1      2021-11-26 [1] CRAN (R 4.1.2)
 mgcv                   1.8-38     2021-10-06 [2] CRAN (R 4.1.2)
 mime                   0.12       2021-09-28 [1] CRAN (R 4.1.2)
 miniUI                 0.1.1.1    2018-05-18 [1] CRAN (R 4.1.2)
 MOFA2                * 1.4.0      2021-10-26 [1] Bioconductor
 munsell                0.5.0      2018-06-12 [1] CRAN (R 4.1.2)
 nlme                   3.1-153    2021-09-07 [2] CRAN (R 4.1.2)
 parallelly             1.32.1     2022-07-21 [1] CRAN (R 4.1.2)
 patchwork              1.1.1      2020-12-17 [1] CRAN (R 4.1.2)
 pbapply                1.5-0      2021-09-16 [1] CRAN (R 4.1.2)
 pheatmap               1.0.12     2019-01-04 [1] CRAN (R 4.1.2)
 pillar                 1.8.0      2022-07-18 [1] CRAN (R 4.1.2)
 pkgbuild               1.3.1      2021-12-20 [1] CRAN (R 4.1.2)
 pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.1.2)
 pkgload                1.3.0      2022-06-27 [1] CRAN (R 4.1.2)
 plotly                 4.10.0     2021-10-09 [1] CRAN (R 4.1.2)
 plyr                 * 1.8.7      2022-03-24 [1] CRAN (R 4.1.2)
 png                    0.1-7      2013-12-03 [1] CRAN (R 4.1.2)
 polyclip               1.10-0     2019-03-14 [1] CRAN (R 4.1.2)
 prettyunits            1.1.1      2020-01-24 [1] CRAN (R 4.1.2)
 princurve            * 2.1.6      2021-01-18 [1] CRAN (R 4.1.2)
 processx               3.7.0      2022-07-07 [1] CRAN (R 4.1.2)
 profvis                0.3.7      2020-11-02 [1] CRAN (R 4.1.2)
 progressr              0.10.1     2022-06-03 [1] CRAN (R 4.1.2)
 promises               1.2.0.1    2021-02-11 [1] CRAN (R 4.1.2)
 ps                     1.7.0      2022-04-23 [1] CRAN (R 4.1.2)
 purrr                  0.3.4      2020-04-17 [1] CRAN (R 4.1.2)
 quadprog             * 1.5-8      2019-11-20 [1] CRAN (R 4.1.2)
 R6                     2.5.1      2021-08-19 [1] CRAN (R 4.1.2)
 RANN                 * 2.6.1      2019-01-08 [1] CRAN (R 4.1.2)
 rappdirs               0.3.3      2021-01-31 [1] CRAN (R 4.1.2)
 RColorBrewer         * 1.1-3      2022-04-03 [1] CRAN (R 4.1.2)
 Rcpp                 * 1.0.9      2022-07-08 [1] CRAN (R 4.1.2)
 RcppAnnoy              0.0.19     2021-07-30 [1] CRAN (R 4.1.2)
 RcppRoll               0.3.0      2018-06-05 [1] CRAN (R 4.1.2)
 RCurl                  1.98-1.8   2022-07-30 [1] CRAN (R 4.1.2)
 remotes                2.4.2      2021-11-30 [1] CRAN (R 4.1.2)
 reshape2               1.4.4      2020-04-09 [1] CRAN (R 4.1.2)
 reticulate             1.25       2022-05-11 [1] CRAN (R 4.1.2)
 rgeos                  0.5-9      2021-12-15 [1] CRAN (R 4.1.2)
 rhdf5                * 2.38.1     2022-03-10 [1] Bioconductor
 rhdf5filters           1.6.0      2021-10-26 [1] Bioconductor
 Rhdf5lib               1.16.0     2021-10-26 [1] Bioconductor
 rlang                  1.0.4      2022-07-12 [1] CRAN (R 4.1.2)
 ROCR                   1.0-11     2020-05-02 [1] CRAN (R 4.1.2)
 rpart                  4.1-15     2019-04-12 [2] CRAN (R 4.1.2)
 Rsamtools              2.10.0     2021-10-26 [1] Bioconductor
 RSpectra             * 0.16-1     2022-04-24 [1] CRAN (R 4.1.2)
 Rtsne                  0.16       2022-04-17 [1] CRAN (R 4.1.2)
 S4Vectors            * 0.32.4     2022-03-24 [1] Bioconductor
 scales               * 1.2.0      2022-04-13 [1] CRAN (R 4.1.2)
 scattermore            0.8        2022-02-14 [1] CRAN (R 4.1.2)
 sctransform            0.3.3      2022-01-13 [1] CRAN (R 4.1.2)
 sessioninfo            1.2.2      2021-12-06 [1] CRAN (R 4.1.2)
 Seurat               * 4.1.1      2022-05-02 [1] CRAN (R 4.1.2)
 SeuratData           * 0.2.1      2022-04-12 [1] Github (satijalab/seurat-data@cf04099)
 SeuratObject         * 4.1.0      2022-05-01 [1] CRAN (R 4.1.2)
 shiny                  1.7.2      2022-07-19 [1] CRAN (R 4.1.2)
 Signac               * 1.7.0      2022-06-01 [1] CRAN (R 4.1.2)
 SingleCellExperiment * 1.16.0     2021-10-26 [1] Bioconductor
 slingshot            * 2.3.1      2022-04-25 [1] Github (kstreet13/slingshot@cf69e1a)
 sp                   * 1.5-0      2022-06-05 [1] CRAN (R 4.1.2)
 sparseMatrixStats    * 1.6.0      2021-10-26 [1] Bioconductor
 spatstat.core          2.4-2      2022-04-01 [1] CRAN (R 4.1.2)
 spatstat.data          2.2-0      2022-04-18 [1] CRAN (R 4.1.2)
 spatstat.geom          2.4-0      2022-03-29 [1] CRAN (R 4.1.2)
 spatstat.random        2.2-0      2022-03-30 [1] CRAN (R 4.1.2)
 spatstat.sparse        2.1-1      2022-04-18 [1] CRAN (R 4.1.2)
 spatstat.utils         2.3-1      2022-05-06 [1] CRAN (R 4.1.2)
 stringi                1.7.8      2022-07-11 [1] CRAN (R 4.1.2)
 stringr              * 1.4.0      2019-02-10 [1] CRAN (R 4.1.2)
 SummarizedExperiment * 1.24.0     2021-10-26 [1] Bioconductor
 survival               3.2-13     2021-08-24 [2] CRAN (R 4.1.2)
 tensor                 1.5        2012-05-05 [1] CRAN (R 4.1.2)
 testthat             * 3.1.4      2022-04-26 [1] CRAN (R 4.1.2)
 tibble                 3.1.8      2022-07-22 [1] CRAN (R 4.1.2)
 tidyr                  1.2.0      2022-02-01 [1] CRAN (R 4.1.2)
 tidyselect             1.1.2      2022-02-21 [1] CRAN (R 4.1.2)
 tiltedCCA            * 1.0.0.001  2022-09-06 [1] local
 TrajectoryUtils      * 1.2.0      2021-10-26 [1] Bioconductor
 urlchecker             1.0.1      2021-11-30 [1] CRAN (R 4.1.2)
 usethis                2.1.6      2022-05-25 [1] CRAN (R 4.1.2)
 utf8                   1.2.2      2021-07-24 [1] CRAN (R 4.1.2)
 uwot                   0.1.13     2022-08-16 [1] CRAN (R 4.1.2)
 vctrs                  0.4.1      2022-04-13 [1] CRAN (R 4.1.2)
 viridisLite            0.4.0      2021-04-13 [1] CRAN (R 4.1.2)
 withr                  2.5.0      2022-03-03 [1] CRAN (R 4.1.2)
 xtable                 1.8-4      2019-04-21 [1] CRAN (R 4.1.2)
 XVector                0.34.0     2021-10-26 [1] Bioconductor
 zlibbioc               1.40.0     2021-10-26 [1] Bioconductor
 zoo                    1.8-10     2022-04-15 [1] CRAN (R 4.1.2)
```

# References
[1] Stuart, T. et al. Comprehensive integration of single-cell data. Cell 177, 1888–1902 (2019).

[2] Triana, S. et al. Single-cell proteo-genomic reference maps of the hematopoietic system enable the purification and massive profiling of precisely defined cell states. Nature immunology 22, 1577–1589 (2021).

[3] Hao, Y. et al. Integrated analysis of multimodal single-cell data. Cell 184, 3573–3587 (2021).

[4] Trevino, A. E. et al. Chromatin and gene-regulatory dynamics of the developing human cerebral cortex at single-cell resolution. Cell (2021).

[5] Argelaguet, R. et al. Decoding gene regulation in the mouse embryo using single-cell multi-omics. bioRxiv (2022).

[6] Hounkpe, B. W., Chenou, F., de Lima, F. & De Paula, E. V. HRT atlas v1. 0 database: Redefining human and mouse housekeeping genes and candidate reference transcripts by mining massive rna-seq datasets. Nucleic Acids Research 49, D947–D955 (2021).

[7] Riba, A. et al. Cell cycle gene regulation dynamics revealed by RNA velocity and deep-learning. Nature Communications 13, 1–13 (2022).