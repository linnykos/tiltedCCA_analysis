# Purpose

The simulations here are meant to be portable, toy-examples of Tilted-CCA, with the intent to demonstrate that an installation of Tilted-CCA was successful as well as to demonstrate how it differs from the methods that seek to find the "union of information".

A successful installation of `tiltedCCA` is required for this. See the 
last section of this README of all the system/package dependencies used when creating these simulations. Both simulations should complete in less than 2 minutes.

# Overview

Broadly speaking, Tilted-CCA is a pipeline of 6 different function calls. The procedure starts by assuming there are two data matrices, `mat_1` and `mat_2`, (either `matrix` or `dgCMatrix`) with the same number of rows (i.e., the same cells), but possibly different features. We advise making sure each row/column of either matrix has non-zero variance prior to using this pipeline.

* **Step 1 (Initializing low-rank representation of each modality, `create_multiSVD`):** We compute the leading PCs of each modality here. There can be a different number of latent dimensions for each modality. Here, the additional parameters dictate how to compute the leading PCs (for example, should the features be scaled/center before or after the computing the PCs? Which PCs should be used downstream? This would be useful since for ATAC data, it is common to not use the first leading PC.)

The output of this step is a `multiSVD` object, which will continually be updated as we feed it as input and get an updated object after each of the following function calls. It currently has elements `svd_1`, `svd_2`, `default_assay`, and `param`. 

* Step 2 (Passing helping information, `form_metacells`): We pass (optional) large clustering structure or meta-cells here. The former is designed to help aid Tilted-CCA assess what the "intersection of information" is, and the latter is designed for handling large datasets with more than 10,000 cells. If no large structuring structure is needed, the practitioner can simply pass in `large_clustering_1=NULL` and `large_clustering_2=NULL`. If no meta-cells are needed, the practitioner can simply pass in `num_metacells=NULL`.

The updated `multiSVD` object will now have an additional element: `metacell_obj`.

* Step 3 (Creating shared nearest neighbor (SNN) graphs, `compute_snns`): We compute the SNN graphs for both modalities and the target common manifold (computed based on both modality's SNN graph). These graphs dictate how Tilted-CCA represents the information in each modality, and each node of the graph represents a cell/meta-cell. Here, the additional parameters dictate the specific details on how these graphs are constructed. The two important parameters `latent_k` (i.e., the dimensionality of the graph Laplacian eigenbases, which is used when optimizing the tilt) and `num_neigh` (i.e., the number of neighbors for each cell in either modality).

The updated `multiSVD` object will now have additional elements: `snn_list` and `laplacian_list`.

* Step 4 (Initializing Tilted-CCA, `tiltedCCA`): We compute the CCA between both modalities and initialize the tilt across all latent dimensions to a particular value. Recall that CCA's solution can be directly computed using the PC's from each modality. If there are different number of latent dimensions for each modality, Tilted-CCA will have a latent dimensionality (i.e., via CCA) equal to the smaller of the two.

The updated `multiSVD` object will now have additional elements: `cca_obj` and `tcca_obj`.

* Step 5 (Tuning the tilt across each latent dimension, `fine_tuning`): We tune the tilt across each latent dimension. This is the main computational bottleneck of Tilted-CCA, as an optimization is performed cyclically (i.e., the tilt for each latent dimension is updated in sequence, with a few epochs cycling through all latent dimensions). 

There is no new obvious element added to `multiSVD`, but the representation of the common and distinct embeddings is updated upon this function's completion. For example, `multiSVD_obj$tcca_obj$tilt_perc` now is a vector, denoting the tilt of the common vector for each latent dimension.

* Step 6 (Completing the decomposition, `tiltedCCA_decomposition`): Given the tilts of common vectors across all latent dimensions, we now recover the full cell-by-feature decomposition of each modality. Here, the additional parameters dictate whether you want the decompositions to be cell-by-feature or cell-by-PC. The latter is desirable especially if you are analyzing ATAC data, where the cell-by-feature would be a dense matrix with many hundred-of-thousand features (which would be too memory intensive).

The updated `multiSVD` object will now have additional elements: `common_mat_1` and `distint_mat_1` (or `common_dimred_1` and `distinct_dimred_1` if `bool_modality_1_full=F`)
and `common_mat_2` and `distint_mat_2` (or `common_dimred_2` and `distinct_dimred_2` if `bool_modality_2_full=F`)



## Simulation 1: Modality 1 separates cells into 3 clusters, and Modality 2 does not

See https://github.com/linnykos/tiltedCCA_analysis/blob/master/simulation/simulation_1.R for this simulation.

|                                                     |
|-----------------------------------------------------|
| ![simulation1_data.png](fig/simulation1_data.png) |

|                                                     |
|-----------------------------------------------------|
| ![simulation1_tcca.png](fig/simulation1_tcca.png) |

|                                                     |
|-----------------------------------------------------|
| ![simulation1_consensuspca.png](fig/simulation1_consensuspca.png) |

## Simulation 2: Both modalities separate the 5 cell-types into 3 clusters in different ways

See https://github.com/linnykos/tiltedCCA_analysis/blob/master/simulation/simulation_2.R for this simulation.

|                                                     |
|-----------------------------------------------------|
| ![simulation2_data.png](fig/simulation2_data.png) |

|                                                     |
|-----------------------------------------------------|
| ![simulation2_tcca.png](fig/simulation2_tcca.png) |

|                                                     |
|-----------------------------------------------------|
| ![simulation2_consensuspca.png](fig/simulation2_consensuspca.png) |

# Setup

The following shows the suggested package versions that the developer (GitHub username: linnykos) used when developing the Tilted-CCA package.

```
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

─ Packages ──────────────────────── package           * version   date (UTC) lib source
 brio                1.1.3     2021-11-30 [1] CRAN (R 4.0.2) cachem              1.0.6     2021-08-19 [1] CRAN (R 4.0.2) callr               3.7.0     2021-04-20 [1] CRAN (R 4.0.2) cli                 3.2.0     2022-02-14 [1] CRAN (R 4.0.5) crayon              1.5.1     2022-03-26 [1] CRAN (R 4.0.5) desc                1.4.2     2022-09-08 [1] CRAN (R 4.0.5) devtools            2.4.3     2021-11-30 [1] CRAN (R 4.0.2) ellipsis            0.3.2     2021-04-29 [1] CRAN (R 4.0.2) fastmap             1.1.0     2021-01-25 [1] CRAN (R 4.0.2) fs                  1.5.2     2021-12-08 [1] CRAN (R 4.0.2) glue                1.6.2     2022-02-24 [1] CRAN (R 4.0.5) irlba               2.3.5     2021-12-06 [1] CRAN (R 4.0.2) lattice             0.20-45   2021-09-22 [1] CRAN (R 4.0.2) lifecycle           1.0.1     2021-09-24 [1] CRAN (R 4.0.2) magrittr            2.0.3     2022-03-30 [1] CRAN (R 4.0.5) MASS                7.3-56    2022-03-23 [1] CRAN (R 4.0.5) Matrix              1.4-1     2022-03-23 [1] CRAN (R 4.0.5) MatrixGenerics      1.2.1     2021-01-30 [1] Bioconductor matrixStats         0.61.0    2021-09-17 [1] CRAN (R 4.0.2) memoise             2.0.1     2021-11-26 [1] CRAN (R 4.0.2) pkgbuild            1.3.1     2021-12-20 [1] CRAN (R 4.0.2) pkgload             1.2.4     2021-11-30 [1] CRAN (R 4.0.2) prettyunits         1.1.1     2020-01-24 [1] CRAN (R 4.0.2) processx            3.5.3     2022-03-25 [1] CRAN (R 4.0.5) ps                  1.6.0     2021-02-28 [1] CRAN (R 4.0.2) purrr               0.3.4     2020-04-17 [1] CRAN (R 4.0.2) quadprog            1.5-8     2019-11-20 [1] CRAN (R 4.0.2) R6                  2.5.1     2021-08-19 [1] CRAN (R 4.0.2) RANN                2.6.1     2019-01-08 [1] CRAN (R 4.0.2) Rcpp                1.0.8.3   2022-03-17 [1] CRAN (R 4.0.5) remotes             2.4.2     2021-11-30 [1] CRAN (R 4.0.2) rlang               1.0.2     2022-03-04 [1] CRAN (R 4.0.5) rprojroot           2.0.3     2022-04-02 [1] CRAN (R 4.0.5) RSpectra            0.16-0    2019-12-01 [1] CRAN (R 4.0.2) sessioninfo         1.2.2     2021-12-06 [1] CRAN (R 4.0.5) sparseMatrixStats   1.2.1     2021-02-02 [1] Bioconductor testthat            3.1.3     2022-03-29 [1] CRAN (R 4.0.5) tiltedCCA         * 1.0.0.001 2022-10-05 [1] local usethis             2.1.6     2022-05-25 [1] CRAN (R 4.0.5) withr               2.5.0     2022-03-03 [1] CRAN (R 4.0.5)
```