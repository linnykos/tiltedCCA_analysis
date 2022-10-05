# Purpose

The simulations here are meant to be portable, toy-examples of Tilted-CCA, with the intent to demonstrate that an installation of Tilted-CCA was successful as well as to demonstrate how it differs from the methods that seek to find the "union of information".

A successful installation of `tiltedCCA` is required for this. See the 
last section of this README of all the system/package dependencies used when creating these simulations. Both simulations should complete in less than 2 minutes.

# Overview

Broadly speaking, Tilted-CCA is a pipeline of 6 different function calls. The procedure starts by assuming there are two data matrices, `mat_1` and `mat_2`, (either `matrix` or `dgCMatrix`) with the same number of rows (i.e., the same cells), but possibly different features. We advise making sure each row/column of either matrix has non-zero variance prior to using this pipeline.

* **Step 1 (Initializing low-rank representation of each modality, `create_multiSVD`):** We compute the leading PCs of each modality here. There can be a different number of latent dimensions for each modality. Here, the additional parameters dictate how to compute the leading PCs (for example, should the features be scaled/center before or after the computing the PCs? Which PCs should be used downstream? This would be useful since for ATAC data, it is common to not use the first leading PC.)

The output of this step is a `multiSVD` object, which will continually be updated as we feed it as input and get an updated object after each of the following function calls. It currently has elements `svd_1`, `svd_2`, `default_assay`, and `param`. 

* **Step 2 (Passing helping information, `form_metacells`):** We pass (optional) large clustering structure or meta-cells here. The former is designed to help aid Tilted-CCA assess what the "intersection of information" is, and the latter is designed for handling large datasets with more than 10,000 cells. If no large structuring structure is needed, the practitioner can simply pass in `large_clustering_1=NULL` and `large_clustering_2=NULL`. If no meta-cells are needed, the practitioner can simply pass in `num_metacells=NULL`.

The updated `multiSVD` object will now have an additional element: `metacell_obj`.

* **Step 3 (Creating shared nearest neighbor (SNN) graphs, `compute_snns`):** We compute the SNN graphs for both modalities and the target common manifold (computed based on both modality's SNN graph). These graphs dictate how Tilted-CCA represents the information in each modality, and each node of the graph represents a cell/meta-cell. Here, the additional parameters dictate the specific details on how these graphs are constructed. The two important parameters `latent_k` (i.e., the dimensionality of the graph Laplacian eigenbases, which is used when optimizing the tilt) and `num_neigh` (i.e., the number of neighbors for each cell in either modality).

The updated `multiSVD` object will now have additional elements: `snn_list` and `laplacian_list`.

* **Step 4 (Initializing Tilted-CCA, `tiltedCCA`):** We compute the CCA between both modalities and initialize the tilt across all latent dimensions to a particular value. Recall that CCA's solution can be directly computed using the PC's from each modality. If there are different number of latent dimensions for each modality, Tilted-CCA will have a latent dimensionality (i.e., via CCA) equal to the smaller of the two.

The updated `multiSVD` object will now have additional elements: `cca_obj` and `tcca_obj`.

* **Step 5 (Tuning the tilt across each latent dimension, `fine_tuning`):** We tune the tilt across each latent dimension. This is the main computational bottleneck of Tilted-CCA, as an optimization is performed cyclically (i.e., the tilt for each latent dimension is updated in sequence, with a few epochs cycling through all latent dimensions). 

There is no new obvious element added to `multiSVD`, but the representation of the common and distinct embeddings is updated upon this function's completion. For example, `multiSVD_obj$tcca_obj$tilt_perc` now is a vector, denoting the tilt of the common vector for each latent dimension.

* **Step 6 (Completing the decomposition, `tiltedCCA_decomposition`):** Given the tilts of common vectors across all latent dimensions, we now recover the full cell-by-feature decomposition of each modality. Here, the additional parameters dictate whether you want the decompositions to be cell-by-feature or cell-by-PC. The latter is desirable especially if you are analyzing ATAC data, where the cell-by-feature would be a dense matrix with many hundred-of-thousand features (which would be too memory intensive).

The updated `multiSVD` object will now have additional elements: `common_mat_1` and `distint_mat_1` (or `common_dimred_1` and `distinct_dimred_1` if `bool_modality_1_full=F`)
and `common_mat_2` and `distint_mat_2` (or `common_dimred_2` and `distinct_dimred_2` if `bool_modality_2_full=F`)



## Simulation 1: Modality 1 separates cells into 3 clusters, and Modality 2 does not

See https://github.com/linnykos/tiltedCCA_analysis/blob/master/simulation/simulation_1.R for this simulation. In this simulation, there are 300 cells across two modalities of 10 features each. Modality 1 has a "high" amount of distinct information -- it separates the 300 cells into 3 obvious clusters. Modality 2 has a "low" amount of distinct information -- all 300 cells are in one amorphous ball.

The following shows each modality based on their leading 2 PCs respectively, where the cells are colored by the true cell-types.

![simulation1_data.png](fig/simulation1_data.png | width=200)

The following shows Tilted-CCA's common (i.e., shared axes of variation between both modalities) and distinct (i.e., axes of variation unique to a modality, after the common axes have been accounted for) axes of variation, where the cells are colored by the true cell-types. Here, observe that the common embedding does not really contain information to separate the cell-types -- this is desirable, as all the cell-type separation information are unique to Modality 1 (hence, appearing in the Modality 1's distinct embedding). The common embedding demonstrates the "intersection of information."

![simulation1_tcca.png](fig/simulation1_tcca.png | width=300)

In contrast with Tilted-CCA's "intersection of information," we demonstrate Consensus PCA, which is a prototypical method illustrating the "union of information" where a low-dimensional embedding is constructed for multimodal data by combining the leading axes of variation from each modality.
These embeddings are useful to visualize the "best of both worlds" (i.e., having cell-type separation that best combines both modalities), but is not useful for to understand the shared and distinct signals between the two modalities. Here, observe that the three cell-types are clearly separated by Consensus PCA.

![simulation1_consensuspca.png](fig/simulation1_consensuspca.png | width=100)

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

─ Packages ────────────────────────
 brio                1.1.3     2021-11-30 [1] CRAN (R 4.0.2)
```