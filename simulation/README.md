# Purpose

The simulations here are meant to be portable, toy-examples of Tilted-CCA, with the intent to demonstrate that an installation of Tilted-CCA was successful as well as to demonstrate how it differs from the methods that seek to find the "union of information".

A successful installation of `tiltedCCA` is required for this. See the 
last section of this README of all the system/package dependencies used when creating these simulations.

# Overview

## Simulation 1: Modality 1 separates cells into 3 clusters, and Modality 2 does not

## Simulation 2: Both modalities separate the 5 cell-types into 3 clusters in different ways

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