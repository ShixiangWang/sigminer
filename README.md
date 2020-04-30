
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sigminer: an easy-to-use and scalable toolkit for genomic alteration signature analysis and visualization in R <img src="man/figures/logo.png" align="right"  height="140" width="120" alt="logo"/>

[![CRAN
status](https://www.r-pkg.org/badges/version/sigminer)](https://cran.r-project.org/package=sigminer)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-blue.svg)](https://www.tidyverse.org/lifecycle/#stable)
![R-CMD-check](https://github.com/ShixiangWang/sigminer/workflows/R-CMD-check/badge.svg)
[![Coverage
status](https://codecov.io/gh/ShixiangWang/sigminer/branch/master/graph/badge.svg)](https://codecov.io/github/ShixiangWang/sigminer?branch=master)
[![](http://cranlogs.r-pkg.org/badges/grand-total/sigminer?color=orange)](https://cran.r-project.org/package=sigminer)
[![Closed
issues](https://img.shields.io/github/issues-closed/ShixiangWang/sigminer.svg)](https://github.com/ShixiangWang/sigminer/issues?q=is%3Aissue+is%3Aclosed)

## Overview

Genomic alterations including single nucleotide substitution (SBS), copy
number alteration (CNA), etc. are the major force for cancer
initialization and development. Due to the specificity of molecular
lesions caused by genomic alterations, we can generate characteristic
alteration spectra, called ‘signature’. This package helps users to
extract, analyze and visualize signatures from genomic alteration
records, thus providing new insight into cancer study.

**SBS signatures**:

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />

**Copy number signatures**:

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

### Feature

  - supports a standard *de novo* pipeline for identification of two
    types of signatures: SBS signature and copy number signature
  - supports quantify exposure for one sample based on *known
    signatures*
  - supports two methods for calling copy number signatures: one is from
    Macintyre et al. 2018 and the other is created by us
  - supports association and group analysis and visualization for
    signatures
  - supports a bayesian variant of NMF algorithm to enable optimal
    inferences for the number of signatures through the automatic
    relevance determination technique from **SignatureAnalyzer** package
  - supports two plot styles for signature profile: ‘default’ (like
    **SignatureAnalyzer** package) and ‘cosmic’ (like **COSMIC**
    database)
  - supports two types of signatrue exposures: relative exposure
    (relative contribution of signatures in each sample) and absolute
    exposure (estimated variation records of signatures in each sample)
  - supports basic summary and visualization for profile of mutation
    (powered by **maftools**) and copy number
  - supports parallel computation by R packages **foreach**, **future**
    and **NMF**
  - efficient code powered by R packages **data.table** and
    **tidyverse**
  - elegant plots powered by R packages **ggplot2**, **ggpubr**,
    **cowplot** and **patchwork**
  - well tested by R package **testthat** and documented by R package
    **roxygen2**, **roxytest**, **pkgdown**, and etc. for both reliable
    and reproducible research

## Installation

You can install the stable release of **sigminer** from CRAN with:

``` r
install.packages("sigminer", dependencies = TRUE)
# Or
BiocManager::install("sigminer", dependencies = TRUE)
```

You can install the development version of **sigminer** from Github
with:

``` r
remotes::install_github("ShixiangWang/sigminer", dependencies = TRUE)
# For Chinese users, run 
remotes::install_git("https://gitee.com/ShixiangWang/sigminer", dependencies = TRUE)
```

## Usage

A complete documentation of **sigminer** can be read online at
<https://shixiangwang.github.io/sigminer-doc/> (For Chinese users, you
can also read it at <https://shixiangwang.gitee.io/sigminer-doc>). All
functions are well organized and documented at
<https://shixiangwang.github.io/sigminer/reference/index.html> (For
Chinese users, you can also read it at
<https://shixiangwang.gitee.io/sigminer/reference/index.html>). For
usage of a specific function `fun`, run `?fun` in your R console to see
its documentation.

## Citation

-----

***Wang, Shixiang, et al. “Copy number signature analyses in prostate
cancer reveal distinct etiologies and clinical outcomes” medRxiv
(2020)***
<https://www.medrxiv.org/content/early/2020/04/29/2020.04.27.20082404>

-----

## Acknowledgments

If you use **NMF** package in R, please also cite:

    Gaujoux, Renaud, and Cathal Seoighe. "A Flexible R Package for 
        Nonnegative Matrix Factorization."" BMC Bioinformatics 11, no. 1 (December 2010).

The method “M” for extracting copy number signatures was based in part
on the source code from paper *Copy number signatures and mutational
processes in ovarian carcinoma*, if you use this feature, please also
cite:

    Macintyre, Geoff, et al. "Copy number signatures and mutational
        processes in ovarian carcinoma." Nature genetics 50.9 (2018): 1262.

The code for extracting SBS signatures was based in part on the source
code of the **maftools** package, if you use this feature, please also
cite:

    Mayakonda, Anand, et al. "Maftools: efficient and comprehensive analysis
        of somatic variants in cancer." Genome research 28.11 (2018): 1747-1756.

The code for extracting mutational signatures was based in part on the
source code of the **SignatureAnalyzer** package, if you use this
feature, please also
    cite:

    Kim, Jaegil, et al. "Somatic ERCC2 mutations are associated with a distinct genomic
        signature in urothelial tumors." Nature genetics 48.6 (2016): 600.

## References

1.  Alexandrov, Ludmil B., et al. “The repertoire of mutational
    signatures in human cancer.” Nature 578.7793 (2020): 94-101.
2.  Macintyre, Geoff, et al. “Copy number signatures and mutational
    processes in ovarian carcinoma.” Nature genetics 50.9 (2018): 1262.
3.  Mayakonda, Anand, et al. “Maftools: efficient and comprehensive
    analysis of somatic variants in cancer.” Genome research 28.11
    (2018): 1747-1756.
4.  Gaujoux, Renaud, and Cathal Seoighe. “A Flexible R Package for
    Nonnegative Matrix Factorization.”" BMC Bioinformatics 11, no. 1
    (December 2010).
5.  H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
    Springer-Verlag New York, 2016.
6.  Tan, Vincent YF, and Cédric Févotte. “Automatic relevance
    determination in nonnegative matrix factorization with the/spl
    beta/-divergence.” IEEE Transactions on Pattern Analysis and Machine
    Intelligence 35.7 (2012): 1592-1605.
7.  Kim, Jaegil, et al. “Somatic ERCC2 mutations are associated with a
    distinct genomic signature in urothelial tumors.” Nature genetics
    48.6 (2016): 600.

## LICENSE

The software is made available for non commercial research purposes only
under the
[MIT](https://github.com/ShixiangWang/sigminer/blob/master/LICENSE.md).
However, notwithstanding any provision of the MIT License, the software
currently may not be used for commercial purposes without explicit
written permission after contacting Shixiang Wang
<wangshx@shanghaitech.edu.cn> or Xue-Song Liu
<liuxs@shanghaitech.edu.cn>.

MIT © 2019-2020 Shixiang Wang, Xue-Song Liu

MIT © 2018 Geoffrey Macintyre

MIT © 2018 Anand Mayakonda

-----

**[Cancer Biology Group](https://github.com/XSLiuLab) @ShanghaiTech**

**Research group led by Xue-Song Liu in ShanghaiTech University**
