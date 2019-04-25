
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sigminer <img src="man/figures/logo.png" align="right"  height="140" width="120" alt="logo"/>

[![CRAN
status](https://www.r-pkg.org/badges/version/sigminer)](https://cran.r-project.org/package=sigminer)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/ShixiangWang/sigminer?branch=master&svg=true)](https://ci.appveyor.com/project/ShixiangWang/sigminer)
[![Travis build
status](https://travis-ci.org/ShixiangWang/sigminer.svg?branch=master)](https://travis-ci.org/ShixiangWang/sigminer)
[![Coverage
status](https://codecov.io/gh/ShixiangWang/sigminer/branch/master/graph/badge.svg)](https://codecov.io/github/ShixiangWang/sigminer?branch=master)
[![emoji-log](https://cdn.rawgit.com/ahmadawais/stuff/ca97874/emoji-log/non-flat-round.svg)](https://github.com/ahmadawais/Emoji-Log/)

## Overview

The goal of **sigminer** is to provide an uniform interface for genomic
variation signature analysis and visualization.

**sigminer** is powered by [NMF](https://github.com/renozao/NMF) package
and [maftools](https://github.com/PoisonAlien/maftools) package.

## Installation

You can install the development version of sigminer from Github with:

``` r
remotes::install_github("ShixiangWang/sigminer")
```

## Citation

If you use **sigminer** in academic field, please
    cite:

    Wang, Shixiang, et al. "The predictive power of tumor mutational burden 
        in lung cancer immunotherapy response is influenced by patients' sex."
        International journal of cancer (2019).

and

    Gaujoux, Renaud, and Cathal Seoighe. "A Flexible R Package for 
        Nonnegative Matrix Factorization."" BMC Bioinformatics 11, no. 1 (December 2010).

## Acknowledgments

The code for extracting copy number signatures was based in part on the
source code from paper *Copy number signatures and mutational processes
in ovarian carcinoma*, if you use this feature, please also cite:

    Macintyre, Geoff, et al. "Copy number signatures and mutational
        processes in ovarian carcinoma." Nature genetics 50.9 (2018): 1262.

The code for extracting mutational signatures was based in part on the
source code of the **maftools** package, if you use this feature, please
also
    cite:

    Mayakonda, Anand, et al. "Maftools: efficient and comprehensive analysis
        of somatic variants in cancer." Genome research 28.11 (2018): 1747-1756.

## LICENSE

MIT © 2019 Shixiang Wang, Geoffrey Macintyre, Xue-Song Liu
