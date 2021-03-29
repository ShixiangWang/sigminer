# sigminer: 癌症基因组突变模式分析工具箱 <img src="man/figures/logo.png" align="right"  height="140" width="120" alt="logo"/>

[![CRAN
status](https://www.r-pkg.org/badges/version/sigminer)](https://cran.r-project.org/package=sigminer)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-blue.svg)](https://www.tidyverse.org/lifecycle/#stable)
![R-CMD-check](https://github.com/ShixiangWang/sigminer/workflows/R-CMD-check/badge.svg)
[![Coverage
status](https://codecov.io/gh/ShixiangWang/sigminer/branch/master/graph/badge.svg)](https://codecov.io/github/ShixiangWang/sigminer?branch=master)
[![](http://cranlogs.r-pkg.org/badges/grand-total/sigminer?color=orange)](https://cran.r-project.org/package=sigminer)
[![Closed
issues](https://img.shields.io/github/issues-closed/ShixiangWang/sigminer.svg)](https://github.com/ShixiangWang/sigminer/issues?q=is%3Aissue+is%3Aclosed)


## 安装

从 CRAN 上安装稳定版本：

```r
install.packages("BiocManager")
BiocManager::install("sigminer", dependencies = TRUE)
```

从 Gitee 上安装开发版本：

```r
remotes::install_git("https://gitee.com/ShixiangWang/sigminer", dependencies = TRUE)
```

使用 conda：

```sh
# 注意：bioconda 版本可能不是最新版本

# 你可以先安装一个独立的 conda 环境
# conda create -n sigminer
# conda activate sigminer
conda install -c bioconda -c conda-forge r-sigminer
```

## 使用文档

完整的文档请阅读 <https://shixiangwang.gitee.io/sigminer-doc>，
整理好的开放的函数列表请阅读 <https://shixiangwang.gitee.io/sigminer/reference/index.html>。

## 引用

-----

***

***Wang, Shixiang, et al. "Copy number signature analyses in prostate cancer reveal distinct etiologies and clinical outcomes" medRxiv (2020)*** <https://www.medrxiv.org/content/early/2020/04/29/2020.04.27.20082404>

***

-----

Currently, **sigminer** has not been published. If you use **sigminer**
in academic field, please cite our recent
    work:

    Wang, Shixiang, et al. "The predictive power of tumor mutational burden 
        in lung cancer immunotherapy response is influenced by patients' sex."
        International journal of cancer (2019).

and the paper implementing NMF in R:

    Gaujoux, Renaud, and Cathal Seoighe. "A Flexible R Package for 
        Nonnegative Matrix Factorization."" BMC Bioinformatics 11, no. 1 (December 2010).
