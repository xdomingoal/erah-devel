# erah

<!-- badges: start -->
[![R-CMD-check](https://github.com/xdomingoal/erah-devel/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/xdomingoal/erah-devel/actions/workflows/R-CMD-check.yaml)
[![Coverage status](https://codecov.io/gh/xdomingoal/erah-devel/branch/master/graph/badge.svg)](https://codecov.io/github/xdomingoal/erah-devel?branch=devel)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/erah)](https://cran.r-project.org/package=erah)
[![](http://cranlogs.r-pkg.org/badges/erah)](http://cran.rstudio.com/web/packages/erah/index.html)
<!-- badges: end -->

### Automated Spectral Deconvolution, Alignment, and Metabolite Identification in GC/MS-Based Untargeted Metabolomics

Automated compound deconvolution, alignment across samples, and identification of metabolites by spectral library matching in Gas Chromatography - Mass spectrometry (GC-MS) untargeted metabolomics. Outputs a table with compound names, matching scores and the integrated area of the compound for each sample.

eRah development version. For downloads, please use CRAN
https://CRAN.R-project.org/package=erah 

Or from R console, execute:

```r
install.packages('erah')
```
To install this development version, execute:

```r
devtools::install_github('xdomingoal/erah-devel',build_vignettes = TRUE)
```
