# erah

<!-- badges: start -->
[![R-CMD-check](https://github.com/xdomingoal/erah-devel/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/xdomingoal/erah-devel/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/xdomingoal/erah-devel/branch/master/graph/badge.svg)](https://app.codecov.io/gh/xdomingoal/erah-devel?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/erah)](https://cran.r-project.org/package=erah)
[![](http://cranlogs.r-pkg.org/badges/erah)](http://cran.rstudio.com/web/packages/erah/index.html)
<!-- badges: end -->

### Automated Spectral Deconvolution, Alignment, and Metabolite Identification in GC/MS-Based Untargeted Metabolomics

Automated compound deconvolution, alignment across samples, and identification of metabolites by spectral library matching in Gas Chromatography - Mass spectrometry (GC-MS) untargeted metabolomics. Outputs a table with compound names, matching scores and the integrated area of the compound for each sample.

[Domingo-Almenara, X., Brezmes, J., Vinaixa, M., Samino, S., Ramirez, N., Ramon-Krauel, M., Lerin, C., Díaz, M., Ibáñez, L., Correig, X. and Perera-Lluna, A., 2016. eRah: a computational tool integrating spectral deconvolution and alignment with quantification and identification of metabolites in GC/MS-based metabolomics. Analytical Chemistry, 88(19), pp.9821-9829.](https://doi.org/10.1021/acs.analchem.6b02927)

eRah development version. For downloads, please use CRAN
https://CRAN.R-project.org/package=erah 

## Installation

From an R console, execute:

```r
install.packages('erah')
```

To install this development version, execute:

```r
remotes::install_github('xdomingoal/erah-devel')
```

## Learn more

The package documentation can be browsed online at <https://xdomingoal.github.io/erah-devel/>. 

If this is your first time using `erah` see the [manual](https://xdomingoal.github.io/erah-devel/articles/erah.html) for information on how to get started.

If you believe you've found a bug in `erah`, please file a bug (and, if
possible, a [reproducible example](https://reprex.tidyverse.org)) at
<https://github.com/xdomingoal/erah-devel/issues>.
