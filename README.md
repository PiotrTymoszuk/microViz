# microViz

Pre-processing, visualization, testing and interpretation tools for whole genome differential gene expression analyses

The `microViz` package provides tools for data pre-processing and variable selection (e.g. by Gini index or variance cutoff), testing for differences between dataset subsets in ahigh throughput manner (t test, Wilcoxon test, ANOVA and negative binomial modeling) and interpretation of the comparison results (subset marker finding and distances between the data subsets).

## Installation

You may easily fetch the package `devtools`: 

```r

devtools::install_github('PiotrTymoszuk/microViz')

## suggested for use of cross-disstance object OOP:

devtools::install_github('PiotrTymoszuk/clustTools')

```

## Terms of use

The package is available under a [GPL-3 license](https://github.com/PiotrTymoszuk/microViz/blob/main/LICENSE).

## Contact

The package maintainer is [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com).

## Acknowledgements

Many thanks to authors, maintainers and contributors of the [tidyverse evironment](https://www.tidyverse.org/), and packages [caret](https://topepo.github.io/caret/), [coxed](https://cran.r-project.org/web/packages/coxed/index.html), [descTools](https://cran.r-project.org/web/packages/DescTools/index.html), [furrr](https://furrr.futureverse.org/) and [future](https://www.futureverse.org/packages-overview.html), [ggrepel](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html), [ggwordcloud](https://lepennec.github.io/ggwordcloud/), [limma](https://kasperdanielhansen.github.io/genbioconductor/html/limma.html), [MASS](https://cran.r-project.org/web/packages/MASS/index.html), [pROC](https://github.com/cran/pROC/tree/master), [rcompanion](https://rcompanion.org/handbook/), [stringi](https://stringi.gagolewski.com/index.html), [utils](https://cran.r-project.org/web/packages/R.utils/index.html), and [rlang](https://rlang.r-lib.org/).


