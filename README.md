# microViz

Pre-processing, visualization, testing and interpretation tools for whole genome differential gene expression analyses

The `microViz` package provides tools for data pre-processing and variable selection (e.g. by Gini index or variance cutoff), testing for differences between dataset subsets in ahigh throughput manner (t test, Wilcoxon test, ANOVA and negative binomial modeling) and interpretation of the comparison results (subset marker finding and distances between the data subsets).

## Installation

You may easily fetch the package via `devtools`: 

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

Many thanks to authors, maintainers and contributors of the [tidyverse evironment](https://www.tidyverse.org/), and packages [furrr](https://furrr.futureverse.org/) and [future](https://www.futureverse.org/packages-overview.html), [ggrepel](https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html), [ggwordcloud](https://lepennec.github.io/ggwordcloud/), [limma](https://kasperdanielhansen.github.io/genbioconductor/html/limma.html), [MASS](https://cran.r-project.org/web/packages/MASS/index.html), [pROC](https://github.com/cran/pROC/tree/master), [rcompanion](https://rcompanion.org/handbook/), [stringi](https://stringi.gagolewski.com/index.html), [utils](https://cran.r-project.org/web/packages/R.utils/index.html), and [rlang](https://rlang.r-lib.org/).

## Usage

### Breast cancer transcriptome

<details>

Let's start with some example whole genome data provided with the `microViz` package. The `brca` data set stores normalized gene counts obtained by RNA sequencing of more than 150 breast carcinomas from [a patient-intitiated study]( https://www.cbioportal.org/study/summary?id=brca_mbcproject_2022). The data set comes along with clinical information concerning time point of sampling, metastatasis, histology and estrogen receptor (ER) status. In this vignette, we'll use both the untransformed counts as well as counts following `log2(x + 1)` transformation:

```r
## some tools

  library(tidyverse)
  library(rlang)
  library(microViz)
  library(trafo)

  library(org.Hs.eg.db)
  library(AnnotationDbi)

  ## for parallelization

  library(furrr)

  select <- dplyr::select
  reduce <- purrr::reduce

## the data set

  data("brca")

  counts <- brca

  ## listing all available genes

  genes <- names(brca)

  genes <- genes[!genes %in% c('sample_id',
                               'patient_id',
                               'timepoint',
                               'metastasis',
                               'histology',
                               'er_status')]

  ## a data frame with log2-transformed variables
  ## increasing by 1 to avoid potential log2(zero)

  log_expression <- counts

  log_expression[genes] <- log_expression[genes] %>%
    map_dfc(~log2(.x + 1))

```
```r
> head(counts)
# A tibble: 6 × 38,957
  sample…¹ patie…² timep…³ metas…⁴ histo…⁵ er_st…⁶ TSPAN6  TNMD  DPM1 SCYL3 C1orf…⁷   FGR   CFH FUCA2  GCLC  NFYA STPG1 NIPAL3 LAS1L ENPP4
  <chr>    <chr>   <fct>   <fct>   <fct>   <fct>    <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl>
1 MBC-MBC… MBCPro… T2      METAST… IDC     NA        3.36  0.17  17.3  6.44    8.23  5.09 25.1   6.11  25.3  17.4  2.3   11.2   16.8  2.05
2 MBC-MBC… MBCPro… T1      METAST… MIXED_… NEGATI…   5.58  0.61  34.4 17.4    22.4   1.71 35.6   7.66  33.9  32.9  1.54   6.79  17.6  2.04
3 MBC-MBC… MBCPro… T2      METAST… ILC     NA        5.42  0     34.8 13.5    13.1   4.14 17.4   3.51  33.9  16.6  0.85   5.74  13.8  1.11
4 MBC-MBC… MBCPro… T1      METAST… IDC     POSITI…   2.96  0.01  30.1  8.5    27.2   0.75  9.74  6.27  37.9  16.2  2.6    6.56  17.9  3.11
5 MBC-MBC… MBCPro… T3      METAST… IDC     NA        4.47  0.55  42.7  8.61   13.4   2.94 51.5   8.02  35.8  10.1  2.14  10.9   12.7  2.87
6 MBC-MBC… MBCPro… T1      METAST… IDC     POSITI…   5.41  0.26  25.9  4.26    6.18  2.49 48.0   7.97  26.9  18.0  2.39  14.6   13.4  2.76
# … with 38,937 more variables: SEMA3F <dbl>, CFTR <dbl>, ANKIB1 <dbl>, CYP51A1 <dbl>, KRIT1 <dbl>, RAD52 <dbl>, MYH16 <dbl>, BAD <dbl>,
#   LAP3 <dbl>, CD99 <dbl>, HS3ST1 <dbl>, AOC1 <dbl>, WNT16 <dbl>, HECW1 <dbl>, MAD1L1 <dbl>, LASP1 <dbl>, SNX11 <dbl>, TMEM176A <dbl>,
#   M6PR <dbl>, KLHL13 <dbl>, CYP26B1 <dbl>, ICA1 <dbl>, DBNDD1 <dbl>, ALS2 <dbl>, CASP10 <dbl>, CFLAR <dbl>, TFPI <dbl>, NDUFAF7 <dbl>,
#   RBM5 <dbl>, MTMR7 <dbl>, SLC7A2 <dbl>, ARF5 <dbl>, SARM1 <dbl>, POLDIP2 <dbl>, PLXND1 <dbl>, AK2 <dbl>, CD38 <dbl>, FKBP4 <dbl>,
#   KDM1A <dbl>, RBM6 <dbl>, CAMKK1 <dbl>, RECQL <dbl>, CCDC132 <dbl>, HSPB6 <dbl>, ARHGAP33 <dbl>, NDUFAB1 <dbl>, PDK4 <dbl>,
#   SLC22A16 <dbl>, ZMYND10 <dbl>, ABCB5 <dbl>, ARX <dbl>, SLC25A13 <dbl>, ST7 <dbl>, CDC27 <dbl>, SLC4A1 <dbl>, CALCR <dbl>, HCCS <dbl>,
#   DVL2 <dbl>, PRSS22 <dbl>, UPF1 <dbl>, SKAP2 <dbl>, SLC25A5 <dbl>, CCDC109B <dbl>, HOXA11 <dbl>, POLR2J <dbl>, DHX33 <dbl>, …
# ℹ Use `colnames()` to see all variable names

> head(log_expression)
# A tibble: 6 × 38,957
  sampl…¹ patie…² timep…³ metas…⁴ histo…⁵ er_st…⁶ TSPAN6   TNMD  DPM1 SCYL3 C1orf…⁷   FGR   CFH FUCA2  GCLC  NFYA STPG1 NIPAL3 LAS1L ENPP4
  <chr>   <chr>   <fct>   <fct>   <fct>   <fct>    <dbl>  <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl>
1 MBC-MB… MBCPro… T2      METAST… IDC     NA        2.12 0.227   4.20  2.90    3.21 2.61   4.71  2.83  4.72  4.20 1.72    3.61  4.15  1.61
2 MBC-MB… MBCPro… T1      METAST… MIXED_… NEGATI…   2.72 0.687   5.15  4.20    4.55 1.44   5.20  3.11  5.13  5.08 1.34    2.96  4.21  1.60
3 MBC-MB… MBCPro… T2      METAST… ILC     NA        2.68 0       5.16  3.86    3.82 2.36   4.20  2.17  5.12  4.14 0.888   2.75  3.89  1.08
4 MBC-MB… MBCPro… T1      METAST… IDC     POSITI…   1.99 0.0144  4.96  3.25    4.82 0.807  3.42  2.86  5.28  4.11 1.85    2.92  4.24  2.04
5 MBC-MB… MBCPro… T3      METAST… IDC     NA        2.45 0.632   5.45  3.26    3.85 1.98   5.71  3.17  5.20  3.47 1.65    3.57  3.78  1.95
6 MBC-MB… MBCPro… T1      METAST… IDC     POSITI…   2.68 0.333   4.75  2.40    2.84 1.80   5.62  3.17  4.80  4.25 1.76    3.97  3.85  1.91
# … with 38,937 more variables: SEMA3F <dbl>, CFTR <dbl>, ANKIB1 <dbl>, CYP51A1 <dbl>, KRIT1 <dbl>, RAD52 <dbl>, MYH16 <dbl>, BAD <dbl>,
#   LAP3 <dbl>, CD99 <dbl>, HS3ST1 <dbl>, AOC1 <dbl>, WNT16 <dbl>, HECW1 <dbl>, MAD1L1 <dbl>, LASP1 <dbl>, SNX11 <dbl>, TMEM176A <dbl>,
#   M6PR <dbl>, KLHL13 <dbl>, CYP26B1 <dbl>, ICA1 <dbl>, DBNDD1 <dbl>, ALS2 <dbl>, CASP10 <dbl>, CFLAR <dbl>, TFPI <dbl>, NDUFAF7 <dbl>,
#   RBM5 <dbl>, MTMR7 <dbl>, SLC7A2 <dbl>, ARF5 <dbl>, SARM1 <dbl>, POLDIP2 <dbl>, PLXND1 <dbl>, AK2 <dbl>, CD38 <dbl>, FKBP4 <dbl>,
#   KDM1A <dbl>, RBM6 <dbl>, CAMKK1 <dbl>, RECQL <dbl>, CCDC132 <dbl>, HSPB6 <dbl>, ARHGAP33 <dbl>, NDUFAB1 <dbl>, PDK4 <dbl>,
#   SLC22A16 <dbl>, ZMYND10 <dbl>, ABCB5 <dbl>, ARX <dbl>, SLC25A13 <dbl>, ST7 <dbl>, CDC27 <dbl>, SLC4A1 <dbl>, CALCR <dbl>, HCCS <dbl>,
#   DVL2 <dbl>, PRSS22 <dbl>, UPF1 <dbl>, SKAP2 <dbl>, SLC25A5 <dbl>, CCDC109B <dbl>, HOXA11 <dbl>, POLR2J <dbl>, DHX33 <dbl>, …
# ℹ Use `colnames()` to see all variable names
```
</details>

### Column and row statistics and variable selection

<details>
  
The package provides a wode range of statistics of central tendency and distribution which somehow are not provided in base R. For sake of speed, they all operate in C++ under the hood:

```r
## geometric mean

>   Gmean(counts$HERC2)
[1] 36.00704

## harmonic mean
> Hmean(counts$HERC2)
[1] 31.7873

## 95% percentile confidence interval
> perCI(counts$HERC2)
[1] 15.002 76.100

## 95% confidence interval computed with the BCA method
> bcaCI(counts$HERC2)
[1] 15.18945 81.06810

## Gini coefficient: with and without sample bias correction

> Gini(counts$HERC2, unbiased = TRUE)
[1] 0.2364145

> Gini(counts$HERC2, unbiased = FALSE)
[1] 0.2349086

## ratio of frequencies of the first most common to the semond most common element
## of a vector
> freqRatio(counts$TSPAN6) 
[1] 1.5

## percentage of unique values
> percUnique(counts$TSPAN6)
[1] 91.71975
```
Each of them has a column- and row-wise counterpart, which may be useful at selection of genes expressed with sufficient variability for differential gene expression analysis of modeling:

```r
## numeric stats for 38K+ genes, here for variable medians:

>   system.time(counts[genes] %>%
+                 colMedians(na.rm = TRUE))
   user  system elapsed 
   0.74    0.03    0.77

## variable minima and maxima

>   colMins(counts[genes]) %>% head
  TSPAN6     TNMD     DPM1    SCYL3 C1orf112      FGR 
    0.32     0.00     0.04     1.03     2.15     0.00

>   counts[genes] %>% colMax %>% head
  TSPAN6     TNMD     DPM1    SCYL3 C1orf112      FGR 
   19.31    13.69    89.72    20.25    57.76    30.55

## geometric and harmonic means

> counts[genes] %>% colGmeans %>% head
   TSPAN6      TNMD      DPM1     SCYL3  C1orf112       FGR 
 4.297067  0.000000 28.100647  6.329931 13.418881  0.000000

> counts[genes] %>% colHmeans %>% head
   TSPAN6      TNMD      DPM1     SCYL3  C1orf112       FGR 
 3.113677  0.000000  5.067529  5.551144 11.331793  0.000000

## variances and Gini coefficients

> counts[genes] %>% colVars %>% head
    TSPAN6       TNMD       DPM1      SCYL3   C1orf112        FGR 
 10.557030   2.244567 211.335627  11.319315  88.910633   9.789100

>  counts[genes] %>% colGini %>% head
   TSPAN6      TNMD      DPM1     SCYL3  C1orf112       FGR 
0.3292340 0.6346546 0.2448090 0.2575666 0.3083520 0.4285603

## ration of frequencies of the most common to the second most common element
## and percentage of unique elements

> counts[genes] %>% colFreqRatios %>% head
  TSPAN6     TNMD     DPM1    SCYL3 C1orf112      FGR 
     1.5      4.5      1.0      1.5      1.0      1.0

>  counts[genes] %>% colPercUniques %>% head
  TSPAN6     TNMD     DPM1    SCYL3 C1orf112      FGR 
91.71975 64.33121 97.45223 94.90446 95.54140 80.25478 

```
The process of variable selection based e.g. on Gini coefficient, frequency ratio, or variance to mean ratio can be facilitated by the `distr_stats()` (variables in columns) and `row_stats()` (variables in rows). As such, those two functions generate ta similar output bunch of stats as `caret::nearZeroVar()()` but are expected to be several times faster:

```r
>   system.time(distr_stats(log_expression[genes]))
   user  system elapsed 
   5.22    0.25    5.48

>   system.time(caret::nearZeroVar(log_expression[genes]))
   user  system elapsed 
  22.92    0.44   23.38

```

```r
## selection of genes with a Gini coefficient cutoff
  
  colStats <- log_expression[genes] %>%
    distr_stats

  variant_genes <- colStats %>%
    filter(gini_coef >= 0.1,
           freqRatio < 5) %>%
    .$variable
```

```r
> head(colStats)

# A tibble: 6 × 9
  variable  mean   var gini_coef var_mean_ratio freqRatio percentUnique zeroVar nzv  
  <chr>    <dbl> <dbl>     <dbl>          <dbl>     <dbl>         <dbl> <lgl>   <lgl>
1 TSPAN6   2.47  0.569    0.171           0.231       1.5          91.7 FALSE   FALSE
2 TNMD     0.718 0.480    0.511           0.669       4.5          64.3 FALSE   FALSE
3 DPM1     4.90  0.550    0.0769          0.112       1            97.5 FALSE   FALSE
4 SCYL3    2.90  0.356    0.116           0.123       1.5          94.9 FALSE   FALSE
5 C1orf112 3.87  0.573    0.110           0.148       1            95.5 FALSE   FALSE
6 FGR      1.71  0.596    0.247           0.349       1            80.3 FALSE   FALSE

> head(variant_genes)
[1] "TSPAN6"   "TNMD"     "SCYL3"    "C1orf112" "FGR"      "CFH"

> length(genes)
[1] 38951

> length(variant_genes)
[1] 13749

```
</details>


