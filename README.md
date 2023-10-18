
<!-- README.md is generated from README.Rmd. Please edit that file -->

# traj2

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/IcatianTrout/traj2.svg?branch=main)](https://travis-ci.com/IcatianTrout/traj2)
[![R-CMD-check](https://github.com/IcatianTrout/traj2/workflows/R-CMD-check/badge.svg)](https://github.com/IcatianTrout/traj2/actions)
<!-- badges: end -->

The goal of traj2 is to implements the three-step procedure proposed by
Leffondree et al. (2004) to identify clusters of individual longitudinal
trajectories. The procedure involves (1) calculating a number of
measures of change capturing various features of the trajectories; (2)
using a PCA-based dimensionality reduction algorithm to select a subset
of measures and (3) using the K-means clustering algorithm to identify
clusters of trajectories.

## Installation

You can install the development version of traj2 from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("IcatianTrout/traj2")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
#library(traj2)
#data = trajdata

#mes = step1measures(data, ID=TRUE)
#sel = step2selection(mes)

#sel$loadings

#sel2 = step2selection(mes,select=c(10,12,8,4))

#clus <- step3(sel2,nclusters=4)$partition

#plot.traj(clus)
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
