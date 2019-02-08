
<!-- README.md is generated from README.Rmd. Please edit that file -->
ICcforest
=========

The goal of ICcforest is to implement the conditional inference forest approach to modeling interval-censored survival data. It also provides functions to tune the parameters and evaluate the model fit.

Installation
------------

You can install the released version of ICcforest from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("ICcforest")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code with miceData
library(ICcforest)
library(survival)
library(icenReg)
#> Loading required package: Rcpp
#> Loading required package: coda
data(miceData)

## For ICcforest to run, Inf should be set to be a large number, for example, 9999999.
idx_inf <- (miceData$u == Inf)
miceData$u[idx_inf] <- 9999999.

## Fit an iterval-censored conditional inference forest
Cforest <- ICcforest(Surv(l, u, type = "interval2") ~ grp, data = miceData)
#> mtry = 1  OOB Brier score = 0.06497173 
#> Searching left ...
#> Searching right ...
```
