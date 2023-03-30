
<!-- README.md is generated from README.Rmd. Please edit that file -->

# InTRM.build

InTRM.build is a rough version of the package that will later be made
publicly available and referred to as InTRM. The Individualized
Treatment Rules following Multiple Imputation (InTRM) package can be
used for implementing the methodology described in our paper,
“Estimation and Evaluation of Individualized Treatment Rules following
Multiple Imputation.” The m-out-of-n bootstrap approach used in our work
was adapted from the original work by Chakraborty et al (2013).

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the development version of InTRM.build from …

``` r
# install.packages("devtools")
devtools::install_github("jshen650/InTRM.build")
```

## Example

This is a basic example which shows you how to simulate data and apply
one of the two approaches for estimating optimal decision rules
following multiple imputation:

``` r
## simulate data from one of the scenarios - say, n=300 with 30% missingness only in Y
example_dat <- InTRM.build::genDat(size=300, propMiss = 30, seedNum = 215, settingNumber = 1)

str(example_dat) ## data set without and with missingness returned
#> List of 2
#>  $ :'data.frame':    300 obs. of  8 variables:
#>   ..$ missBin: int [1:300] 0 1 1 0 0 0 0 1 1 0 ...
#>   ..$ Y      : num [1:300] -8.01 14.27 1.23 -9.03 9.9 ...
#>   ..$ A_2    : num [1:300] -1 1 -1 1 -1 1 -1 -1 -1 -1 ...
#>   ..$ X1     : num [1:300] -5.228 -0.612 -0.743 -3.329 -3.196 ...
#>   ..$ X2     : num [1:300] -0.318 -1.958 -0.357 -2.257 0.226 ...
#>   ..$ X3     : num [1:300] -5.027 -2.312 0.53 -1.769 -0.905 ...
#>   ..$ X4     : num [1:300] -2.14 2.13 -5.01 -2.53 4.51 ...
#>   ..$ X5     : num [1:300] 4.995 3.6925 -0.8913 -0.0932 1.4059 ...
#>  $ :'data.frame':    300 obs. of  7 variables:
#>   ..$ Y  : num [1:300] -8.01 NA NA -9.03 9.9 ...
#>   ..$ A_2: num [1:300] -1 1 -1 1 -1 1 -1 -1 -1 -1 ...
#>   ..$ X1 : num [1:300] -5.228 -0.612 -0.743 -3.329 -3.196 ...
#>   ..$ X2 : num [1:300] -0.318 -1.958 -0.357 -2.257 0.226 ...
#>   ..$ X3 : num [1:300] -5.027 -2.312 0.53 -1.769 -0.905 ...
#>   ..$ X4 : num [1:300] -2.14 2.13 -5.01 -2.53 4.51 ...
#>   ..$ X5 : num [1:300] 4.995 3.6925 -0.8913 -0.0932 1.4059 ...

missDat <- example_dat[[2]] ## data set containing missingness
```

``` r

## estimate Value and variance with the data splitting approach
example_Val <- split_Val(missDat, 5, seedNum = 215)
example_Val
#>   Value_IPW  Var_IPW Value_AIPW Var_AIPW
#> 1  7.055345 2.498996   7.559298 2.346575

## estimate a 95% confidence interval with the m-out-of-n bootstrap
example_mBoot <- run_mBoot(missDat, numImp=5, seedNum=215, reps=10, eta=0.05)
example_mBoot[[3]] ## final 95% confidence interval from running m-out-of-n bootstrap on data sets imputed r=5 times
#>      0.025    0.975
#> 1 6.752622 9.626606
```

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
