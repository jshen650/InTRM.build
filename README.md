
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

This is a basic example which shows how to generate a data set from one
of the scenarios described in the paper and apply data splitting or the
m-out-of-n bootstrap approach with multiple imputation:

Data can be simulated with the `genDat` function. Data set size,
proportion of row-wise missingness, seed number, and missingness setting
number can be specified with “size,” “propMiss,” “seedNum,” and
“settingNumber,” respectively. Three missingness settings were
considered in our paper. Setting 1 refers to MAR missingness only in
$Y$; Setting 2 refers to MAR missingness in $Y$ and in covariate $X_4$;
and Setting 3 refers to MAR missingness in $Y$ and $X_4$ with MCAR
missingness in the other four covariates.

``` r
## Simulate data from one of the scenarios - say, n=300 with 30% missingness only in Y (Setting 1)
example_dat <- InTRM.build::genDat(size=300, propMiss = 30, seedNum = 215, settingNumber = 1)

str(example_dat) ## data set without and with missingness returned
#> List of 2
#>  $ fullDat:'data.frame': 300 obs. of  7 variables:
#>   ..$ Y  : num [1:300] -8.01 14.27 1.23 -9.03 9.9 ...
#>   ..$ A_2: num [1:300] -1 1 -1 1 -1 1 -1 -1 -1 -1 ...
#>   ..$ X1 : num [1:300] -5.228 -0.612 -0.743 -3.329 -3.196 ...
#>   ..$ X2 : num [1:300] -0.318 -1.958 -0.357 -2.257 0.226 ...
#>   ..$ X3 : num [1:300] -5.027 -2.312 0.53 -1.769 -0.905 ...
#>   ..$ X4 : num [1:300] -2.14 2.13 -5.01 -2.53 4.51 ...
#>   ..$ X5 : num [1:300] 4.995 3.6925 -0.8913 -0.0932 1.4059 ...
#>  $ missDat:'data.frame': 300 obs. of  7 variables:
#>   ..$ Y  : num [1:300] -8.01 NA NA -9.03 9.9 ...
#>   ..$ A_2: num [1:300] -1 1 -1 1 -1 1 -1 -1 -1 -1 ...
#>   ..$ X1 : num [1:300] -5.228 -0.612 -0.743 -3.329 -3.196 ...
#>   ..$ X2 : num [1:300] -0.318 -1.958 -0.357 -2.257 0.226 ...
#>   ..$ X3 : num [1:300] -5.027 -2.312 0.53 -1.769 -0.905 ...
#>   ..$ X4 : num [1:300] -2.14 2.13 -5.01 -2.53 4.51 ...
#>   ..$ X5 : num [1:300] 4.995 3.6925 -0.8913 -0.0932 1.4059 ...

missDat <- example_dat[[2]] ## data set containing missingness
```

Value and the variance of the Value can be estimated with the data
splitting approach using `split_Val`. Alternatively, the $m$-out-of-$n$
bootstrap can be run with `mn_Val`, which allows for a user-specified
value of $\alpha$ that is used for estimating the size of the resamples,
$m$. To determine a recommended $\alpha$, the double bootstrap procedure
can be run with `mn_double_Val`.

``` r

## estimate Value and variance with the data splitting approach
## on missing data set, specifying 5 imputations
example_Val <- split_Val(missDat, numImp=5, seedNum = 215)
example_Val
#>   Value_IPW  Var_IPW lower_IPW upper_IPW Value_AIPW Var_AIPW lower_AIPW
#> 1  7.055345 2.498996  2.157402  11.95329   7.559298 2.346575   4.556921
#>   upper_AIPW
#> 1   10.56168

## estimate a (1-eta)% confidence interval with the m-out-of-n bootstrap on data set with missingness
## specifying r=5 imputations
## and performing 10 repetitions of the m-out-of-n bootstrap that estimates m with a specified alpha value
example_mBoot <- mn_Val(missDat, numImp=5, seedNum=215, reps=10, eta=0.05, alpha=0.025)
example_mBoot[[3]] ## AIPW estimate of Value
#> [1] 8.459346
example_mBoot[[4]] ## final 95% confidence interval from running m-out-of-n bootstrap on data sets imputed r=5 times
#>      0.025    0.975
#> 1 6.752622 9.626606


## to determine a recommended alpha, the double bootstrap procedure can be used for obtaining (1-eta)% coverage
# candidate values for alpha
alpha_see <- seq(from= 0.025, to=1 ,by=0.025)
# running the double bootstrap procedure on data sets with missingness by specifying the number of imputations,
# a vector of candidate values of alpha, the number of outer (B1) and conditional (B2) bootstraps to perform,
# and the number of cores to use with mclapply
example_db <- mn_double_Val(missDat, numImp=5, alpha_see=alpha_see, seedNum=215, B1=5, B2=5, eta=0.05, nCores=4 )
example_db
#>        alpha coverage
#> db_Res 0.275        1
```

<!-- NOTE TO JENNY: If you get an error when knitting about not being able to find a function, make sure it is included in NAMESPACE as export(functionName)

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
