factorEx: Design and Analysis for Factorial Experiments
=======================================================

Authors:

-   [Naoki Egami](https://scholar.princeton.edu/negami/)
-   [Brandon de la Cuesta](https://www.brandondelacuesta.com//)
-   [Kosuke Imai](https://imai.fas.harvard.edu/)

References:

-   de la Cuesta, Egami, and Imai. (2019+). [Improving the External Validity of Conjoint Analysis: The Essential Role of Profile Distribution.](https://scholar.princeton.edu/sites/default/files/negami/files/conjoint_profile.pdf) (Working Paper)

-   Egami and Imai. (2019). [Causal Interaction in Factorial Experiments: Application to Conjoint Analysis.](https://scholar.princeton.edu/sites/default/files/negami/files/causalint.pdf) *Journal of the American Statistical Association*, Vol.114, No.526 (June), pp. 529–540.

Installation Instructions
-------------------------

`factorEx` is available on CRAN and can be installed using:

``` r
install.packages("factorEx")
```

You can also install the most recent development version using the `devtools` package. First you have to install `devtools` using the following code. Note that you only have to do this once:

``` r
if(!require(devtools)) install.packages("devtools")
```

Then, load `devtools` and use the function `install_github()` to install `factorEx`:

``` r
library(devtools)
install_github("naoki-egami/factorEx", dependencies=TRUE)
```

Example
-------

`factorEx` implements the design-based and model-based approaches to estimate the pAMCE.

### 1. Design-based Confirmatory Analysis

``` r
## Load the package and data
library(factorEx)
data("OnoBurden")
OnoBurden_data <- OnoBurden$OnoBurden_data
target_dist_marginal <- OnoBurden$target_dist_marginal

out_design_mar <- design_pAMCE(formula = Y ~ gender + age + family + race + experience + party + pos_security,
                           data = OnoBurden_data,
                           pair_id = OnoBurden_data$pair_id,
                           cluster_id = OnoBurden_data$id,
                           target_dist  = target_dist_marginal,
                           target_type = "marginal")
```

    ## Estimaing the pAMCEs for gender...age...family...race...experience...party...pos_security...

``` r
summary(out_design_mar, factor_name = "gender",  sample = TRUE)
```

    ## ----------------
    ## Population AMCEs:
    ## ----------------
    ##  target_dist factor  level    Estimate  Std. Error p value  
    ##       target gender Female  0.01782688 0.016824048   0.289  
    ##  sample AMCE gender Female -0.01250802 0.005752017   0.030 *
    ## ---
    ## Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
target_dist_partial <- OnoBurden$target_dist_partial
partial_joint_name  <- list(c("gender", "age", "family"), "race", "party", c("experience", "pos_security"))
out_design_par <- design_pAMCE(formula = Y ~ gender + age + family + race + experience + party + pos_security,
                           data = OnoBurden_data,
                           pair_id = OnoBurden_data$pair_id,
                           cluster_id = OnoBurden_data$id,
                           target_dist  = target_dist_partial,
                           target_type = "partial_joint",
                           partial_joint_name = partial_joint_name)
```

    ## Estimaing the pAMCEs for gender...age...family...race...experience...party...pos_security...

``` r
summary(out_design_par, factor_name = "gender",  sample = TRUE)
```

    ## ----------------
    ## Population AMCEs:
    ## ----------------
    ##  target_dist factor  level    Estimate  Std. Error p value  
    ##       target gender Female  0.01332267 0.018309857   0.467  
    ##  sample AMCE gender Female -0.01250802 0.005752017   0.030 *
    ## ---
    ## Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
## Load the package and data
library(factorEx)
data("OnoBurden")
OnoBurden_data_pr <- OnoBurden$OnoBurden_data_pr
target_dist_marginal <- OnoBurden$target_dist_marginal

out_design_mar_pr <- design_pAMCE(formula = Y ~ gender + age + family + race + experience + party + pos_security, 
                                  factor_name  = c("gender"),
                                  data = OnoBurden_data_pr,
                                  pair_id = OnoBurden_data_pr$pair_id,
                                  cluster_id = OnoBurden_data_pr$id,
                                  target_dist  = target_dist_marginal,
                                  target_type = "marginal")
```

    ## Estimaing the pAMCEs for gender...

``` r
summary(out_design_mar_pr, factor_name = "gender")
```

    ## ----------------
    ## Population AMCEs:
    ## ----------------
    ##  target_dist factor  level   Estimate  Std. Error p value    
    ##       target gender Female 0.02477776 0.006076128       0 ***
    ## ---
    ## Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(out_design_mar, factor_name = "gender")
```

    ## ----------------
    ## Population AMCEs:
    ## ----------------
    ##  target_dist factor  level   Estimate Std. Error p value 
    ##       target gender Female 0.01782688 0.01682405   0.289 
    ## ---
    ## Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
