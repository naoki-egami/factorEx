factorEx: Design and Analysis for Factorial Experiments
=======================================================

Description: R package `factorEx` provides design-based and model-based estimators for the population average marginal component effects (the pAMCE) in factorial experiments, including conjoint analysis. The package also implements a series of recommendations offered in de la Cuesta, Egami, and Imai (2019+).

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

Example of Design-based Confirmatory Analysis
---------------------------------------------

We use the conjoint experiment data based the marginal population randomization design.

### Case 1. Target Profile Distributions as Combination of Marginal Distributions

When using marginal distributions, `target_dist` should be a list and each element should have a factor name. Within each list, a `numeric` vector should have the same level names as those in `data`.

``` r
## Load the package and data
library(factorEx)
data("OnoBurden")
OnoBurden_data_pr <- OnoBurden$OnoBurden_data_pr
target_dist_marginal <- OnoBurden$target_dist_marginal

target_dist_marginal
```

    ## $gender
    ##      Male    Female 
    ## 0.6778243 0.3221757 
    ## 
    ## $age
    ## 36 years old 44 years old 52 years old 60 years old 68 years old 
    ##   0.05020921   0.13807531   0.23012552   0.22594142   0.25104603 
    ## 76 years old 
    ##   0.10460251 
    ## 
    ## $family
    ## Single (never married)      Single (divorced)     Married (no child) 
    ##             0.07729469             0.03864734             0.12560386 
    ## Married (two children) 
    ##             0.75845411 
    ## 
    ## $race
    ##          White       Hispanic Asian American          Black 
    ##      0.6725664      0.1283186      0.0000000      0.1991150 
    ## 
    ## $experience
    ##      None   4 years   8 years  12 years 
    ## 0.1966527 0.2259414 0.1548117 0.4225941 
    ## 
    ## $party
    ## Dem Rep 
    ##   1   0 
    ## 
    ## $pos_security
    ##     Cut military budget Maintain strong defense 
    ##              0.98557692              0.01442308

We can estimate the pAMCE with `design_pAMCE` with `target_type = "marginal"`.

``` r
out_design_mar <- 
  design_pAMCE(formula = Y ~ gender + age + family + race + experience + pos_security,
               data = OnoBurden_data_pr,
               pair_id = OnoBurden_data_pr$pair_id,
               cluster_id = OnoBurden_data_pr$id,
               target_dist  = target_dist_marginal, target_type = "marginal")
```

    ## Estimaing the pAMCEs for gender...age...family...race...experience...pos_security...

``` r
summary(out_design_mar)
```

    ## 
    ## ----------------
    ## Population AMCEs:
    ## ----------------
    ##  target_dist       factor                   level    Estimate  Std. Error
    ##       target       gender                  Female  0.02477776 0.006076128
    ##       target          age            44 years old  0.07306318 0.014829073
    ##       target          age            52 years old -0.01754745 0.014031192
    ##       target          age            60 years old  0.05025204 0.014063363
    ##       target          age            68 years old  0.06144573 0.013900972
    ##       target          age            76 years old -0.02126525 0.015390251
    ##       target       family      Married (no child) -0.01260551 0.012630608
    ##       target       family  Married (two children)  0.01366118 0.010457583
    ##       target       family       Single (divorced)  0.03244985 0.017206778
    ##       target         race                   Black -0.03719006 0.007298216
    ##       target         race                Hispanic -0.02519223 0.008688875
    ##       target   experience                12 years  0.05751447 0.007740096
    ##       target   experience                 4 years  0.02368536 0.008638306
    ##       target   experience                 8 years  0.05293234 0.009273195
    ##       target pos_security Maintain strong defense  0.04101219 0.023043993
    ##  p value    
    ##    0.000 ***
    ##    0.000 ***
    ##    0.211    
    ##    0.000 ***
    ##    0.000 ***
    ##    0.167    
    ##    0.318    
    ##    0.191    
    ##    0.059   .
    ##    0.000 ***
    ##    0.004  **
    ##    0.000 ***
    ##    0.006  **
    ##    0.000 ***
    ##    0.075   .
    ## ---
    ## Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Case 2. Target Profile Distributions as Combination of Marginal and Partial Joint Distributions

The use of partial joint distributions is useful because it can relax the assumption of no three-way or higher-order interactions (see de la Cuesta, Egami, and Imai (2019+)).

When using a combination of marginal and partial joint distributions, `target_dist` should be a list and each element should be a `numeric` vector (if marginal) or an `array`/`table` (if partial joint). Then, use argument `partial_joint_name` to specify which factors are marginal and partial joints. In the following example, `c("gender", "age", "family")` has the partial joint distributions over the three factors. `race` and `party` are based on the marginal distributions, respectively. `c("experience", "pos_security")` has the partial joint distributions over the two factors. Within each list, a `numeric` vector or an `array`/`table` should have the same level names as those in `data`.

``` r
target_dist_partial <- OnoBurden$target_dist_partial
target_dist_partial
```

    ## $`gender:age:family`
    ## , , family = Single (never married)
    ## 
    ##         age
    ## gender   36 years old 44 years old 52 years old 60 years old 68 years old
    ##   Male    0.004184100  0.004184100  0.004184100  0.004184100  0.008368201
    ##   Female  0.000000000  0.004184100  0.004184100  0.008368201  0.016736402
    ##         age
    ## gender   76 years old
    ##   Male    0.004184100
    ##   Female  0.004184100
    ## 
    ## , , family = Single (divorced)
    ## 
    ##         age
    ## gender   36 years old 44 years old 52 years old 60 years old 68 years old
    ##   Male    0.004184100  0.000000000  0.004184100  0.004184100  0.004184100
    ##   Female  0.000000000  0.004184100  0.004184100  0.004184100  0.000000000
    ##         age
    ## gender   76 years old
    ##   Male    0.000000000
    ##   Female  0.004184100
    ## 
    ## , , family = Married (no child)
    ## 
    ##         age
    ## gender   36 years old 44 years old 52 years old 60 years old 68 years old
    ##   Male    0.008368201  0.008368201  0.025104603  0.008368201  0.029288703
    ##   Female  0.000000000  0.000000000  0.004184100  0.008368201  0.012552301
    ##         age
    ## gender   76 years old
    ##   Male    0.000000000
    ##   Female  0.004184100
    ## 
    ## , , family = Married (two children)
    ## 
    ##         age
    ## gender   36 years old 44 years old 52 years old 60 years old 68 years old
    ##   Male    0.025104603  0.079497908  0.117154812  0.092050209  0.092050209
    ##   Female  0.004184100  0.020920502  0.050209205  0.062761506  0.033472803
    ##         age
    ## gender   76 years old
    ##   Male    0.041841004
    ##   Female  0.037656904
    ## 
    ## 
    ## $race
    ##          White       Hispanic Asian American          Black 
    ##      0.6725664      0.1283186      0.0000000      0.1991150 
    ## 
    ## $party
    ## Dem Rep 
    ##   1   0 
    ## 
    ## $`experience:pos_security`
    ##           pos_security
    ## experience Cut military budget Maintain strong defense
    ##   None             0.066945607             0.000000000
    ##   4 years          0.221757322             0.004184100
    ##   8 years          0.154811715             0.000000000
    ##   12 years         0.414225941             0.008368201

``` r
partial_joint_name  <- list(c("gender", "age", "family"), "race", "party", c("experience", "pos_security"))
```

We can estimate the pAMCE with `design_pAMCE` with `target_type = "partial_joint"` and appropriate `partial_joint_name`. The function can use `factor_name` to specify for which factors we estimate the pAMCE.

``` r
out_design_par <- 
      design_pAMCE(formula = Y ~ gender + age + family + race + experience + party + pos_security,
                   factor_name = c("gender", "race"),
                   data = OnoBurden_data_pr,
                   pair_id = OnoBurden_data_pr$pair_id,
                   cluster_id = OnoBurden_data_pr$id,
                   target_dist  = target_dist_partial, target_type = "partial_joint",
                   partial_joint_name = partial_joint_name)
```

    ## Estimaing the pAMCEs for gender...race...

``` r
summary(out_design_par, factor_name = c("gender", "race"))
```

    ## 
    ## ----------------
    ## Population AMCEs:
    ## ----------------
    ##  target_dist factor    level    Estimate  Std. Error p value    
    ##       target gender   Female  0.03123110 0.006559877   0.000 ***
    ##       target   race    Black -0.03499248 0.008122832   0.000 ***
    ##       target   race Hispanic -0.02745158 0.009568282   0.004  **
    ## ---
    ## Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
