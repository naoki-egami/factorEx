factorEx: Design and Analysis for Factorial Experiments
=======================================================

**Description:**

R package `factorEx` provides design-based and model-based estimators for the population average marginal component effects (the pAMCE) in factorial experiments, including conjoint analysis. The package also implements a series of recommendations offered in de la Cuesta, Egami, and Imai (2019+) and Egami and Imai (2019, JASA).

**Authors:**

-   [Naoki Egami](https://scholar.princeton.edu/negami/)
-   [Brandon de la Cuesta](https://www.brandondelacuesta.com/)
-   [Kosuke Imai](https://imai.fas.harvard.edu/)

**References:**

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

Examples
--------

-   [Design-based Confirmatory Analysis](#design)
-   [Case 1: Use Marginal Distributions for Target Profile Distribution](#designMar)
-   [Case 2: Use Combination of Marginal and Partial Joint Distributions for Target Profile Distribution](#designJoint)

-   [Model-based Exploratory Analysis](#model)

(1) Design-based Confirmatory Analysis
--------------------------------------

Here, we use the conjoint experiment that randomized profiles according to the marginal population randomization design.

### Case 1: Use Marginal Distributions for Target Profile Distributions

When using marginal distributions, `target_dist` should be a list and each element should have a factor name. Within each list, a `numeric` vector should have the same level names as those in `data`.

``` r
## Load the package and data
library(factorEx)
data("OnoBurden")

OnoBurden_data_pr <- OnoBurden$OnoBurden_data_pr # randomization based on marginal population design

# we focus on target profile distributions based on Democratic legislators. 
# See de la Cuesta, Egami, and Imai (2019+) for details.
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

We can estimate the pAMCE with `design_pAMCE` with `target_type = "marginal"`. Use `factor_name` to specify for which factors we estimate the pAMCE.

``` r
out_design_mar <- 
  design_pAMCE(formula = Y ~ gender + age + family + race + experience + party + pos_security,
               factor_name = c("gender", "age", "experience"),
               data = OnoBurden_data_pr,
               pair_id = OnoBurden_data_pr$pair_id,
               cluster_id = OnoBurden_data_pr$id,
               target_dist  = target_dist_marginal, target_type = "marginal")
```

    ## Estimaing the pAMCEs for gender...age...experience...

``` r
summary(out_design_mar)
```

    ## 
    ## ----------------
    ## Population AMCEs:
    ## ----------------
    ##  target_dist     factor        level     Estimate  Std. Error p value    
    ##       target     gender       Female  0.027987587 0.005861738   0.000 ***
    ##       target        age 44 years old  0.019219282 0.014421828   0.183    
    ##       target        age 52 years old -0.008792916 0.013765415   0.523    
    ##       target        age 60 years old -0.006826945 0.013875303   0.623    
    ##       target        age 68 years old  0.011247969 0.013569292   0.407    
    ##       target        age 76 years old -0.052741541 0.014775629   0.000 ***
    ##       target experience     12 years  0.041672460 0.007627281   0.000 ***
    ##       target experience      4 years  0.046173813 0.008868432   0.000 ***
    ##       target experience      8 years  0.040752213 0.009313376   0.000 ***
    ## ---
    ## Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Use `plot` to visualize the estimated pAMCEs.

``` r
plot(out_design_mar, factor_name = c("gender", "experience"))
```

### Case 2: Use Combination of Marginal and Partial Joint Distributions for Target Profile Distribution

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
                   factor_name = c("gender", "age", "race"),
                   data = OnoBurden_data_pr,
                   pair_id = OnoBurden_data_pr$pair_id,
                   cluster_id = OnoBurden_data_pr$id,
                   target_dist  = target_dist_partial, target_type = "partial_joint",
                   partial_joint_name = partial_joint_name)
```

    ## Estimaing the pAMCEs for gender...age...race...

``` r
summary(out_design_par)
```

    ## 
    ## ----------------
    ## Population AMCEs:
    ## ----------------
    ##  target_dist factor        level     Estimate  Std. Error p value    
    ##       target gender       Female  0.024756315 0.006362147   0.000 ***
    ##       target    age 44 years old  0.024750351 0.015045579   0.100    
    ##       target    age 52 years old -0.006198274 0.014335803   0.665    
    ##       target    age 60 years old -0.001011886 0.014397430   0.944    
    ##       target    age 68 years old  0.016337413 0.014132614   0.248    
    ##       target    age 76 years old -0.046107728 0.015464360   0.003  **
    ##       target   race        Black -0.025770076 0.008043842   0.001  **
    ##       target   race     Hispanic -0.028217748 0.009332710   0.002  **
    ## ---
    ## Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(2) Model-based Exploratory Analysis
------------------------------------

Here, we use the conjoint experiment that randomized profiles according to the uniform distribution and incorporate the target profile distribution in the analysis stage.

``` r
OnoBurden_data <- OnoBurden$OnoBurden_data # randomization based on uniform 

# due to large sample size, focus on "congressional candidates" for this example
OnoBurden_data_cong <- OnoBurden_data[OnoBurden_data$office == "Congress", ]

out_model <- 
      model_pAMCE(formula = Y ~ gender + age + family + race + experience + party + pos_security,
                   data = OnoBurden_data_cong, 
                  reg =  TRUE,
                   pair_id = OnoBurden_data_cong$pair_id,
                   cluster_id = OnoBurden_data_cong$id,
                   target_dist  = target_dist_marginal, target_type = "marginal")
```

    ## Note: suggest 'boot' greater than 500 for final results
    ## Cross-Validation: 20%..40%..60%..80%..100%..
    ## Bootstrap (100):

``` r
summary(out_model, factor_name = c("gender"))
```

    ## 
    ## ----------------
    ## Population AMCEs:
    ## ----------------
    ##  target_dist factor  level   Estimate Std. Error p value 
    ##     target_1 gender Female 0.02485328 0.01783633   0.163 
    ## ---
    ## Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

When `sample = TRUE`, the function also reports the AMCE based on the in-sample profile distributions (`sample AMCE`), which is the uniform AMCE in this example.

``` r
summary(out_model, factor_name = c("gender"), sample = TRUE)
```

    ## 
    ## ----------------
    ## Population AMCEs:
    ## ----------------
    ##  target_dist factor  level     Estimate  Std. Error p value 
    ##  sample AMCE gender Female -0.002290771 0.008321458   0.783 
    ##     target_1 gender Female  0.024853283 0.017836332   0.163 
    ## ---
    ## Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Use `plot` to visualize the estimated pAMCEs. When `diagnose = TRUE`, it provides two diagnostic checks; specification tests and the check of bootstrap distributions.

``` r
plot(out_model, factor_name = c("gender"), diagnose = TRUE)
```

In the model-based analysis, we can also decompose the difference between the pAMCE and the uniform AMCE. Use `effect_name` to specify which pAMCE we want to decompose. `effect_name` has two elements; the first is a factor name and the second is a level name of interest.

``` r
decompose_pAMCE(out_model, effect_name = c("gender", "Female"))
```

    ##                type       factor      estimate          se      low.95ci
    ## 1 target_1 - sample          age -4.476601e-03 0.002526321 -9.542598e-03
    ## 2 target_1 - sample       family -1.028249e-03 0.002956149 -6.693040e-03
    ## 3 target_1 - sample         race  5.505289e-03 0.007778271 -9.474694e-03
    ## 4 target_1 - sample   experience  6.965927e-05 0.000791340 -1.264228e-03
    ## 5 target_1 - sample        party  1.061621e-02 0.007640463 -6.015014e-03
    ## 6 target_1 - sample pos_security  1.685740e-02 0.008586040 -2.671033e-05
    ##       high.95ci
    ## 1 -0.0003348427
    ## 2  0.0058319109
    ## 3  0.0219185310
    ## 4  0.0015099003
    ## 5  0.0247779725
    ## 6  0.0315124642

Or use `plot_decompose` to visualize the decomposition.

``` r
plot_decompose(out_model, effect_name = c("gender", "Female"))
```
