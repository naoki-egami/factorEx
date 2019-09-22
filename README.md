fastLink: Fast Probabilistic Record Linkage [![Build Status](https://travis-ci.org/kosukeimai/fastLink.svg?branch=master)](https://travis-ci.org/kosukeimai/fastLink) [![CRAN Version](http://www.r-pkg.org/badges/version/fastLink)](https://CRAN.R-project.org/package=fastLink) ![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/fastLink)
===========================================================================================================================================================================================================================================================================================================================================================

Authors:

-   [Ted Enamorado](https://www.tedenamorado.com/)
-   [Ben Fifield](https://www.benfifield.com/)
-   [Kosuke Imai](https://imai.fas.harvard.edu/)

For a detailed description of the method see:

-   [Using a Probabilistic Model to Assist Merging of Large-scale Administrative Records](http://imai.fas.harvard.edu/research/files/linkage.pdf)

Applications of the method:

-   [Validating Self-reported Turnout by Linking Public Opinion Surveys with Administrative Records](http://imai.fas.harvard.edu/research/files/turnout.pdf)

Technical reports:

-   [User’s Guide and Codebook for the ANES 2016 Time Series Voter Validation Supplemental Data](https://www.electionstudies.org/wp-content/uploads/2018/03/anes_timeseries_2016voteval_userguidecodebook.pdf)

-   [User’s Guide and Codebook for the CCES 2016 Voter Validation Supplemental Data](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/2NNA4L)

Data:

-   [ANES 2016 Time Series Voter Validation Supplemental Data](http://www.electionstudies.org/studypages/download/datacenter_all_NoData.php)

-   [CCES 2016 Voter Validation Supplemental Data](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/2NNA4L)

Installation Instructions
-------------------------

`fastLink` is available on CRAN and can be installed using:

``` r
install.packages("fastLink")
```

You can also install the most recent development version of `fastLink` using the `devtools` package. First you have to install `devtools` using the following code. Note that you only have to do this once:

``` r
if(!require(devtools)) install.packages("devtools")
```

Then, load `devtools` and use the function `install_github()` to install `fastLink`:

``` r
library(devtools)
install_github("kosukeimai/fastLink",dependencies=TRUE)
```

Simple usage example
--------------------

The linkage algorithm can be run either using the `fastLink()` wrapper, which runs the algorithm from start to finish, or step-by-step. We will outline the workflow from start to finish using both examples. In both examples, we have two dataframes called `dfA` and `dfB` that we want to merge together, and they have seven commonly named fields:

-   `firstname`

-   `middlename`

-   `lastname`

-   `housenum`

-   `streetname`

-   `city`

-   `birthyear`

### Running the algorithm using the `fastLink()` wrapper

The `fastLink` wrapper runs the entire algorithm from start to finish, as seen below:

``` r
## Load the package and data
library(fastLink)
data(samplematch)

matches.out <- fastLink(
  dfA = dfA, dfB = dfB, 
  varnames = c("firstname", "middlename", "lastname", "housenum", "streetname", "city", "birthyear"),
  stringdist.match = c("firstname", "middlename", "lastname", "streetname", "city"),
  partial.match = c("firstname", "lastname", "streetname")
)
```

    ## 
    ## ==================== 
    ## fastLink(): Fast Probabilistic Record Linkage
    ## ==================== 
    ## 
    ## Calculating matches for each variable.
    ## Getting counts for parameter estimation.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Running the EM algorithm.
    ## Getting the indices of estimated matches.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Deduping the estimated matches.
    ## Getting the match patterns for each estimated match.

-   `varnames` should be a vector of variable names to be used for matching. These variable names should exist in both `dfA` and `dfB`

-   `stringdist.match` should be a vector of variable names present in `varnames`. For those variables included in `stringdist.match`, agreement will be calculated using Jaro-Winkler distance.

-   `partial.match` is another vector of variable names present in both `stringdist.match` and `varnames`. A variable included in `partial.match` will have a partial agreement category calculated in addition to disagreement and absolute agreement, as a function of Jaro-Winkler distance.

Other arguments that can be provided include:

-   `cut.a`: Lower bound for full string-distance match, ranging between 0 and 1. Default is 0.92.

-   `cut.p`: Lower bound for partial string-distance match, ranging between 0 and 1. Default is 0.88.

-   `priors.obj`: The output from `calcMoversPriors()`, allowing the inclusion of auxiliary information on moving behavior to aid matching. We will discuss this option further at the end of this vignette.

-   `w.lambda`: The user-specified weighting of the MLE and prior estimate for the *λ* parameter, a number between 0 and 1. We will discuss this option further at the end of this vignette.

-   `w.pi`: The user-specified weighting of the MLE and prior estimate for the *π* parameter, a number between 0 and 1. We will discuss this option further at the end of this vignette.

-   `address.field`: The name of the address field, to be specified when providing a prior on the probability of moving in-state through `priors.obj`. The variable listed in `address.field` must be listed in `varnames`. We will discuss this option further at the end of this vignette.

-   `gender.field`: The name of the gender field, if matching on gender. If provided, the EM algorithm will implement a prior that enforces near-perfect blocking on gender, so that no matches that disagree on gender will be in the matched set. Can be used in conjunction with movers priors, if the user does not want to specify the same prior for both genders when blocking.

-   `estimate.only`: Whether to stop running the algorithm after running the EM estimation step. Can be used when running the algorithm on a random sample, and then applying those estimates to the full data set.

-   `em.obj`: An EM object, either from an `estimate.only = TRUE` run of `fastLink` or from `emlinkMARmov()`. If provided, the algorithm will skip the EM estimation step and proceed to apply the estimates from the EM object to the full data set. To be used when the EM has been estimated on a random sample of data and should be applied to the full data set.

-   `dedupe.matches`: Whether to dedupe the matches returned by the algorithm, ensuring that each observation in dataset A is matched to at most one observation in dataset B (and vice versa). Can be done either using Winkler's linear assignment solution (recommended) or by iteratively selecting the maximum posterior value for a given observation (if N size makes linear assignment solution prohibitively slow). Default is `TRUE`.

-   `linprog.dedupe`: Whether to use Winkler's linear programming solution to the deduplication problem (recommended when N size is not prohibitively large). Default is `FALSE`.

-   `n.cores`: The number of registered cores to parallelize over. If left unspecified. the function will estimate this on its own.

-   `tol.em`: Convergence tolerance for the EM algorithm. Default is 1e-04

-   `threshold.match`: Lower bound for the posterior probability of a match that will be accepted. Default is 0.85.

-   `verbose`: Whether to print out runtime for each step and EM output. Default is FALSE.

The output from `fastLink()` when `estimate.only = FALSE` will be a list of length 4 with two entries:

-   `matches`: A matrix where each row is a match with the relevant indices of `dfA` (column 1) and `dfB` (column 2).

-   `EM`: The output from the EM algorithm.

-   `nobs.a`: The number of observations in dataset A.

-   `nobs.b`: The number of observations in dataset B.

When `estimate.only = TRUE`, `fastLink()` outputs the EM object.

The datasets can then be subsetted down to the matches as follows:

``` r
dfA.match <- dfA[matches.out$matches$inds.a,]
dfB.match <- dfB[matches.out$matches$inds.b,]
```

or using the `getMatches()` function:

``` r
matched_dfs <- getMatches(
  dfA = dfA, dfB = dfB, 
  fl.out = matches.out, threshold.match = 0.85
)
```

We can also examine the EM object:

``` r
matches.out$EM
```

    ## $zeta.j
    ##                        [,1]
    ##  [1,] 4.865161101204859e-50
    ##  [2,] 8.798168149596711e-46
    ##  [3,] 5.821547331784233e-40
    ##  [4,] 2.545028663970180e-42
    ##  [5,] 4.602435492959898e-38
    ##  [6,] 3.045326664389509e-32
    ##  [7,] 1.704645594591585e-44
    ##  [8,] 3.082684882309823e-40
    ##  [9,] 2.039742324334265e-34
    ## [10,] 3.216699987320368e-43
    ## [11,] 1.127061026026548e-37
    ## [12,] 1.269642525513900e-37
    ## [13,] 2.296024365274058e-33
    ## [14,] 6.641664177632459e-30
    ## [15,] 4.448548553278978e-32
    ## [16,] 2.647658248164910e-41
    ## [17,] 9.276812986933977e-36
    ## [18,] 1.110043527195273e-25
    ## [19,] 1.540835829269232e-39
    ## [20,] 1.843731075704416e-29
    ## [21,] 8.689024055350509e-44
    ## [22,] 1.571325041547613e-39
    ## [23,] 1.039709965466667e-33
    ## [24,] 3.044443188856715e-38
    ## [25,] 5.744924574388899e-37
    ## [26,] 2.751884122100472e-33
    ## [27,] 2.118104573173345e-38
    ## [28,] 1.108007882961469e-30
    ## [29,] 7.421373217528215e-33
    ## [30,] 2.116291733758300e-43
    ## [31,] 3.827106675377243e-39
    ## [32,] 2.532309257526213e-33
    ## [33,] 1.107059563228856e-35
    ## [34,] 7.415021426378293e-38
    ## [35,] 1.340934122941629e-33
    ## [36,] 8.872655458246708e-28
    ## [37,] 1.399229224179350e-36
    ## [38,] 4.902591883813271e-31
    ## [39,] 5.522805772881025e-31
    ## [40,] 1.935069843462120e-25
    ## [41,] 1.151702307046043e-34
    ## [42,] 6.702466908554155e-33
    ## [43,] 3.779630191117644e-37
    ## [44,] 1.324299405600061e-31
    ## [45,] 9.213522648473533e-32
    ## [46,] 3.228215976137059e-26
    ## [47,] 3.689647247212229e-43
    ## [48,] 6.672366283128951e-39
    ## [49,] 4.414952689215923e-33
    ## [50,] 1.930102171080890e-35
    ## [51,] 1.292771358383097e-37
    ## [52,] 2.337850598044997e-33
    ## [53,] 1.546902455118760e-27
    ## [54,] 9.628731611879719e-31
    ## [55,] 7.035361581688835e-29
    ## [56,] 1.168541094037208e-32
    ## [57,] 6.589593442002521e-37
    ## [58,] 1.604956096484750e-36
    ## [59,] 1.920456010964904e-26
    ## [60,] 8.395732813361129e-29
    ## [61,] 5.623413659843774e-31
    ## [62,] 2.866400894039056e-30
    ## [63,] 9.999999999999982e-01
    ## [64,] 9.999999999997602e-01
    ## 
    ## $p.m
    ## [1] 0.0002857142857142756
    ## 
    ## $p.u
    ## [1] 0.9997142857142858
    ## 
    ## $p.gamma.k.m
    ## $p.gamma.k.m[[1]]
    ## [1] 3.340503443908832e-34 1.420086535562352e-26 1.000000000000000e+00
    ## 
    ## $p.gamma.k.m[[2]]
    ## [1] 1.340060306002473e-27 1.000000000000000e+00
    ## 
    ## $p.gamma.k.m[[3]]
    ## [1] 9.805415945148028e-33 7.449605841015272e-27 1.000000000000000e+00
    ## 
    ## $p.gamma.k.m[[4]]
    ## [1] 1.683733467905158e-26 1.000000000000000e+00
    ## 
    ## $p.gamma.k.m[[5]]
    ## [1] 7.865517568809706e-32 1.647762580084093e-26 1.000000000000000e+00
    ## 
    ## $p.gamma.k.m[[6]]
    ## [1] 2.377880188696157e-27 1.000000000000000e+00
    ## 
    ## $p.gamma.k.m[[7]]
    ## [1] 1.774578514503502e-26 1.000000000000000e+00
    ## 
    ## 
    ## $p.gamma.k.u
    ## $p.gamma.k.u[[1]]
    ## [1] 0.986041726207487867 0.011603315232923702 0.002354958559588454
    ## 
    ## $p.gamma.k.u[[2]]
    ## [1] 0.993302076356329500 0.006697923643670463
    ## 
    ## $p.gamma.k.u[[3]]
    ## [1] 0.9990225778793940803 0.0006516147470706017 0.0003258073735352958
    ## 
    ## $p.gamma.k.u[[4]]
    ## [1] 0.999833906071019518 0.000166093928980537
    ## 
    ## $p.gamma.k.u[[5]]
    ## [1] 0.998445270077164881 0.001166047442126322 0.000388682480708784
    ## 
    ## $p.gamma.k.u[[6]]
    ## [1] 0.8836924835667334 0.1163075164332666
    ## 
    ## $p.gamma.k.u[[7]]
    ## [1] 0.98807659331237496 0.01192340668762505
    ## 
    ## 
    ## $p.gamma.j.m
    ##                        [,1]
    ##  [1,] 7.259778058130710e-47
    ##  [2,] 1.544916893261878e-44
    ##  [3,] 2.074686714834686e-39
    ##  [4,] 2.560811692458025e-41
    ##  [5,] 5.449534699907367e-39
    ##  [6,] 7.318243002739926e-34
    ##  [7,] 2.560818952236063e-41
    ##  [8,] 5.449550149076333e-39
    ##  [9,] 7.318263749607016e-34
    ## [10,] 3.130767951581248e-43
    ## [11,] 1.104349174487971e-37
    ## [12,] 6.178644160313666e-38
    ## [13,] 1.314846220406412e-35
    ## [14,] 2.179452881696326e-32
    ## [15,] 2.179459060340499e-32
    ## [16,] 6.563159332676735e-42
    ## [17,] 2.315093198591687e-36
    ## [18,] 6.616032975475065e-29
    ## [19,] 2.315118806781198e-36
    ## [20,] 6.616106158112530e-29
    ## [21,] 1.514220338838310e-43
    ## [22,] 3.222336224138526e-41
    ## [23,] 4.327312481408308e-36
    ## [24,] 5.341270918352221e-38
    ## [25,] 6.530051567014713e-40
    ## [26,] 4.828797734565817e-33
    ## [27,] 1.230393957253638e-38
    ## [28,] 4.340087543772796e-33
    ## [29,] 4.340099847712333e-33
    ## [30,] 4.156313173645719e-41
    ## [31,] 8.844841239272702e-39
    ## [32,] 1.187783931548381e-33
    ## [33,] 1.466096523525085e-35
    ## [34,] 1.466100679838268e-35
    ## [35,] 3.119935195495635e-33
    ## [36,] 4.189796958963440e-28
    ## [37,] 1.792403566140903e-37
    ## [38,] 6.322536288955084e-32
    ## [39,] 3.537350579198452e-32
    ## [40,] 1.247767401617662e-26
    ## [41,] 3.757490294705181e-36
    ## [42,] 1.325434292636939e-30
    ## [43,] 8.669099649771470e-38
    ## [44,] 3.057943990050060e-32
    ## [45,] 7.044158337016904e-33
    ## [46,] 2.484761108059193e-27
    ## [47,] 6.643863436043056e-42
    ## [48,] 1.413847197074074e-39
    ## [49,] 1.898671708082991e-34
    ## [50,] 2.343554173954209e-36
    ## [51,] 2.343560817817626e-36
    ## [52,] 4.987214028916716e-34
    ## [53,] 6.697387241319161e-29
    ## [54,] 5.654452201201138e-33
    ## [55,] 2.118682269622400e-31
    ## [56,] 2.118705705230598e-31
    ## [57,] 1.385754917404610e-38
    ## [58,] 3.803694396993602e-36
    ## [59,] 1.087013152405612e-28
    ## [60,] 1.341713893780707e-30
    ## [61,] 1.341717697475093e-30
    ## [62,] 7.933619144462027e-33
    ## [63,] 4.999992912630213e-01
    ## [64,] 5.000007087369787e-01
    ## 
    ## $p.gamma.j.u
    ##                        [,1]
    ##  [1,] 2.857221956840588e-01
    ##  [2,] 3.362255594698992e-03
    ##  [3,] 6.823889000857922e-04
    ##  [4,] 1.926649999827365e-03
    ##  [5,] 2.267198642176678e-05
    ##  [6,] 4.601408620302853e-06
    ##  [7,] 2.876488456838862e-01
    ##  [8,] 3.384927581120759e-03
    ##  [9,] 6.869903087060948e-04
    ## [10,] 1.863623519939928e-04
    ## [11,] 1.876190097925842e-04
    ## [12,] 9.318147581124404e-05
    ## [13,] 1.096519427265986e-06
    ## [14,] 6.283309209696910e-07
    ## [15,] 9.380980673221370e-05
    ## [16,] 4.746449684758977e-05
    ## [17,] 4.778455414178877e-05
    ## [18,] 1.141236131965105e-07
    ## [19,] 2.876966302380280e-01
    ## [20,] 6.871044323192913e-04
    ## [21,] 3.336844252887254e-04
    ## [22,] 3.926654431255759e-06
    ## [23,] 7.969368547073538e-07
    ## [24,] 3.359344888386911e-04
    ## [25,] 2.176457246231368e-07
    ## [26,] 3.359902946589156e-04
    ## [27,] 1.112281342427360e-04
    ## [28,] 7.500211326127721e-07
    ## [29,] 1.119781553753490e-04
    ## [30,] 3.760543258673840e-02
    ## [31,] 4.425245151260446e-04
    ## [32,] 8.981286777065927e-05
    ## [33,] 2.535767531580321e-04
    ## [34,] 3.785900933989642e-02
    ## [35,] 4.455084970142911e-04
    ## [36,] 9.041848387010354e-05
    ## [37,] 2.452814996691950e-05
    ## [38,] 2.469354544841644e-05
    ## [39,] 1.226411444362296e-05
    ## [40,] 1.234681245045480e-05
    ## [41,] 6.247057328507991e-06
    ## [42,] 3.786529852168424e-02
    ## [43,] 4.391799919637748e-05
    ## [44,] 4.421414214369571e-05
    ## [45,] 1.463933207568891e-05
    ## [46,] 1.473804638478754e-05
    ## [47,] 3.447892452700791e-03
    ## [48,] 4.057331164370739e-05
    ## [49,] 8.234584410845844e-06
    ## [50,] 2.324944331852401e-05
    ## [51,] 3.471141896019317e-03
    ## [52,] 4.084690106746290e-05
    ## [53,] 8.290110940787820e-06
    ## [54,] 1.124447844914257e-06
    ## [55,] 5.766300485921191e-07
    ## [56,] 3.471718526067909e-03
    ## [57,] 4.026666562541004e-06
    ## [58,] 4.537956419029361e-04
    ## [59,] 1.083797876466995e-06
    ## [60,] 3.059984091542979e-06
    ## [61,] 4.568556259944790e-04
    ## [62,] 5.299712106872824e-07
    ## [63,] 1.541304070454287e-19
    ## [64,] 2.301166982865473e-17
    ## 
    ## $patterns.w
    ##    gamma.1 gamma.2 gamma.3 gamma.4 gamma.5 gamma.6 gamma.7 counts
    ## 1        0       0       0       0       0       0       0 100827
    ## 43       1       0       0       0       0       0       0    261
    ## 52       2       0       0       0       0       0       0   1193
    ## 19       0       2       0       0       0       0       0    690
    ## 48       1       2       0       0       0       0       0      1
    ## 57       2       2       0       0       0       0       0      9
    ## 25       0      NA       0       0       0       0       0  48376
    ## 49       1      NA       0       0       0       0       0    101
    ## 59       2      NA       0       0       0       0       0    563
    ## 13       0       0       1       0       0       0       0     37
    ## 39       0      NA       1       0       0       0       0     11
    ## 16       0       0       2       0       0       0       0     64
    ## 47       1       0       2       0       0       0       0      1
    ## 24       0       2       2       0       0       0       0      1
    ## 41       0      NA       2       0       0       0       0     36
    ## 11       0       0       0       2       0       0       0     15
    ## 33       0      NA       0       2       0       0       0      8
    ## 62       2      NA       0       2       0       0       0      1
    ## 35       0      NA       0      NA       0       0       0    322
    ## 63       2      NA       0      NA       0       0       0      4
    ## 5        0       0       0       0       1       0       0    124
    ## 46       1       0       0       0       1       0       0      2
    ## 56       2       0       0       0       1       0       0      4
    ## 29       0      NA       0       0       1       0       0     49
    ## 15       0       0       1       0       1       0       0      1
    ## 38       0      NA       0      NA       1       0       0      1
    ## 9        0       0       0       0       2       0       0     33
    ## 23       0       2       0       0       2       0       0      1
    ## 31       0      NA       0       0       2       0       0     28
    ## 3        0       0       0       0       0       2       0  12998
    ## 45       1       0       0       0       0       2       0     30
    ## 54       2       0       0       0       0       2       0    154
    ## 21       0       2       0       0       0       2       0     75
    ## 27       0      NA       0       0       0       2       0   6689
    ## 51       1      NA       0       0       0       2       0      9
    ## 61       2      NA       0       0       0       2       0     75
    ## 14       0       0       1       0       0       2       0      7
    ## 40       0      NA       1       0       0       2       0      1
    ## 18       0       0       2       0       0       2       0      8
    ## 42       0      NA       2       0       0       2       0      3
    ## 12       0       0       0       2       0       2       0      4
    ## 37       0      NA       0      NA       0       2       0     20
    ## 7        0       0       0       0       1       2       0     13
    ## 30       0      NA       0       0       1       2       0      8
    ## 10       0       0       0       0       2       2       0      2
    ## 32       0      NA       0       0       2       2       0      4
    ## 2        0       0       0       0       0       0       2   1199
    ## 44       1       0       0       0       0       0       2      6
    ## 53       2       0       0       0       0       0       2     19
    ## 20       0       2       0       0       0       0       2     10
    ## 26       0      NA       0       0       0       0       2    592
    ## 50       1      NA       0       0       0       0       2      1
    ## 60       2      NA       0       0       0       0       2      5
    ## 17       0       0       2       0       0       0       2      1
    ## 34       0      NA       0       2       0       0       2      1
    ## 36       0      NA       0      NA       0       0       2      3
    ## 6        0       0       0       0       1       0       2      1
    ## 4        0       0       0       0       0       2       2    149
    ## 55       2       0       0       0       0       2       2      3
    ## 22       0       2       0       0       0       2       2      3
    ## 28       0      NA       0       0       0       2       2     92
    ## 8        0       0       0       0       1       2       2      1
    ## 58       2       2       2       2       2       2       2     43
    ## 64       2      NA       2       2       2       2       2      7
    ##                weights           p.gamma.j.m           p.gamma.j.u
    ## 1  -104.98641482887176 7.259778058130710e-47 2.857221956840588e-01
    ## 43  -95.18363075131680 1.544916893261878e-44 3.362255594698992e-03
    ## 52  -81.78109763737348 2.074686714834686e-39 6.823889000857922e-04
    ## 19  -87.21369200154120 2.560811692458025e-41 1.926649999827365e-03
    ## 48  -77.41090792398624 5.449534699907367e-39 2.267198642176678e-05
    ## 57  -64.00837481004292 7.318243002739926e-34 4.601408620302853e-06
    ## 25  -92.21964687360556 2.560818952236063e-41 2.876488456838862e-01
    ## 49  -82.41686279605058 5.449550149076333e-39 3.384927581120759e-03
    ## 59  -69.01432968210726 7.318263749607016e-34 6.869903087060948e-04
    ## 13  -89.28206302194938 3.130767951581248e-43 1.863623519939928e-04
    ## 39  -76.51529506668318 1.104349174487971e-37 1.876190097925842e-04
    ## 16  -76.39617306449129 6.178644160313666e-38 9.318147581124404e-05
    ## 47  -66.59338898693633 1.314846220406412e-35 1.096519427265986e-06
    ## 24  -58.62345023716075 2.179452881696326e-32 6.283309209696910e-07
    ## 41  -63.62940510922508 2.179459060340499e-32 9.380980673221370e-05
    ## 11  -84.87157325195780 6.563159332676735e-42 4.746449684758977e-05
    ## 33  -72.10480529669159 2.315093198591687e-36 4.778455414178877e-05
    ## 62  -48.89948810519331 6.616032975475065e-29 1.141236131965105e-07
    ## 35  -80.80775361988402 2.315118806781198e-36 2.876966302380280e-01
    ## 63  -57.60243642838574 6.616106158112530e-29 6.871044323192913e-04
    ## 5   -90.59094347398857 1.514220338838310e-43 3.336844252887254e-04
    ## 46  -80.78815939643360 3.222336224138526e-41 3.926654431255759e-06
    ## 56  -67.38562628249028 4.327312481408308e-36 7.969368547073538e-07
    ## 29  -77.82417551872236 5.341270918352221e-38 3.359344888386911e-04
    ## 15  -74.88659166706618 6.530051567014713e-40 2.176457246231368e-07
    ## 38  -66.41228226500083 4.828797734565817e-33 3.359902946589156e-04
    ## 9   -78.18697192285298 1.230393957253638e-38 1.112281342427360e-04
    ## 23  -60.41424909552242 4.340087543772796e-33 7.500211326127721e-07
    ## 31  -65.42020396758677 4.340099847712333e-33 1.119781553753490e-04
    ## 3   -89.70075363244200 4.156313173645719e-41 3.760543258673840e-02
    ## 45  -79.89796955488704 8.844841239272702e-39 4.425245151260446e-04
    ## 54  -66.49543644094371 1.187783931548381e-33 8.981286777065927e-05
    ## 21  -71.92803080511146 1.466096523525085e-35 2.535767531580321e-04
    ## 27  -76.93398567717580 1.466100679838268e-35 3.785900933989642e-02
    ## 51  -67.13120159962082 3.119935195495635e-33 4.455084970142911e-04
    ## 61  -53.72866848567752 4.189796958963440e-28 9.041848387010354e-05
    ## 14  -73.99640182551963 1.792403566140903e-37 2.452814996691950e-05
    ## 40  -61.22963387025342 6.322536288955084e-32 2.469354544841644e-05
    ## 18  -61.11051186806154 3.537350579198452e-32 1.226411444362296e-05
    ## 42  -48.34374391279533 1.247767401617662e-26 1.234681245045480e-05
    ## 12  -69.58591205552804 3.757490294705181e-36 6.247057328507991e-06
    ## 37  -65.52209242345427 1.325434292636939e-30 3.786529852168424e-02
    ## 7   -75.30528227755882 8.669099649771470e-38 4.391799919637748e-05
    ## 30  -62.53851432229259 3.057943990050060e-32 4.421414214369571e-05
    ## 10  -62.90131072642322 7.044158337016904e-33 1.463933207568891e-05
    ## 32  -50.13454277115701 2.484761108059193e-27 1.473804638478754e-05
    ## 2   -89.14488815080966 6.643863436043056e-42 3.447892452700791e-03
    ## 44  -79.34210407325470 1.413847197074074e-39 4.057331164370739e-05
    ## 53  -65.93957095931138 1.898671708082991e-34 8.234584410845844e-06
    ## 20  -71.37216532347912 2.343554173954209e-36 2.324944331852401e-05
    ## 26  -76.37812019554346 2.343560817817626e-36 3.471141896019317e-03
    ## 50  -66.57533611798850 4.987214028916716e-34 4.084690106746290e-05
    ## 60  -53.17280300404518 6.697387241319161e-29 8.290110940787820e-06
    ## 17  -60.55464638642921 5.654452201201138e-33 1.124447844914257e-06
    ## 34  -56.26327861862950 2.118682269622400e-31 5.766300485921191e-07
    ## 36  -64.96622694182192 2.118705705230598e-31 3.471718526067909e-03
    ## 6   -74.74941679592646 1.385754917404610e-38 4.026666562541004e-06
    ## 4   -73.85922695437991 3.803694396993602e-36 4.537956419029361e-04
    ## 55  -50.65390976288163 1.087013152405612e-28 1.083797876466995e-06
    ## 22  -56.08650412704936 1.341713893780707e-30 3.059984091542979e-06
    ## 28  -61.09245899911370 1.341717697475093e-30 4.568556259944790e-04
    ## 8   -59.46375559949671 7.933619144462027e-33 5.299712106872824e-07
    ## 58   42.62333931176212 4.999992912630213e-01 1.541304070454287e-19
    ## 64   37.61738443969779 5.000007087369787e-01 2.301166982865473e-17
    ## 
    ## $iter.converge
    ## [1] 3
    ## 
    ## $nobs.a
    ## [1] 500
    ## 
    ## $nobs.b
    ## [1] 350
    ## 
    ## $varnames
    ## [1] "firstname"  "middlename" "lastname"   "housenum"   "streetname"
    ## [6] "city"       "birthyear" 
    ## 
    ## attr(,"class")
    ## [1] "fastLink"    "fastLink.EM"

which is a list of parameter estimates for different fields. These fields are:

-   `zeta.j`: The posterior match probabilities for each unique pattern.

-   `p.m`: The posterior probability of a pair matching.

-   `p.u`: The posterior probability of a pair not matching.

-   `p.gamma.k.m`: The posterior of the matching probability for a specific matching field.

-   `p.gamma.k.u`: The posterior of the non-matching probability for a specific matching field.

-   `p.gamma.j.m`: The posterior probability that a pair is in the matched set given a particular agreement pattern.

-   `p.gamma.j.u`: The posterior probability that a pair is in the unmatched set given a particular agreement pattern.

-   `patterns.w`: Counts of the agreement patterns observed (2 = match, 1 = partial match, 0 = non-match), along with the Felligi-Sunter Weights.

-   `iter.converge`: The number of iterations it took the EM algorithm to converge.

-   `nobs.a`: The number of observations in dataset A.

-   `nobs.b`: The number of observations in dataset B.

Lastly, we can summarize the accuracy of the match using the `summary()` function:

``` r
summary(matches.out)
```

    ##                   95%     85%     75%   Exact
    ## 1 Match Count      50      50      50      43
    ## 2  Match Rate 14.286% 14.286% 14.286% 12.286%
    ## 3         FDR      0%      0%      0%        
    ## 4         FNR      0%      0%      0%

where each column gives the match count, match rate, false discovery rate (FDR) and false negative rate (FNR) under different cutoffs for matches based on the posterior probability of a match. Other arguments include:

-   `num.comparisons`: The number of comparisons attempted for each observation in an across-state merge. For instance, if matching each state's voter file to every other state's voter file to try and find movers, `num.comparisons` = 49. Default is 1.

-   `thresholds`: A vector of thresholds between 0 and 1 to summarize the match.

-   `weighted`: Whether to weight the FDR and FNR calculations when doing across-state matches, so that the pooled FDR and FNR calculations are the sum of the within and across-geography FDR and FNR. Default is TRUE.

-   `digits`: Number of digits to include in the summary object. Default is 3.

### Preprocessing Matches via Blocking

In order to reduce the number of pairwise comparisons that need to be conducted, researchers will often block similar observations from dataset A and dataset B together so that comparisons are only made between these maximally similar groups. Here, we implement a form of this clustering that uses word embedding, a common preprocessing method for textual data, to form maximally similar groups.

In , the function `blockData()` can block two data sets using a single variable or combinations of variables using several different blocking techniques. The basic functionality is similar to that of `fastLink()`, where the analyst inputs two data sets and a vector of variable names that they want to block on. A simple example follows, where we are blocking the two sample data sets by gender:

    ## 
    ## ==================== 
    ## fastLink(): Fast Probabilistic Record Linkage
    ## ==================== 
    ## 
    ## Calculating matches for each variable.
    ## Getting counts for parameter estimation.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Running the EM algorithm.
    ## Getting the indices of estimated matches.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Deduping the estimated matches.
    ## Getting the match patterns for each estimated match.

``` r
blockgender_out <- blockData(dfA, dfB, varnames = "gender")
```

    ## 
    ## ==================== 
    ## blockData(): Blocking Methods for Record Linkage
    ## ==================== 
    ## 
    ## Blocking variables.
    ##     Blocking variable gender using exact blocking.
    ## 
    ## Combining blocked variables for final blocking assignments.

``` r
names(blockgender_out)
```

    ## [1] "block.1" "block.2"

In its simplest usage, takes two data sets and a single variable name for the argument, and it returns the indices of the member observations for each block. Data sets can then be subsetted as follows and the match can then be run within each block separately:

``` r
## Subset dfA into blocks
dfA_block1 <- dfA[blockgender_out$block.1$dfA.inds,]
dfA_block2 <- dfA[blockgender_out$block.2$dfA.inds,]

## Subset dfB into blocks
dfB_block1 <- dfB[blockgender_out$block.1$dfB.inds,]
dfB_block2 <- dfB[blockgender_out$block.2$dfB.inds,]

## Run fastLink on each
fl_out_block1 <- fastLink(
  dfA_block1, dfB_block1,
  varnames = c("firstname", "lastname", "housenum",
               "streetname", "city", "birthyear")
)
```

    ## 
    ## ==================== 
    ## fastLink(): Fast Probabilistic Record Linkage
    ## ==================== 
    ## 
    ## Calculating matches for each variable.
    ## Getting counts for parameter estimation.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Running the EM algorithm.
    ## Getting the indices of estimated matches.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Deduping the estimated matches.
    ## Getting the match patterns for each estimated match.

``` r
fl_out_block2 <- fastLink(
  dfA_block2, dfB_block2,
  varnames = c("firstname", "lastname", "housenum",
               "streetname", "city", "birthyear")
)
```

    ## 
    ## ==================== 
    ## fastLink(): Fast Probabilistic Record Linkage
    ## ==================== 
    ## 
    ## Calculating matches for each variable.
    ## Getting counts for parameter estimation.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Running the EM algorithm.
    ## Getting the indices of estimated matches.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Deduping the estimated matches.
    ## Getting the match patterns for each estimated match.

`blockData()` also implements other methods of blocking other than exact blocking. Analysts commonly use *window blocking* for numeric variables, where a given observation in dataset A will be compared to all observations in dataset B where the value of the blocking variable is within ±*K* of the value of the same variable in dataset A. The value of *K* is the size of the window --- for instance, if we wanted to compare observations where birth year is within ±1 year, the window size is 1. Below, we block `dfA` and `dfB` on gender and birth year, using exact blocking on gender and window blocking with a window size of 1 on birth year:

``` r
## Exact block on gender, window block (+/- 1 year) on birth year
blockdata_out <- blockData(dfA, dfB, varnames = c("gender", "birthyear"),
                           window.block = "birthyear", window.size = 1)
```

    ## 
    ## ==================== 
    ## blockData(): Blocking Methods for Record Linkage
    ## ==================== 
    ## 
    ## Blocking variables.
    ##     Blocking variable gender using exact blocking.
    ##     Blocking variable birthyear using window blocking.
    ## 
    ## Combining blocked variables for final blocking assignments.

`blockData()` also allows users to block variables using k-means clustering, so that similar values of string and numeric variables are blocked together. When applying k-means blocking to string variables such as name, the algorithm orders observations so that alphabetically close names are grouped together in a block. In the following example, we block `dfA` and `dfB` on gender and first name, again using exact blocking on gender and k-means blocking on first name while specifying 2 clusters for the k-means algorithm:

``` r
## Exact block on gender, k-means block on first name with 2 clusters
blockdata_out <- blockData(dfA, dfB, varnames = c("gender", "firstname"),
                           kmeans.block = "firstname", nclusters = 2)
```

    ## 
    ## ==================== 
    ## blockData(): Blocking Methods for Record Linkage
    ## ==================== 
    ## 
    ## Blocking variables.
    ##     Blocking variable gender using exact blocking.
    ##     Blocking variable firstname using k-means blocking.
    ## 
    ## Combining blocked variables for final blocking assignments.

Using Auxiliary Information to Inform `fastLink`
------------------------------------------------

The `fastLink` algorithm also includes several ways to incorporate auxiliary information on migration behavior to inform the matching of data sets over time. Auxiliary information is incorporated into the estimation as priors on two parameters of the model:

-   
    *λ*
    : The probability that a randomly selected pair of observations from dataset A and dataset B are a true match. When matching, for example, the same state to itself in subsequent years, the prior for this quantity is equal to the number of non-movers to the number of in-state movers, divided by the size of the cross-product of A and B. When matching two different states in subsequent years to find movers, the numerator is the size of the outflow from state A to state B, divided by the size of the cross-product of A and B.

-   
    *π*<sub>*k*, *l*</sub>
    : The probability that an address field does not match conditional on being in the matched set. Specified when trying to find movers within the same geography over time. For example, when trying to find movers within the same state over time, this quantity is equal to the estimated number of in-state movers divided by the number of in-state movers and non-movers.

The functions `calcMoversPriors()` can be used to calculate estimates for the corresponding prior distributions using the IRS Statistics of Income Migration Data.

Below, we show an example where we incorporate the auxiliary moving information for California into our estimates. First, we use `calcMoversPriors()` to estimate optimal parameter values for the priors:

``` r
priors.out <- calcMoversPriors(geo.a = "CA", geo.b = "CA", year.start = 2014, year.end = 2015)
names(priors.out)
```

    ## [1] "lambda.prior" "pi.prior"

where the `lambda.prior` entry is the estimate of the match rate, while `pi.prior` is the estimate of the in-state movers rate.

The `calcMoversPriors()` function accepts the following functions:

-   `geo.a`: The state name or county name of dataset A

-   `geo.b`: The state name or county name of dataset B

-   `year.start`: The year of dataset A

-   `year.end`: The year of dataset B

-   `county`: Boolean, whether the geographies in `geo.a` or `geo.b` refer to counties or states. Default is FALSE

-   `state.a`: If `county = TRUE`, the name of the state for `geo.a`

-   `state.b`: If `county = TRUE`, the name of the state for `geo.b`

-   `matchrate.lambda`: If TRUE, then returns the match rate for lambda (the expected share of observations in dataset A that can be found in dataset B). If FALSE, then returns the expected share of matches across all pairwise comparisons of datasets A and B. Default is FALSE.

-   `remove.instate`: If TRUE, then for calculating cross-state movers rates assumes that successful matches have been subsetted out. The interpretation of the prior is then the match rate conditional on being an out-of-state or county mover. Default is TRUE.

### Incorporating Auxiliary Information with `fastLink()` Wrapper

We can re-run the full match above while incorporating auxiliary information as follows:

``` r
## Reasonable prior estimates for this dataset
priors.out <- list(lambda.prior = 50/(nrow(dfA) * nrow(dfB)), pi.prior = 0.02)

matches.out.aux <- fastLink(
  dfA = dfA, dfB = dfB, 
  varnames = c("firstname", "middlename", "lastname", "housenum", "streetname", "city", "birthyear"),
  stringdist.match = c("firstname", "middlename", "lastname", "streetname", "city"),
  partial.match = c("firstname", "lastname", "streetname"),
  priors.obj = priors.out, 
  w.lambda = .5, w.pi = .5, 
  address.field = "streetname"
)
```

    ## 
    ## ==================== 
    ## fastLink(): Fast Probabilistic Record Linkage
    ## ==================== 
    ## 
    ## Calculating matches for each variable.
    ## Getting counts for parameter estimation.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Running the EM algorithm.
    ## Getting the indices of estimated matches.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Deduping the estimated matches.
    ## Getting the match patterns for each estimated match.

where `priors.obj` is an input for the the optimal prior parameters. This can be calculated by `calcMoversPriors()`, or can be provided by the user as a list with two entries named `lambda.prior` and `pi.prior`. `w.lambda` and `w.pi` are user-specified weights between 0 and 1 indicating the weighting between the MLE estimate and the prior, where a weight of 0 indicates no weight being placed on the prior. `address_field` is a vector of booleans of the same length as `varnames`, where `TRUE` indicates an address-related field used for matching.

Aggregating Multiple Matches Together
-------------------------------------

Often, we run several different matches for a single data set - for instance, when blocking by gender or by some other criterion to reduce the number of pairwise comparisons. Here, we walk through how to aggregate those multiple matches into a single summary. Here, we run `fastLink()` on the subsets of data defined by blocking on gender in the previous section:

``` r
## Run fastLink on each
link.1 <- fastLink(
  dfA_block1, dfB_block1,
  varnames = c("firstname", "lastname", "housenum",
               "streetname", "city", "birthyear")
)
```

    ## 
    ## ==================== 
    ## fastLink(): Fast Probabilistic Record Linkage
    ## ==================== 
    ## 
    ## Calculating matches for each variable.
    ## Getting counts for parameter estimation.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Running the EM algorithm.
    ## Getting the indices of estimated matches.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Deduping the estimated matches.
    ## Getting the match patterns for each estimated match.

``` r
link.2 <- fastLink(
  dfA_block2, dfB_block2,
  varnames = c("firstname", "lastname", "housenum",
               "streetname", "city", "birthyear")
)
```

    ## 
    ## ==================== 
    ## fastLink(): Fast Probabilistic Record Linkage
    ## ==================== 
    ## 
    ## Calculating matches for each variable.
    ## Getting counts for parameter estimation.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Running the EM algorithm.
    ## Getting the indices of estimated matches.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Deduping the estimated matches.
    ## Getting the match patterns for each estimated match.

To aggregate the two matches into a single summary, we use the `aggregateEM()` function as follows:

``` r
agg.out <- aggregateEM(em.list = list(link.1, link.2))
```

`aggregateEM()` accepts two arguments:

-   `em.list`: A list of either `fastLink` or `fastLink.EM` objects to be aggregated together.

-   `within.geo`: A vector of booleans the same length of `em.list`, to be used if the user wants to aggregate together within-geography matches (for instance, CA 2015 voter file to CA 2016 voter file) and across-geography matches (for instance, CA 2015 voter file to NY 2016 voter file). For entry `i` in `em.list`, `within.geo = TRUE` if it is a within-geography match, and `FALSE` if an across-geogpraphy match. Default is `NULL` (assumes all matches are within-geography).

We can then summarize the aggregated output as done previously:

``` r
summary(agg.out)
```

    ##                   95%     85%     75%   Exact
    ## 1 Match Count      50      50      50      50
    ## 2  Match Rate 14.286% 14.286% 14.286% 14.286%
    ## 3         FDR      0%      0%      0%        
    ## 4         FNR      0%      0%      0%

If we assume that the first `fastLink` run was for a within-geography match and the second was an across-geography match, the call to `aggregateEM()` would be:

``` r
agg.out <- aggregateEM(em.list = list(link.1, link.2), within.geo = c(TRUE, FALSE))
summary(agg.out)
```

    ##                                 95%     85%     75%   Exact
    ## 1  Match Count          All      50      50      50      50
    ## 2              Within-State      24      24      24      24
    ## 3              Across-State      26      26      26      26
    ## 4   Match Rate          All 27.322% 27.322% 27.322% 27.322%
    ## 5              Within-State 13.115% 13.115% 13.115% 13.115%
    ## 6              Across-State 14.208% 14.208% 14.208% 14.208%
    ## 7          FDR          All      0%      0%      0%        
    ## 8              Within-State      0%      0%      0%        
    ## 9              Across-State      0%      0%      0%        
    ## 10         FNR          All      0%      0%      0%        
    ## 11             Within-State      0%      0%      0%        
    ## 12             Across-State      0%      0%      0%

Random Sampling with `fastLink`
-------------------------------

The probabilistic modeling framework of `fastLink` is especially flexible in that it allows us to run the matching algorithm on a random smaller subset of data to be matched, and then apply those estimates to the full sample of data. This may be desired, for example, when using blocking along with a prior. We may want to block in order to reduce the number of pairwise comparisons, but may also be uncomfortable making the assumption that the same prior applies to all blocks uniformly. Random sampling allows us to run the EM algorithm with priors on a random sample from the full dataset, and the estimates can then be applied to each block separately to get matches for the entire dataset.

This functionality is incorporated into the `fastLink()` wrapper, which we show in the following example:

``` r
## Take 30% random samples of dfA and dfB
dfA.s <- dfA[sample(1:nrow(dfA), nrow(dfA) * .3),]
dfB.s <- dfB[sample(1:nrow(dfB), nrow(dfB) * .3),]

## Run the algorithm on the random samples
rs.out <- fastLink(
  dfA = dfA.s, dfB = dfB.s, 
  varnames = c("firstname", "middlename", "lastname", "housenum", "streetname", "city", "birthyear"),
  stringdist.match = c("firstname", "middlename", "lastname", "streetname", "city"),
  partial.match = c("firstname", "lastname", "streetname"),
  estimate.only = TRUE
)
```

    ## 
    ## ==================== 
    ## fastLink(): Fast Probabilistic Record Linkage
    ## ==================== 
    ## 
    ## Calculating matches for each variable.
    ## Getting counts for parameter estimation.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Running the EM algorithm.

``` r
class(rs.out)
```

    ## [1] "fastLink"    "fastLink.EM"

``` r
## Apply to the whole dataset
fs.out <- fastLink(
  dfA = dfA, dfB = dfB, 
  varnames = c("firstname", "middlename", "lastname", "housenum", "streetname", "city", "birthyear"),
  stringdist.match = c("firstname", "middlename", "lastname", "streetname", "city"),
  partial.match = c("firstname", "lastname", "streetname"),
  em.obj = rs.out
)
```

    ## 
    ## ==================== 
    ## fastLink(): Fast Probabilistic Record Linkage
    ## ==================== 
    ## 
    ## Calculating matches for each variable.
    ## Getting counts for parameter estimation.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Imputing matching probabilities using provided EM object.
    ## Getting the indices of estimated matches.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Deduping the estimated matches.
    ## Getting the match patterns for each estimated match.

``` r
summary(fs.out)
```

    ##                   95%     85%     75%   Exact
    ## 1 Match Count      50      50      50      43
    ## 2  Match Rate 14.286% 14.286% 14.286% 12.286%
    ## 3         FDR      0%      0%      0%        
    ## 4         FNR      0%      0%      0%

In the first run of `fastLink()`, we specify `estimate.only = TRUE`, which runs the algorithm only through the EM estimation step and returns the EM object. In the second run of `fastLink()`, we provide the EM object from the first stage as an argument to `em.obj`. Then, using the parameter values calculated in the previous EM stage, we estimate posterior probabilities of belonging to the matched set for all matching patterns in the full dataset that were not present in the random sample.

Finding Duplicates within a Dataset via `fastLink`
--------------------------------------------------

The following lines of code represent an example on how to find duplicates withing a dataset via `fastLink`. As before, we use `fastLink()` (the wrapper function) to do the merge. `fastLink()` will automatically detect that two datasets are identical, and will use the probabilistic match algorithm to indicate duplicated entries in the `dedupe.ids` covariate in the returned data frame.

``` r
## Add duplicates
dfA <- rbind(dfA, dfA[sample(1:nrow(dfA), 10, replace = FALSE),])

## Run fastLink
fl_out_dedupe <- fastLink(
  dfA = dfA, dfB = dfA,
  varnames = c("firstname", "lastname", "housenum",
               "streetname", "city", "birthyear")
)
```

    ## 
    ## ==================== 
    ## fastLink(): Fast Probabilistic Record Linkage
    ## ==================== 
    ## 
    ## dfA and dfB are identical, assuming deduplication of a single data set.
    ## Setting return.all to FALSE.
    ## 
    ## Calculating matches for each variable.
    ## Getting counts for parameter estimation.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Running the EM algorithm.
    ## Getting the indices of estimated matches.
    ##     Parallelizing calculation using OpenMP. 1 threads out of 8 are used.
    ## Calculating the posterior for each pair of matched observations.
    ## Getting the match patterns for each estimated match.

``` r
## Run getMatches
dfA_dedupe <- getMatches(dfA = dfA, dfB = dfA, fl.out = fl_out_dedupe)

## Look at the IDs of the duplicates
names(table(dfA_dedupe$dedupe.ids)[table(dfA_dedupe$dedupe.ids) > 1])
```

    ##  [1] "501" "502" "503" "504" "505" "506" "507" "508" "509" "510"

``` r
## Show duplicated observation
dfA_dedupe[dfA_dedupe$dedupe.ids == 501,]
```

    ##      firstname middlename lastname housenum streetname      city birthyear
    ## 297    anthony          r labrecue     2040    iris ct Livermore      1944
    ## 2971   anthony          r labrecue     2040    iris ct Livermore      1944
    ##      gender dedupe.ids
    ## 297       F        501
    ## 2971      F        501
