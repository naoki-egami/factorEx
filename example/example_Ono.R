rm(list=ls())

## Please install "genlasso" by yourself
library("genlasso")
library(PopCon)
library(arm)
library(sandwich)
library("igraph")
library(prodlim)

# load data
load("example/ono_data_cleaned.RData")

formula_u <- Y ~ age + family + race + experience + trait + party +
  policy_expertise + pos_security + pos_immigrants + pos_abortion +
  pos_deficit + fav_rating + gender

marginal_dist <- list(distList$Exp, distList$`Rep Pop`)

ameOut <- AME_estimate_full(formula = as.formula(formula_u),
                            data = dfOnoRep,
                            pair = TRUE, pair_id = dfOnoRep$pair_id,
                            cluster = dfOnoRep$id,
                            marginal_dist = marginal_dist, marginal_type = c("Exp", "Rep Pop"),
                            boot = 100)

# Figure 1: Estimates of AMCE
plot_AME(ameOut, factor_name = c("gender", "race"),
         plot_difference = "none",
         plot_type = c("Exp", "Rep Pop"), plot_name = c("Exp", "Rep Pop"),
         col = c("black", "red"))

# Figure 2: Decompose Bias
plot_decompose(ameOut, factor_name = "gender", level_name = "Male",
               marginal_diff = c("Rep Pop", "Exp"))

plot_decompose(ameOut, factor_name = "race", level_name = "Black",
               marginal_diff = c("Rep Pop", "Exp"))

# Note: the last row is the across-profile interaction with the same factor

# Note: trait has zero effect because "Exp" and "Rep Pop" use the same marginals
marginal_dist[[1]][marginal_dist[[1]]$fac == "trait", ]
marginal_dist[[2]][marginal_dist[[1]]$fac == "trait", ]


# Figure 3: Conditional Effect plots
cAME_gender <- cAME_from_boot(ameOut, factor_name = "gender", level_name = "Male", difference = TRUE)
plot_cAME(cAME_gender, factor_name = c(setdiff(all.vars(formula_u), c("Y", "gender")), "gender"),
          marginal_prop = c("Exp", "Rep Pop"), marginal_effect = "Exp", col = c("black", "red"))
## the last row is the across-profile interaction with the same factor

cAME_raceB <- cAME_from_boot(ameOut, factor_name = "race", level_name = "Black", difference = TRUE)
plot_cAME(cAME_raceB, factor_name = c(setdiff(all.vars(formula_u), c("Y", "race")), "race"),
          marginal_prop = c("Exp", "Rep Pop"), marginal_effect = "Exp", col = c("black", "red"))

## Note: if there are too many, focus on a few factors that matter
