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

joint_dist <- createDist(formula_u, target_data = dfOnoRep, exp_data  = dfOnoRep, type = "joint")

exp_marginal <- createDist(formula_u, target_data = dfOnoRep, exp_data  = dfOnoRep, type = "marginal")

load("example/ono_data_joint_distribution.RData")
full_joint <- dfJoint

is.element(all.vars(formula_u), colnames(full_joint))

#
marginal_rep <- tapply(distList$`Rep Pop`$prop, distList$`Rep Pop`$factor, function(x) x)
marginal_rep_n <- tapply(distList$`Rep Pop`$levels, distList$`Rep Pop`$factor, function(x) x)
for(i in 1:length(marginal_rep)){
  names(marginal_rep[[i]])  <- marginal_rep_n[[i]]
}

marginal_dem <- tapply(distList$`Dem Pop`$prop, distList$`Dem Pop`$factor, function(x) x)
marginal_dem_n <- tapply(distList$`Dem Pop`$levels, distList$`Dem Pop`$factor, function(x) x)
for(i in 1:length(marginal_dem)){
  names(marginal_dem[[i]])  <- marginal_dem_n[[i]]
}

# Experimental Marginal Distribution
exp_marginal <- createDist(formula_u, target_data = dfOnoRep, exp_data  = dfOnoRep, type = "marginal")

# Experimental Joint Distribution
exp_joint <- createDist(formula_u, target_data = dfOnoRep, exp_data  = dfOnoRep, type = "joint")

# Target Data
target_data <- full_joint[, is.element(colnames(full_joint), all.vars(formula_u)) == TRUE]

## Need to have the same level names
levels(target_data$age)  <-  levels(dfOnoRep$age)
levels(target_data$race) <-  levels(dfOnoRep$race)
levels(target_data$experience)   <- levels(dfOnoRep$experience)
levels(target_data$pos_security) <- levels(dfOnoRep$pos_security)
levels(target_data$pos_deficit)  <- levels(dfOnoRep$pos_deficit)

target_dist2 <- list(exp_marginal, exp_joint, marginal_rep, marginal_dem,  target_data)
names(target_dist2) <- c("Exp_Mar", "Exp-Joint", "Rep_Mar", "Dem_Mar", "Full-Joint")

ameOut <- AME_estimate_full(formula = as.formula(formula_u),
                            data = dfOnoRep, type = "No-Reg",
                            pair = TRUE, pair_id = dfOnoRep$pair_id,
                            cluster = dfOnoRep$id,
                            target_dist = target_dist2,
                            target_type = c("marginal", "joint", "marginal", "marginal" , "target_data"),
                            boot = 100)

# Figure 1: Estimates of AMCE
plot_AME(ameOut, factor_name = c("gender", "race"),
         plot_difference = "none",
         plot_type  = c("Exp_Mar", "Exp-Joint", "Rep_Mar", "Dem_Mar", "Full-Joint"),
         col = c("black",  "gray", "red", "blue", "green"))

# Figure 2: Decompose Bias
plot_decompose(ameOut, factor_name = "gender", level_name = "Male",
               marginal_diff = c("Exp", "Exp-Joint"))

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
