rm(list=ls())

## Please install "genlasso" by yourself
library(PopCon)
library(arm)
library("genlasso")
library(sandwich)
library("igraph")
library(prodlim)

# load data
load("example/ono_data_cleaned.RData")
load("example/ono_data_joint_distribution.RData")
full_joint <- dfJoint

# I changed the format of marginal_distritbution. Please check this
## Republican
marginal_rep <- tapply(distList$`Rep Pop`$prop, distList$`Rep Pop`$factor, function(x) x)
marginal_rep_n <- tapply(distList$`Rep Pop`$levels, distList$`Rep Pop`$factor, function(x) x)
for(i in 1:length(marginal_rep)){
  names(marginal_rep[[i]])  <- marginal_rep_n[[i]]
}
marginal_rep

## Democrat
marginal_dem <- tapply(distList$`Dem Pop`$prop, distList$`Dem Pop`$factor, function(x) x)
marginal_dem_n <- tapply(distList$`Dem Pop`$levels, distList$`Dem Pop`$factor, function(x) x)
for(i in 1:length(marginal_dem)){
  names(marginal_dem[[i]])  <- marginal_dem_n[[i]]
}

# formula
formula_u <- Y ~ age + family + race + experience + trait + party +
  policy_expertise + pos_security + pos_immigrants + pos_abortion +
  pos_deficit + fav_rating + gender

# Create  Experimental Marginal Distribution
exp_marginal <- createDist(formula_u, target_data = dfOnoRep, exp_data  = dfOnoRep, type = "marginal")

# Create Experimental Joint Distribution
exp_joint <- createDist(formula_u, target_data = dfOnoRep, exp_data  = dfOnoRep, type = "joint")

# Target Data (columns should be restricted to what in 'data')
target_data <- full_joint[, is.element(colnames(full_joint), all.vars(formula_u)) == TRUE]

## Need to have the same level names
levels(target_data$age)  <-  levels(dfOnoRep$age)
levels(target_data$race) <-  levels(dfOnoRep$race)
levels(target_data$experience)   <- levels(dfOnoRep$experience)
levels(target_data$pos_security) <- levels(dfOnoRep$pos_security)
levels(target_data$pos_deficit)  <- levels(dfOnoRep$pos_deficit)

target_dist <- list(exp_marginal, exp_joint, marginal_rep, marginal_dem,  target_data)
names(target_dist) <- c("Exp_Mar", "Exp-Joint", "Rep_Mar", "Dem_Mar", "Full-Joint")

## Need to specify "target_type"; one of "marginal", "joint", and "target_data"
## For us, we will only compare "marginal" and "target"
## But for illustration, I show everything here

# Without regularization
ameOut <- AME_estimate_full(formula = as.formula(formula_u),
                            formula_three = ~  age*family*race + family*experience*party, ## Please insert plausible three-way interactions
                            data = dfOnoRep, type = "No-Reg",
                            pair = TRUE, pair_id = dfOnoRep$pair_id,
                            cluster = dfOnoRep$id,
                            target_dist = target_dist,
                            target_type = c("marginal", "joint", "marginal", "marginal" , "target_data"),
                            boot = 100)

# With regularization
ameOut_reg <- AME_estimate_full(formula = as.formula(formula_u),
                                formula_three = ~  age*family*race + family*experience*party,
                                data = dfOnoRep, type = "genlasso",
                                pair = TRUE, pair_id = dfOnoRep$pair_id,
                                cluster = dfOnoRep$id,
                                target_dist = target_dist,
                                target_type = c("marginal", "joint", "marginal", "marginal" , "target_data"),
                                boot = 100)

## Checking AMEs
plot_AME(ameOut, factor_name = c("gender", "race"),
         plot_difference = "none",
         plot_type  = c("Exp_Mar", "Exp-Joint", "Rep_Mar", "Dem_Mar", "Full-Joint"),
         col = c("black",  "gray", "red", "blue", "green"))

plot_AME(ameOut_reg, factor_name = c("gender", "race"),
         plot_difference = "none",
         plot_type  = c("Exp_Mar", "Exp-Joint", "Rep_Mar", "Dem_Mar", "Full-Joint"),
         col = c("black",  "gray", "red", "blue", "green"))

## Naoki is working on extending "decomposition" and "conditional_effect_plot". Please wait.
