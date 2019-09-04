rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


## Please install "genlasso" by yourself
library(PopCon)
library(arm)
library("genlasso")
library(sandwich)
library("igraph")
library(prodlim)

# load data
load("ono_data_cleaned_07302019.RData")
load("ono_data_joint_distribution.RData")
full_joint <- dfJoint

# I changed the format of marginal_distritbution. Please check this
## Republican
marginal_rep <- tapply(marDist$`Rep Pop`$prop, marDist$`Rep Pop`$factor, function(x) x)
marginal_rep_n <- tapply(marDist$`Rep Pop`$levels, marDist$`Rep Pop`$factor, function(x) x)
for(i in 1:length(marginal_rep)){
  names(marginal_rep[[i]])  <- marginal_rep_n[[i]]
}
marginal_rep

## Democrat
marginal_dem <- tapply(marDist$`Dem Pop`$prop, marDist$`Dem Pop`$factor, function(x) x)
marginal_dem_n <- tapply(marDist$`Dem Pop`$levels, marDist$`Dem Pop`$factor, function(x) x)
for(i in 1:length(marginal_dem)){
  names(marginal_dem[[i]])  <- marginal_dem_n[[i]]
}

# formula
formula_u <- Y ~ gender + age + family + race + experience + party + trait +
  policy_expertise + pos_security + pos_immigrants + pos_abortion +
  pos_deficit + fav_rating

# Create  Experimental Marginal Distribution
exp_marginal <- createDist(formula_u, target_data = dfOno, exp_data  = dfOno, type = "marginal")

# Target Data (columns should be restricted to what in 'data')
target_data <- full_joint[, is.element(colnames(full_joint), all.vars(formula_u)) == TRUE]
target_data_rep <- full_joint[full_joint$party == "Rep", is.element(colnames(full_joint), all.vars(formula_u)) == TRUE]
target_data_dem <- full_joint[full_joint$party == "Dem", is.element(colnames(full_joint), all.vars(formula_u)) == TRUE]


## Need to have the same level names
target_dist <- list(exp_marginal, marginal_rep, marginal_dem,
                    target_data_rep, target_data_dem)
names(target_dist) <- c("Exp_Mar", "Rep_Mar","Dem_Mar",
                        "Rep_Joint", "Dem_Joint")
target_use <- target_dist[c("Exp_Mar", "Dem_Mar", "Dem_Joint", "Rep_Mar", "Rep_Joint")]

## from NAOKI: Please check levels carefully. Still not the same as in the data...
for(i in c(3,5)){
  levels(target_use[[i]]$age)  <-  levels(dfOno$age)
  levels(target_use[[i]]$race) <-  levels(dfOno$race)
  levels(target_use[[i]]$experience)   <- levels(dfOno$experience)
  levels(target_use[[i]]$pos_security) <- levels(dfOno$pos_security)
  levels(target_use[[i]]$pos_deficit)  <- levels(dfOno$pos_deficit)
}

# Bootstrap Size (just for illustration, I use 10)
boot_size <- 10  # Change this to 1000

## Three-ways between Gender x Party x One of Policy Positions
three_ways <- as.formula(~ gender*party*pos_abortion + gender*party*pos_deficit + gender*party*pos_abortion
                         + gender*party*pos_immigrants)

### Full sample, group 1
dfOno_full <- na.omit(dfOno[, c(all.vars(formula_u), "pair_id", "id")])
ameOut_reg_all <- AME_estimate_full(formula = as.formula(formula_u),
                                    formula_three = three_ways,
                                    data = dfOno_full, type = "genlasso",
                                    pair = TRUE, pair_id = dfOno_full$pair_id,
                                    cluster = dfOno_full$id,
                                    target_dist = target_use,
                                    target_type = c("marginal",
                                                    "marginal", "target_data",
                                                    "marginal", "target_data"),
                                    boot = boot_size, seed = 1234)

## Checking AMEs
plot_AME(ameOut_reg_all, factor_name = c("gender"),
         plot_difference = "none",
         plot_type  = c("STD", "Exp_Mar", "Rep_Mar", "Rep_Joint", "Dem_Mar",  "Dem_Joint"),
         col = c("black", "gray", "tomato", "red", "skyblue1", "blue"),
         main = "All Candidates: \n Include Gender x Party x Policy Positions",
         xlim = c(-0.1, 0.1))

### Congress, group 2
dfOno_c <- dfOno[dfOno$office == "Congress", ]
dfOno_c <- na.omit(dfOno_c[, c(all.vars(formula_u), "pair_id", "id")])
ameOut_reg_c <- AME_estimate_full(formula = as.formula(formula_u),
                                  formula_three = three_ways,
                                  data = dfOno_c, type = "genlasso",
                                  pair = TRUE, pair_id = dfOno_c$pair_id,
                                  cluster = dfOno_c$id,
                                  target_dist = target_use,
                                  target_type = c("marginal",
                                                  "marginal", "target_data",
                                                  "marginal", "target_data"),
                                  boot = boot_size, seed = 1234)

## Checking AMEs
plot_AME(ameOut_reg_c, factor_name = c("gender"),
         plot_difference = "none",
         plot_type  = c("STD", "Exp_Mar", "Rep_Mar", "Rep_Joint", "Dem_Mar",  "Dem_Joint"),
         col = c("black", "gray", "tomato", "red", "skyblue1", "blue"),
         main = "All Candidates: \n Include Gender x Party x Policy Positions",
         xlim = c(-0.1, 0.1))

### Presidential, group 3
dfOno_p <- dfOno[dfOno$office == "President", ]
dfOno_p <- na.omit(dfOno_p[, c(all.vars(formula_u), "pair_id", "id")])
ameOut_reg_p <- AME_estimate_full(formula = as.formula(formula_u),
                                  formula_three = three_ways,
                                  data = dfOno_p, type = "genlasso",
                                  pair = TRUE, pair_id = dfOno_p$pair_id,
                                  cluster = dfOno_p$id,
                                  target_dist = target_use,
                                  target_type = c("marginal",
                                                  "marginal", "target_data",
                                                  "marginal", "target_data"),
                                  boot = boot_size, seed = 1234)

## Checking AMEs
plot_AME(ameOut_reg_p, factor_name = c("gender"),
         plot_difference = "none",
         plot_type  = c("STD", "Exp_Mar", "Rep_Mar", "Rep_Joint", "Dem_Mar",  "Dem_Joint"),
         col = c("black", "gray", "tomato", "red", "skyblue1", "blue"),
         main = "All Candidates: \n Include Gender x Party x Policy Positions",
         xlim = c(-0.1, 0.1))
