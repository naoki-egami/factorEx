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
load("ono_data_cleaned.RData")
load("ono_data_joint_distribution.RData")
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
target_data_rep <- full_joint[full_joint$party == "Rep", is.element(colnames(full_joint), all.vars(formula_u)) == TRUE]
target_data_dem <- full_joint[full_joint$party == "Dem", is.element(colnames(full_joint), all.vars(formula_u)) == TRUE]



## Need to have the same level names

levels(target_data$age)  <-  levels(dfOnoRep$age)
levels(target_data$race) <-  levels(dfOnoRep$race)
levels(target_data$experience)   <- levels(dfOnoRep$experience)
levels(target_data$pos_security) <- levels(dfOnoRep$pos_security)
levels(target_data$pos_deficit)  <- levels(dfOnoRep$pos_deficit)

levels(target_data_rep$age)  <-  levels(dfOnoRep$age)
levels(target_data_rep$race) <-  levels(dfOnoRep$race)
levels(target_data_rep$experience)   <- levels(dfOnoRep$experience)
levels(target_data_rep$pos_security) <- levels(dfOnoRep$pos_security)
levels(target_data_rep$pos_deficit)  <- levels(dfOnoRep$pos_deficit)

levels(target_data_dem$age)  <-  levels(dfOnoRep$age)
levels(target_data_dem$race) <-  levels(dfOnoRep$race)
levels(target_data_dem$experience)   <- levels(dfOnoRep$experience)
levels(target_data_dem$pos_security) <- levels(dfOnoRep$pos_security)
levels(target_data_dem$pos_deficit)  <- levels(dfOnoRep$pos_deficit)

## check
rep_check <- lapply(target_data_rep, function(x) c(prop.table(table(x))))

check <- cbind(unlist(marginal_rep), c(unlist(rep_check[match(names(marginal_rep), names(rep_check))]), rep(NA, 6)))
colnames(check) <- c("distList$`Dem Pop", "target_data")

target_dist <- list(exp_marginal, exp_joint, marginal_rep, marginal_dem,  target_data, target_data_rep, target_data_dem)
names(target_dist) <- c("Exp_Mar", "Exp-Joint", "Rep_Mar", "Dem_Mar", "Full-Joint", "Rep-Joint", "Dem-Joint")


## Need to specify "target_type"; one of "marginal", "joint", and "target_data"
## For us, we will only compare "marginal" and "target"
## But for illustration, I show everything here

### Full sample

dfUse <- dfOno
ameOut_reg <- AME_estimate_full(formula = as.formula(formula_u),
                                #formula_three = ~  age*family*race + family*experience*party,
                                data = dfUse, type = "genlasso",
                                pair = TRUE, pair_id = dfUse$pair_id,
                                cluster = dfUse$id,
                                target_dist = target_dist,
                                target_type = c("marginal", "joint", "marginal", "marginal" ,
                                                "target_data", "target_data", "target_data"),
                                boot = 500, seed = 3000)
## Checking AMEs
plot_AME(ameOut_reg, factor_name = c("gender", "race"),
         plot_difference = "none",
         plot_type  = c("Rep_Mar", "Dem_Mar", "Rep-Joint", "Dem-Joint"),
         col = c("black",  "gray", "red", "blue"))

###### Congressional Candidates

dfUse <- dfOno[dfOno$office == "Congress", ]
ameOut_reg_congress <- AME_estimate_full(formula = as.formula(formula_u),
                                         #formula_three = ~  age*family*race + family*experience*party,
                                         data = dfUse, type = "genlasso",
                                         pair = TRUE, pair_id = dfUse$pair_id,
                                         cluster = dfUse$id,
                                         target_dist = target_dist,
                                         target_type = c("marginal", "joint", "marginal", "marginal" ,
                                                         "target_data", "target_data", "target_data"),
                                         boot = 500, seed = 3000)
## Checking AMEs
plot_AME(ameOut_reg_congress, factor_name = c("gender", "race"),
         plot_difference = "none",
         plot_type  = c("Rep_Mar", "Dem_Mar", "Rep-Joint", "Dem-Joint"),
         col = c("black",  "gray", "red", "blue"))


###### Presidential Candidates

dfUse <- dfOno[dfOno$office == "President", ]
ameOut_reg_president <- AME_estimate_full(formula = as.formula(formula_u),
                                          #formula_three = ~  age*family*race + family*experience*party,
                                          data = dfUse, type = "genlasso",
                                          pair = TRUE, pair_id = dfUse$pair_id,
                                          cluster = dfUse$id,
                                          target_dist = target_dist,
                                          target_type = c("marginal", "joint", "marginal", "marginal" ,
                                                          "target_data", "target_data", "target_data"),
                                          boot = 500, seed = 3000)
## Checking AMEs
plot_AME(ameOut_reg_president, factor_name = c("gender", "race"),
         plot_difference = "none",
         plot_type  = c("Rep_Mar", "Dem_Mar", "Rep-Joint", "Dem-Joint"),
         col = c("black",  "gray", "red", "blue"))


#### Compare all side-by-side

par(mfrow = c(1, 3))
fac_inc <- c("gender")
plot_order <- c("Exp-Joint", "Full-Joint", "Rep_Mar", "Rep-Joint", "Dem_Mar", "Dem-Joint")
col_order <- c("gray", "black", "tomato", "red", "skyblue1", "blue")
plot_AME(ameOut_reg, factor_name = fac_inc,
         plot_difference = "none",
         plot_type  = plot_order,
         col = col_order,
         main = "Full Sample")
plot_AME(ameOut_reg_congress, factor_name = fac_inc,
         plot_difference = "none",
         plot_type  = plot_order,
         col = col_order,
         main = "Congressional Candidates")
plot_AME(ameOut_reg_president, factor_name = fac_inc,
         plot_difference = "none",
         plot_type  = plot_order,
         col = col_order,
         main = "Presidential Candidates")




