rm(list=ls())

library(PopCon)
# library(arm)
# library("genlasso")
# library(sandwich)
# library("igraph")
# library(prodlim)



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
target_data <- target_data[target_data$party == "Dem", ]

## Need to have the same level names
levels(target_data$age)  <-  levels(dfOnoRep$age)
levels(target_data$race) <-  levels(dfOnoRep$race)
levels(target_data$experience)   <- levels(dfOnoRep$experience)
levels(target_data$pos_security) <- levels(dfOnoRep$pos_security)
levels(target_data$pos_deficit)  <- levels(dfOnoRep$pos_deficit)

target_dist <- list(marginal_rep, marginal_dem,  target_data)
names(target_dist) <- c("Rep_Mar", "Dem_Mar", "Full-Joint")

target_dist_m <- list(marginal_rep, marginal_dem)
names(target_dist_m) <- c("Rep_Mar", "Dem_Mar")

## Need to specify "target_type"; one of "marginal", "joint", and "target_data"
## For us, we will only compare "marginal" and "target"
## But for illustration, I show everything here



# Without regularization
ameOut <- model_pAMCE(formula = as.formula(formula_u),
                      data = dfOnoRep, reg = FALSE,
                      pair = TRUE, pair_id = dfOnoRep$pair_id,
                      cluster = dfOnoRep$id,
                      target_dist = target_dist,
                      target_type = c("marginal", "marginal" , "target_data"),
                      boot = 100)

ameOut2 <- model_pAMCE(formula = as.formula(formula_u),
                      formula_three = ~  gender*family*race + gender*experience*party,
                      data = dfOnoRep, reg = FALSE,
                      pair = TRUE, pair_id = dfOnoRep$pair_id,
                      cluster = dfOnoRep$id,
                      target_dist = target_dist,
                      target_type = c("marginal", "marginal", "target_data"),
                      boot = 100)

# With regularization
ameOut_reg <- model_pAMCE(formula = as.formula(formula_u),
                            formula_three = ~  gender*family*race + family*experience*party,
                          data = dfOnoRep,
                          pair = TRUE, pair_id = dfOnoRep$pair_id,
                          cluster = dfOnoRep$id,
                          target_dist = target_dist,
                          target_type = c("marginal", "marginal" , "target_data"),
                          boot = 10)

ameOut_reg <- model_pAMCE(formula = as.formula(formula_u),
                          data = dfOnoRep,
                          pair_id = dfOnoRep$pair_id, cluster_id = dfOnoRep$id,
                          target_dist = target_dist_m, target_type = c("marginal", "marginal"),
                          boot = 100)


ameOut_reg$AMCE

# summary function
summary(ameOut_reg)
a <- summary(ameOut_reg, factor_name = c("gender"))

# $gender
# type factor level    estimate          se     low.95ci    high.95ci
# 361  Sample gender  Male  0.01876321 0.009760814 -0.003481056 0.0358871347
# 106  sample gender  Male  0.02090719 0.008678511  0.005512503 0.0371071930
# 107 Rep_Mar gender  Male  0.05618836 0.025618133  0.004679538 0.0960902242
# 108 Dem_Mar gender  Male -0.03848595 0.023489553 -0.086385371 0.0005897246

#   ----------------
#   Population AMCEs:
#   ----------------
#   target_dist factor level    Estimate  Std. Error p value
#   Sample gender  Male  0.01963569 0.008592828   0.022 *
#   Rep_Mar gender  Male  0.03343790 0.027948377   0.232
#   Dem_Mar gender  Male -0.05161072 0.023424141   0.028 *
#   Full-Joint gender  Male -0.05708989 0.026911846   0.034 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# plot function
plot(ameOut_reg, mar = 15, xlim = c(-0.1, 0.3), diagnose = TRUE)

x <- ameOut_reg

decompose_pAMCE(ameOut_reg, effect_name = c("gender", "Male"))

plot_decompose(ameOut_reg,
               effect_name = c("gender", "Male"),
               target_diff = c("Dem_Mar",  "sample"))

# summary function
summary(ameOut_reg, factor_name = c("gender"), sample = TRUE)

# plot function
plot(x = ameOut_reg, factor_name = c("gender", "race"),
     diagnose = TRUE, legend_pos = "topleft", xlim = c(-0.2, 0.2))

plot_decompose(ameOut_reg, effect_name = c("gender", "Male"), target_diff = c("Dem_Mar",  "Sample_Mar"))

x <- ameOut_reg

AMCE_table <- do.call("rbind", x$AMCE)
AMCE_table_u <- AMCE_table[AMCE_table$type %in% names(x$input$target_dist), ]
AMCE_table_u$pv <- round(2*(1- pnorm(abs(AMCE_table_u$estimate/AMCE_table_u$se))), 3)
AMCE_print <- AMCE_table_u[, c(1, 2, 3, 4, 5, 8)]

colnames(AMCE_print) <- c("target_dist", "factor", "level", "Estimate", "Std. Error", "p value")
print(AMCE_print, row.names = FALSE)



## Checking AMEs


## Checking AMEs
plot(x = ameOut_reg, factor_name = c("gender", "race"),
     target_dist = c("Sample_Mar", "Dem_Mar", "Rep_Mar"),
     diagnose = TRUE, xlim = c(-0.1,  0.1))

plot_decompose(ameOut, effect_name = c("race", "Black"), target_diff = c("Dem_Mar",  "Sample_Mar"))

plot_diagnose(x  = ameOut_reg, factor_name = c("gender"), legend_pos = "topleft", target_dist = c("Dem_Mar",  "Sample_Mar"))

decom_c <- plot_decompose(ameOut2, factor_name = "gender", level_name = "Male",
                          target_diff = c("Rep_Mar", "Sample_Mar"))

decom_c <- plot_decompose(ameOut, factor_name = "gender", level_name = "Male",
                          target_diff = c("Rep_Mar", "Sample_Mar"))

decom_c <- decompose_pAMCE(ameOut, factor_name = "gender", level_name = "Male",
                         target_diff = c("Rep_Mar", "Sample_Mar"))

decom_c2 <- decompose_pAMCE(ameOut2, factor_name = "gender", level_name = "Male",
                           target_diff = c("Dem_Mar", "Sample_Mar"))

a <- createDist(formula = as.formula(formula_u),
           target_data = target_dist[[4]],
           exp_data = dfOnoRep, type  = "marginal")

a
marginal_dem



plot(ameOut, factor_name = c("gender", "race", "age"),
     target_dist = c("Rep_Mar", "Exp_Mar"),
     xlim = c(-0.6, 0.6))

plot(ameOut, factor_name = c("gender", "race", "age"),
     target_dist = c("Exp_Mar", "Rep_Mar"),
     xlim = c(-0.6, 0.6))

plot_AME(ameOut_reg, factor_name = c("gender", "race"),
         plot_difference = "none",
         plot_type  = c("Exp_Mar", "Exp-Joint", "Rep_Mar", "Dem_Mar", "Full-Joint"),
         col = c("black",  "gray", "red", "blue", "green"))

## Naoki is working on extending "decomposition" and "conditional_effect_plot". Please wait.


## Separate Marginals

formula = as.formula(formula_u)
data = dfOnoRep
dfOnoRep$cluster <- dfOnoRep$id

check.old <- AME.fit.STD(formula  = as.formula(formula_u),
                         data  = dfOnoRep,
                         pair = FALSE,
                         marginal_dist = list(distList$`Rep Pop`),
                         marginal_dist_u_base = marginal_dist_u_base)
check.se.old <- AME.fit.STD.se(formula  = as.formula(formula_u),
                               data  = dfOnoRep,
                               pair = FALSE,
                               marginal_dist = list(distList$`Rep Pop`),
                               marginal_dist_u_base = marginal_dist_u_base)


check <- AME.fit.STD.sep(formula  = as.formula(formula_u),
                         data  = dfOnoRep,
                         pair = FALSE,
                         marginal_dist = list(distList$`Rep Pop`),
                         marginal_dist_u_base = marginal_dist_u_base)
check.se <- AME.fit.STD.se.sep(formula  = as.formula(formula_u),
                               data  = dfOnoRep,
                               pair = FALSE,
                               marginal_dist = list(distList$`Rep Pop`),
                               marginal_dist_u_base = marginal_dist_u_base)

plot(check[, 4], check.se[, 4])
abline(0, 1)

plot(check.old[, 4], check.se.old[, 4])
abline(0, 1)

plot(check[,4], check.old[,4])
abline(0, 1)

plot(check.se[,4], check.se.old[,4])
abline(0, 1)

plot(check.se[,5], check.se.old[,5])
abline(0, 1)
