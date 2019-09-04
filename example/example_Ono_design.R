rm(list=ls())

## Please install "genlasso" by yourself
library(PopCon)
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



# ########################################
# Ono Example
# ########################################
sample_s <- 2000

# ########################
# Interested in One Factor
# ########################
## Unif
ess_unif <- weightDesign(factor_name = "gender",
                     target_dist = marginal_rep, randomize_dist = exp_marginal,
                     sample = sample_s, )
ess_unif$effective_ss

## Pop
ess_pop <- weightDesign(factor_name = "gender",
                         target_dist = marginal_rep, randomize_dist = marginal_rep,
                         sample = sample_s)
ess_pop$effective_ss

## Mixed
mixed_marginal <- marginal_rep
mixed_marginal[["gender"]]  <- exp_marginal[["gender"]]
ess_mix <- weightDesign(factor_name = "gender",
                        target_dist = marginal_rep, randomize_dist = mixed_marginal,
                        sample = sample_s)
ess_mix$effective_ss

# ########################
# Interested in K Factors
# ########################
ess_unif2 <- AverageESS(factor_all = c("gender",  "race"),  type  =  "Unif", sample_s = 2000,
                        target_marginal = marginal_rep, exp_marginal  = exp_marginal)

ess_pop2 <- AverageESS(factor_all = c("gender",  "race"),  type  =  "Pop", sample_s = 2000,
                       target_marginal = marginal_rep, exp_marginal  = exp_marginal)

ess_mixed2 <- AverageESS(factor_all = c("gender",  "race"),  type  =  "Mixed", sample_s = 2000,
                         target_marginal = marginal_rep, exp_marginal  = exp_marginal)

ess_unif2
ess_pop2
ess_mixed2

ess_unif3 <- AverageESS(factor_all = c("gender",  "race",  "pos_security"),  type  =  "Unif", sample_s = 2000,
                        target_marginal = marginal_rep, exp_marginal  = exp_marginal)

ess_pop3 <- AverageESS(factor_all = c("gender",  "race", "pos_security"),  type  =  "Pop", sample_s = 2000,
                       target_marginal = marginal_rep, exp_marginal  = exp_marginal)

ess_mixed3 <- AverageESS(factor_all = c("gender",  "race", "pos_security"),  type  =  "Mixed", sample_s = 2000,
                         target_marginal = marginal_rep, exp_marginal  = exp_marginal)

ess_unif3
ess_pop3
ess_mixed3

AverageESS <- function(factor_all, type =  "Unif",  sample_s, target_marginal, exp_marginal){

  ## Unif
  if(type == "Unif"){
    ess <- c()
    for(i in  1:length(factor_all)){
      ess[i] <- weightDesign(factor_name = factor_all[i],
                             target_dist = target_marginal, randomize_dist = exp_marginal,
                             sample = sample_s)$effective_ss
    }
  }

  if(type == "Pop"){
    ess <- c()
    for(i in  1:length(factor_all)){
      ess[i] <- weightDesign(factor_name = factor_all[i],
                             target_dist = target_marginal, randomize_dist = target_marginal,
                             sample = sample_s)$effective_ss
    }
  }

  if(type == "Mixed"){
    ess <- c()
    mixed_marginal <- target_marginal
    mixed_marginal[factor_all]  <- exp_marginal[factor_all]
    for(i in  1:length(factor_all)){
      ess[i] <- weightDesign(factor_name = factor_all[i],
                             target_dist = target_marginal, randomize_dist = mixed_marginal,
                             sample = sample_s)$effective_ss
    }
  }
  return(ess)
}







# Uniform Randomized Design


# Population Randomized Design
ess2 <- weightDesign(factor_name = "gender", target_dist = marginal_rep, randomize_dist = marginal_rep,
                     sample = sample_s)
plot(density(ess2$ss))

ess2_race <- weightDesign(factor_name = "race", target_dist = marginal_rep, randomize_dist = marginal_rep,
                          sample = sample_s)

ess2_race$effective_ss
ess2$effective_ss

ess2_race$effective_ss*(4/6) + ess2$effective_ss*(2/6)

# Mixed Randomized Design (only Gender)
mixed_marginal <- marginal_rep
mixed_marginal[["gender"]]  <- exp_marginal[["gender"]]

ess3 <- weightDesign(factor_name = "gender", target_dist = marginal_rep, randomize_dist = mixed_marginal,
                     sample = sample_s)
plot(density(ess3$ss))

# Mixed Randomized Design (Gender + Race)
mixed_marginal2 <- marginal_rep
mixed_marginal2[["gender"]]  <- exp_marginal[["gender"]]
# mixed_marginal2[["race"]]    <- exp_marginal[["race"]]
mixed_marginal2[["pos_security"]]    <- exp_marginal[["pos_security"]]

ess4 <- weightDesign(factor_name = "gender", target_dist = marginal_rep, randomize_dist = mixed_marginal2,
                     sample = sample_s)
plot(density(ess4$ss))

ess4_race <- weightDesign(factor_name = "race", target_dist = marginal_rep, randomize_dist = mixed_marginal2,
                          sample = sample_s)
plot(density(ess4_race$ss))

ess4_race$effective_ss*(4/6) + ess4$effective_ss*(2/6)

# Mixed Randomized Design (Gender + Race + Abortion Policy)
mixed_marginal3 <- marginal_rep
mixed_marginal3[["gender"]]  <- exp_marginal[["gender"]]
mixed_marginal3[["pos_security"]]    <- exp_marginal[["pos_security"]]
mixed_marginal3[["pos_abortion"]]    <- exp_marginal[["pos_abortion"]]

ess5 <- weightDesign(factor_name = "gender", target_dist = marginal_rep, randomize_dist = mixed_marginal3,
                     sample = sample_s)
plot(density(ess5$ss))

ess4_race <- weightDesign(factor_name = "race", target_dist = marginal_rep, randomize_dist = mixed_marginal2,
                          sample = sample_s)
plot(density(ess4_race$ss))

ess4_race$effective_ss*(4/6) + ess4$effective_ss*(2/6)

# Gender Effect
ess_mean <- c(ess3$effective_ss, ess4$effective_ss, ess5$effective_ss, ess1$effective_ss, ess2$effective_ss)

plot(seq(1:length(ess_mean)), ess_mean, ylim = c(0, 1),  pch = 19, type  = "o")
abline(h = 1, lty =  2)

sd_prop <- 1/sqrt(ess_mean/sample_s)

plot(seq(1:length(ess_mean)), sd_prop, ylim = c(0, 2.5),  pch = 19, type  = "o")
abline(h = 1, lty =  2)

# ########################################
# Ono Example
# ########################################
sample_s <- 1583*20

## ##########################
## Gender
## ##########################

# Uniform Randomized Design
ess_unif_gender <- weightDesign(factor_name = "gender", target_dist = marginal_rep, randomize_dist = exp_marginal,
                                sample = sample_s)
plot(density(ess_unif_gender$ss))

# Population Randomized Design
ess_pop_gender <- weightDesign(factor_name = "gender", target_dist = marginal_rep, randomize_dist = marginal_rep,
                               sample = sample_s)
plot(density(ess_pop_gender$ss))

# Mixed Randomized Design (only Gender)
mixed_marginal_g1 <- marginal_rep
mixed_marginal_g1[["gender"]]  <- exp_marginal[["gender"]]

ess_mix1_gender <- weightDesign(factor_name = "gender", target_dist = marginal_rep,
                              randomize_dist = mixed_marginal_g1,
                              sample = sample_s)
plot(density(ess_mix1_gender$ss))

# Mixed Randomized Design (Gender + pos_security)
mixed_marginal_g2 <- marginal_rep
mixed_marginal_g2[["gender"]]  <- exp_marginal[["gender"]]
mixed_marginal_g2[["pos_security"]]    <- exp_marginal[["pos_security"]]

ess_mix2_gender <- weightDesign(factor_name = "gender", target_dist = marginal_rep,
                              randomize_dist = mixed_marginal_g2,
                              sample = sample_s)
plot(density(ess_mix2_gender$ss))

# Mixed Randomized Design (Gender + pos_security + pos_immigrants)
mixed_marginal_g3 <- marginal_rep
mixed_marginal_g3[["gender"]]  <- exp_marginal[["gender"]]
mixed_marginal_g3[["pos_security"]]    <- exp_marginal[["pos_security"]]
mixed_marginal_g3[["pos_immigrants"]]    <- exp_marginal[["pos_immigrants"]]

ess_mix3_gender <- weightDesign(factor_name = "gender", target_dist = marginal_rep,
                              randomize_dist = mixed_marginal_g3,
                              sample = sample_s)
plot(density(ess_mix3_gender$ss))

# Effect of Gender
ess_mean_gender <- c(ess_mix1_gender$effective_ss, ess_mix2_gender$effective_ss, ess_mix3_gender$effective_ss,
                     ess_pop_gender$effective_ss, ess_unif_gender$effective_ss)/sample_s

ess_mean_gender <- c(ess_mix1_gender$effective_ss, ess_mix2_gender$effective_ss, ess_mix3_gender$effective_ss,
                     ess_pop_gender$effective_ss, ess_unif_gender$effective_ss)/ess_pop_gender$effective_ss

plot(seq(1:length(ess_mean_gender)), ess_mean_gender, ylim = c(0, 5),  pch = 19, type  = "o")
abline(h = 1, lty =  2)

sd_prop_gender <- 1/sqrt(ess_mean_gender)

plot(seq(1:length(sd_prop_gender)), sd_prop_gender, ylim = c(0, 6.5),  pch = 19, type  = "o")
abline(h = 1, lty =  2)
abline(h = 1.5, lty =  2)

## ##########################
## Race
## ##########################

# Uniform Randomized Design
ess_unif_race <- weightDesign(factor_name = "race", target_dist = marginal_rep, randomize_dist = exp_marginal,
                              sample = sample_s)
plot(density(ess_unif_race$ss))

# Population Randomized Design
ess_pop_race <- weightDesign(factor_name = "race", target_dist = marginal_rep, randomize_dist = marginal_rep,
                             sample = sample_s)
plot(density(ess_pop_race$ss))

# Mixed Randomized Design (only race)
mixed_marginal_r1 <- marginal_rep
mixed_marginal_r1[["race"]]  <- exp_marginal[["race"]]

ess_mix1_race <- weightDesign(factor_name = "race", target_dist = marginal_rep,
                              randomize_dist = mixed_marginal_r1,
                     sample = sample_s)
plot(density(ess_mix1_race$ss))

# Mixed Randomized Design (Race + Gender)
mixed_marginal_r2 <- marginal_rep
mixed_marginal_r2[["gender"]]  <- exp_marginal[["gender"]]
mixed_marginal_r2[["race"]]    <- exp_marginal[["race"]]

ess_mix2_race <- weightDesign(factor_name = "race", target_dist = marginal_rep,
                              randomize_dist = mixed_marginal_r2,
                              sample = sample_s)
plot(density(ess_mix2_race$ss))

# Mixed Randomized Design (Gender + Race + Abortion Policy)
mixed_marginal_r3 <- marginal_rep
mixed_marginal_r3[["gender"]]  <- exp_marginal[["gender"]]
mixed_marginal_r3[["race"]]    <- exp_marginal[["race"]]
mixed_marginal_r3[["pos_immigrants"]]    <- exp_marginal[["pos_immigrants"]]

ess_mix3_race <- weightDesign(factor_name = "race", target_dist = marginal_rep,
                              randomize_dist = mixed_marginal_r3,
                              sample = sample_s)
plot(density(ess_mix3_race$ss))

# Effect of Race
ess_mean_race <- c(ess_mix1_race$effective_ss, ess_mix2_race$effective_ss, ess_mix3_race$effective_ss,
                   ess_pop_race$effective_ss, ess_unif_race$effective_ss)/sample_s

plot(seq(1:length(ess_mean_race)), ess_mean_race, ylim = c(0, 1),  pch = 19, type  = "o")
abline(h = 1, lty =  2)

sd_prop_race <- 1/sqrt(ess_mean_race)

plot(seq(1:length(sd_prop_race)), sd_prop_race, ylim = c(0, 5),  pch = 19, type  = "o")
abline(h = 1, lty =  2)
