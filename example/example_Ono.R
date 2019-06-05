rm(list=ls())

## Please install "genlasso" by yourself
library("genlasso")
library(PopCon)
library(arm)
library(sandwich)
library("igraph")

library(prodlim)

load("example/ono_data_cleaned.RData")
data <- dfOno
data <- dfOnoRep
pair_id <- data$pair_id

lm_1 <- lm(Y ~ pos_abortion, data = data)
summary(lm_1)

data  <- data[order(data$pair_id), ]
side <- rep(c(1,0), times = nrow(data)/2)
data1 <- data[side == 1, ]
data0 <- data[side == 0, ]
colnames(data0) <- paste(colnames(data0), "_right", sep = "")
data_long <- cbind(data1, data0)
table(data_long$Y, data_long$Y_left)

lm_2 <- lm(Y ~ pos_abortion*pos_abortion_right, data = data_long)
summary(lm_2)


# Specify Formula
formula_u <- Y ~ gender + age + family + race + experience + trait + party +
  policy_expertise + pos_security + pos_immigrants + pos_abortion +
  pos_deficit + fav_rating

factor_l <- length(all.vars(formula)[-1])
combMat <- combn(factor_l,2); intNames <- c()
for(k in 1:ncol(combMat)){
  intNames[k] <- paste(all.vars(formula)[-1][combMat[1,k]], "*", all.vars(formula)[-1][combMat[2,k]], sep = "")
}
formula_full <- as.formula(paste(all.vars(formula)[1], "~", paste(intNames, collapse = "+"), sep=""))

formula_full


exp_dist <- createDist(formula_u, data = data)

coefInt["age44 years old:age60 years old"]
coefInt["age60 years old:age44 years old"]

coefInt["fav_rating70%:fav_rating43%"]
coefInt["fav_rating43%:fav_rating70%"]

start_time <- Sys.time()

ameOut_no1 <- AME_estimate_full(formula = as.formula(formula_u),
                                data = data,
                                pair = TRUE, pair_id = data$pair_id,
                                difference = FALSE,
                                cluster = data$id, cross_int = TRUE,
                                marginal_dist = exp_dist,
                                marginal_type = "Exp",
                                type = "genlasso")

ameOut_no1$AME

coef1 <- do.call("rbind", ameOut_no1$AME)
std <- coef1[coef1$type ==  "STD",  ]
exp <- coef1[coef1$type ==  "Exp",  ]

plot(coef1[coef1$type ==  "STD",  "estimate"], seq(1:nrow(std)), col = "black", pch = 19)
segments(coef1[coef1$type ==  "STD",  "low.95ci"], seq(1:nrow(std)),
         coef1[coef1$type ==  "STD",  "high.95ci"], seq(1:nrow(std)), col = "black")
points(coef1[coef1$type ==  "Exp",  "estimate"], seq(1:nrow(std)), col = "red", pch = 15)
segments(coef1[coef1$type ==  "Exp",  "low.95ci"], seq(1:nrow(std)),
         coef1[coef1$type ==  "Exp",  "high.95ci"], seq(1:nrow(std)), col = "red")
abline(v = 0, lty = 2)

ameOut_no2 <- AME_estimate_full(formula = as.formula(formula_u),
                                data = data,
                                pair = TRUE, pair_id = data$pair_id,
                                difference = FALSE,
                                cluster = data$id, cross_int = FALSE,
                                marginal_dist = exp_dist,
                                marginal_type = "Exp",
                                type = "No-Reg")

coef2 <- do.call("rbind", ameOut_no2$AME)

plot(coef2[coef2$type ==  "STD",  "estimate"], coef2[coef2$type ==  "Exp",  "estimate"])
abline(0, 1)

ameOut_no1 <- AME_estimate_full(formula = as.formula(formula_u),
                               data = OnoRep,
                               pair = TRUE, pair_id = OnoRep$pair_id,
                               difference = FALSE,
                               cluster = OnoRep$id,
                               marginal_dist = list(exp_dist),
                               marginal_type = "Exp",
                               type = "No-Reg")

ameOut_no2 <- AME_estimate_full(formula = as.formula(formula_u),
                               data = OnoRep, pair  = FALSE,
                               # pair = TRUE, pair_id = OnoRep$pair_id,
                               difference = FALSE,
                               cluster = OnoRep$id,
                               marginal_dist = exp_dist,
                               marginal_type = "Exp",
                               type = "No-Reg")

ameOut_no1$AME$race
ameOut_no2$AME$race

decompose_AME(ameOut_no1, factor_name = "race",  level_name = "Black", marginal_diff = c("Exp"))

ameOut_no1$AME$gender
ameOut_no2$AME$gender

ameOut_no2 <- AME_estimate_full(formula = as.formula(formula_u),
                               data = OnoRep,
                               pair = TRUE, pair_id = OnoRep$pair_id, pair_var = c("experience", "gender"),
                               difference = FALSE,
                               cluster = OnoRep$id,
                               marginal_dist = exp_dist,
                               marginal_type = "Exp",
                               type = "No-Reg")
end_time <- Sys.time()
end_time - start_time

ameOut_no$AME$experience
ameOut_no2$AME$experience

plot(ameOut_no$AME$experience[,4], ameOut_no2$AME$experience[,4])
abline(0,1)

data1 <- OnoRep
formula <- formula_full

# X1 <- model.matrixBayes(formula, data=OnoRep)
# attr(X1, "assign")
#
# Xf <- model.frame(formula, data=data1)
# Xf <- attr(terms(Xf), "factors")
# Xf_ind <- c(apply(matrix(Xf[pair_var, ], ncol = ncol(Xf)), 1, function(x) which(x == 1)))
# X01 <- model.matrix(formula, data=data1)
# X02 <- model.matrix(formula, data=data2)
# pair_var_ind <- is.element(attr(X01, "assign"), Xf_ind)
# X0 <- X01 - X02
# X0[, pair_var_ind == TRUE] <- X01[, pair_var_ind == TRUE]
# X <- cbind(1, X0[,-1])

start_time <- Sys.time()
ameOut_gash <- AME_estimate_full(formula = as.formula(formula_u),
                                 data = OnoRep,
                                 pair = TRUE, pair_id = OnoRep$pair_id,
                                 pair_var = c("gender","experience"),
                                 difference = FALSE,
                                 cluster = OnoRep$id,
                                 marginal_dist = exp_dist,
                                 marginal_type = "Exp",
                                 cv.collapse.cost = c(0.01, 0.05),
                                 type = "gash-anova",
                                 boot = 10, cv.type = "cv.1Std",
                                 family = "binomial")
end_time <- Sys.time()
end_time - start_time

ameOut_gash$AME$experience
ameOut_gash$AME$experience

start_time <- Sys.time()
ameOut_gen <- AME_estimate_full(formula = as.formula(formula_u),
                                 data = OnoRep,
                                 pair = TRUE, pair_id = OnoRep$pair_id, pair_var = c("gender","experience"),
                                 difference = FALSE,
                                 cluster = OnoRep$id,
                                 marginal_dist = exp_dist,
                                 marginal_type = "Exp",
                                 type = "genlasso",
                                 boot = 100, cv.type = "cv.min", seed = 200)
end_time <- Sys.time()
end_time - start_time

ameOut_gen$AME$experience
ameOut_gen$AME$experience

## ---------------------------------------------
## Investigate Conditional Effects Directly
## ----------------------------------------------
cAME_gen <- cAME_from_boot(ameOut_gen, factor_name = "gender", level_name = "Male", difference = FALSE)

plot_cAME(cAME_gen,
          factor_name = c("age", "family", "race",
                          "experience", "trait",
                          "policy_expertise", "pos_security", "pos_immigrants",
                          "pos_abortion","pos_deficit", "fav_rating",
                          "party"),
          col = c("black"),
          legend_pos = "topright",
          plot_all = FALSE,
          plot_difference = "no",
          cex = 1, mar = 10)
