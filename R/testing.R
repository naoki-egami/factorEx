library("genlasso")
library(PopCon)
library(arm)
require(FindIt)
require(quadprog)
require(prodlim)
require(Hmisc)
require(arm)
require(clubSandwich)
data("OnoData")
OnoRep <- OnoData$OnoRep # Ono data only with Republican subjects
marginal_dist <- OnoData$marginal_dist

formula_u <- Y ~ gender + age + family + race + experience

seed <- 1234
nfolds <- 2
cv.type <- "cv.1Std"
eps <- 0.0001
formula <- as.formula(formula_u)
data = OnoRep
pair = TRUE
pair_id = OnoRep$pair_id
difference = TRUE
cluster = OnoRep$id
marginal_dist = marginal_dist
marginal_type = names(marginal_dist)
type = "No-Reg"


ameOut_no <- AME_estimate_full(formula = as.formula(formula_u),
                                data = OnoRep,
                                pair = TRUE, pair_id = OnoRep$pair_id,
                                difference = TRUE,
                                cluster = OnoRep$id,
                                marginal_dist = marginal_dist,
                                marginal_type = names(marginal_dist),
                                type = "No-Reg",
                                boot = 10, cv.type = "cv.1Std")



ameOut_gen <- AME_estimate_full(formula = as.formula(formula_u),
                                data = OnoRep,
                                pair = TRUE, pair_id = OnoRep$pair_id,
                                difference = TRUE,
                                cluster = OnoRep$id,
                                marginal_dist = marginal_dist,
                                marginal_type = names(marginal_dist),
                                type = "genlasso",
                                boot = 100, cv.type = "cv.1Std")


ameOut_gen_rubin <- AME_estimate_collapse_genlasso_rubinSE(formula = as.formula(formula_u),
                                                     data = OnoRep,
                                                     pair = TRUE, pair_id = OnoRep$pair_id,
                                                     difference = TRUE,
                                                     cluster = OnoRep$id,
                                                     marginal_dist = marginal_dist,
                                                     marginal_type = names(marginal_dist),
                                                     boot = 100, cv.type = "cv.1Std")


ameOut_no_all   <- do.call("rbind", ameOut_no$AME)
ameOut_gen_all  <- do.call("rbind", ameOut_gen$AME)
ameOut_gen_rubin_all  <- do.call("rbind", ameOut_gen_rubin$AME)


plot(ameOut_no_all$estimate, ameOut_gen_rubin_all$estimate, main = "Estimate: No-Reg vs genlasso")
abline(0,1, lty = 2, col = "red")

plot(ameOut_no_all$se, ameOut_gen_rubin_all$se, main = "Estimate: No-Reg vs genlasso")
abline(0,1, lty = 2, col = "red")

plot(ameOut_no_all$estimate, ameOut_gen_all$estimate, main = "Estimate: No-Reg vs genlasso")
abline(0,1, lty = 2, col = "red")

plot(ameOut_gen_all$estimate, ameOut_gen_rubin_all$estimate, main = "Estimate: No-Reg vs genlasso")
abline(0,1, lty = 2, col = "red")
plot(ameOut_gen_all$se, ameOut_gen_rubin_all$se, main = "SE: No-Reg vs genlasso", col = "blue")
abline(0,1, lty = 2, col = "red")
