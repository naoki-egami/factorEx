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


{formula_orig <- formula
  fac_size <- length(all.vars(formula)[-1])
  fac_name_orig <- all.vars(formula)[-1]
  fac_name <- paste("f_", seq(1:fac_size), "_f", sep = "")
  formula  <- as.formula(paste("Y ~", paste(fac_name, collapse = "+"), sep = ""))
  rename_fac <- cbind(all.vars(formula_orig), c("Y", fac_name))
  colnames(rename_fac) <- c("original", "internal")

  # Rename levels
  original_level <- lapply(model.frame(formula_orig, data = data)[,-1], FUN = function(x) levels(x))
  internal_level <- list()
  for(i in 1:length(original_level)){
    internal_level[[i]] <- paste("x_", i, "_", seq(1:length(original_level[[i]])), "_x", sep = "")
  }
  names(internal_level) <- fac_name

  # Rename data
  data_orig <- data
  for(i in 1:fac_size){
    levels(data[,fac_name_orig[i]]) <- internal_level[[i]]
  }
  colnames(data)[match(all.vars(formula_orig), colnames(data))] <- c("Y", fac_name)

  # Rename marginal_dist
  marginal_dist_orig <- marginal_dist
  for(z in 1:length(marginal_dist)){
    for(i in 1:fac_size){
      temp <- marginal_dist[[z]][marginal_dist[[z]]$factor == rename_fac[(i+1),1],]
      temp$levels <- internal_level[[i]][match(temp$levels, original_level[[i]])]
      temp$factor <- rename_fac[(i+1),2]
      marginal_dist[[z]][marginal_dist[[z]]$factor == rename_fac[(i+1),1],] <- temp
    }
  }
}




factor_l <- length(all.vars(formula)[-1])
combMat <- combn(factor_l,2); intNames <- c()
for(k in 1:ncol(combMat)){
  intNames[k] <- paste(all.vars(formula)[-1][combMat[1,k]], "*", all.vars(formula)[-1][combMat[2,k]], sep = "")
}
formula_full <- as.formula(paste(all.vars(formula)[1], "~", paste(intNames, collapse = "+"), sep=""))

formula_full




coefInt["age44 years old:age60 years old"]
coefInt["age60 years old:age44 years old"]

coefInt["fav_rating70%:fav_rating43%"]
coefInt["fav_rating43%:fav_rating70%"]

start_time <- Sys.time()

formula_u <- Y ~ age + family + race + experience + trait + party +
  policy_expertise + pos_security + pos_immigrants + pos_abortion +
  pos_deficit + fav_rating + gender
data <- dfOnoRep
exp_dist <- createDist(formula_u, data = data)
marginal_dist <- list(exp_dist)

formula_three <- ~ age*race*family + pos_security*experience*age
formula_c <- paste(as.character(formula)[2], as.character(formula)[1], as.character(formula)[3], collapse = "")

factor_l <- length(all.vars(formula)[-1])
combMat <- combn(factor_l,2); intNames <- c()
for(k in 1:ncol(combMat)){
  intNames[k] <- paste(all.vars(formula)[-1][combMat[1,k]], "*", all.vars(formula)[-1][combMat[2,k]], sep = "")
}
formula_full_c <- paste(all.vars(formula)[1], "~", paste(intNames, collapse = "+"), sep="")
formula_aug <- paste(formula_full_c, as.character(formula_three)[2], sep = "+")

terms(as.formula(formula_full))

ameOut_no1 <- AME_estimate_full(formula = as.formula(formula_u),
                                data = data,
                                pair = TRUE, pair_id = data$pair_id,
                                difference = FALSE,
                                cluster = data$id, cross_int = TRUE,
                                marginal_dist = exp_dist,
                                marginal_type = "Exp",
                                type = "genlasso", boot = 30)

ameOut_no1$AME
ameOut_no1$boot_coef

ameOut_no12 <- AME_estimate_full(formula = as.formula(formula_u),
                                 formula_three = ~ age*trait*race,
                                 data = data,
                                 pair = TRUE, pair_id = data$pair_id,
                                 difference = FALSE,
                                 cluster = data$id, cross_int = TRUE,
                                 marginal_dist = exp_dist,
                                 marginal_type = "Exp",
                                 type = "genlasso", boot = 30)

ameOut_no12$AME
ameOut_no12$Ftest
dim(ameOut_no1$boot_coef)
dim(ameOut_no12$boot_coef)

ameOut_no2 <- AME_estimate_full(formula = as.formula(formula_u),
                                data = data,
                                pair = TRUE, pair_id = data$pair_id,
                                difference = FALSE,
                                cluster = data$id, cross_int = TRUE,
                                marginal_dist = exp_dist,
                                marginal_type = "Exp",
                                type = "No-Reg")

ameOut_no2$AME
ameOut_no2$boot_coef

ameOut_no22 <- AME_estimate_full(formula = as.formula(formula_u),
                                 formula_three = ~ age*trait*race,
                                 data = data,
                                 pair = TRUE, pair_id = data$pair_id,
                                 difference = FALSE,
                                 cluster = data$id, cross_int = TRUE,
                                 marginal_dist = exp_dist,
                                 marginal_type = "Exp",
                                 type = "No-Reg")

ameOut_no22$AME
ameOut_no22$boot_coef

out1 <- do.call("rbind", ameOut_no1$AME)
out12 <- do.call("rbind", ameOut_no12$AME)
out2 <- do.call("rbind", ameOut_no2$AME)
out22 <- do.call("rbind", ameOut_no22$AME)

plot(out1$estimate, out2$estimate)
abline(0, 1)

plot(out2$estimate, out22$estimate)
abline(0, 1)

plot(out12$estimate, out22$estimate)
abline(0, 1)

plot(out1$estimate, out12$estimate)
abline(0, 1)

plot(out1$se, out12$se)
abline(0, 1)

formula = as.formula(formula_u);
formula_three = ~ age*trait*race;
data = data;
pair = TRUE; pair_id = data$pair_id;
difference = FALSE;
cluster = data$id; cross_int = TRUE;
marginal_dist = list(exp_dist);
marginal_type = "Exp";
type = "genlasso"

names(ameOut_no1)
out <- ameOut_no1
## Name Back
for(i in 1:length(out$AME)){
  match_level <- match(out$AME[[i]]$level, internal_level[[out$AME[[i]]$factor[1]]])
  out$AME[[i]]$factor <- rename_fac[,"original"][match(out$AME[[i]]$factor[1], rename_fac[,"internal"])]
  orignal_use <- original_level[[out$AME[[i]]$factor[1]]]
  out$AME[[i]]$level <- orignal_use[match_level]
}
names(out$AME) <- rename_fac[,"original"][match(names(out$AME), rename_fac[,"internal"])]
for(i in 1:fac_size){
  colnames(out$boot_coef) <- gsub(rename_fac[(i+1), "internal"], rename_fac[(i+1), "original"], colnames(out$boot_coef))
  for(j in 1:length(internal_level[[i]])){
    colnames(out$boot_coef) <- gsub(internal_level[[i]][j], original_level[[i]][j], colnames(out$boot_coef))
  }
}




names(ameOut_no1$input)


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
