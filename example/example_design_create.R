# Replication Code for "Improving the External Validity of Conjoint Analysis" by de la Cuesta, Egami, and Imai
# Our re-analysis of the data from Ono and Burden (2018)
# Code 1: Design-based Analysis
# written on 09/21/2019

rm(list=ls())
library("PopCon")

# Setup
load("../data/ono_data_analysis.rdata")

dfCong <- dfOno[dfOno$office == "Congress", ]
dfPres <- dfOno[dfOno$office == "President", ]
dfList <- list(dfOno, dfCong, dfPres)
names(dfList) <- c("All candidates", "Congressional Candidates", "Presidential Candidates")

# #################################
# Data Generating Process
# #################################
formula_u <- Y ~ gender + age + family + race + experience + party + pos_security

marUse[[1]]  <- marUse[[1]][all.vars(formula_u)[-1]]

featureName <- all.vars(formula_u)[-1]

# Build full two-way interaction formula
factor_l <- length(all.vars(formula_u)[-1])
combMat <- combn(factor_l,2); intNames <- c()
for(k in 1:ncol(combMat)){
  intNames[k] <- paste(all.vars(formula_u)[-1][combMat[1,k]], "*", all.vars(formula_u)[-1][combMat[2,k]], sep = "")
}
formula_full <- as.formula(paste(all.vars(formula_u)[1], "~", paste(intNames, collapse = "+"), sep=""))
fullEq0 <- as.formula(paste("~", paste(intNames, collapse = "+"), sep=""))
formula_u_s <- Y ~ gender + age + family + race + experience + pos_security

# Allow for Cross interactions
base_fac <- all.vars(formula_u)[-1]
for_cross <- paste("~", paste(paste(base_fac, paste(base_fac,"_rp",sep=""), sep="*"),
                              collapse = "+"))
for_cross0 <- paste("~", paste(base_fac, collapse = "+"))

# Function for potential outcome models
TruePF_pair <- function(Data_1, Data_2, true.coef){
  baseX_1 <- model.matrix(as.formula(fullEq0), data= Data_1) # ncol = 135
  Ypred_1 <- cbind(baseX_1) %*% true.coef[1:135]
  baseX_2 <- model.matrix(as.formula(fullEq0), data= Data_2)  # ncol = 135
  Ypred_2 <- cbind(baseX_2) %*% true.coef[1:135]

  # Cross-profile
  Data_2u <- Data_2[,base_fac]; colnames(Data_2u) <- paste(base_fac,"_rp",sep="")
  Data_c <- cbind(Data_1[,base_fac], Data_2u)
  data_cross0 <- model.matrix(as.formula(for_cross0), data = Data_c)[,-1]
  sing <- 2*ncol(data_cross0)
  data_cross <- model.matrix(as.formula(for_cross), data = Data_c)
  data_cross <- data_cross[,-1]
  X_cross <- data_cross[, c((sing + 1): ncol(data_cross))]
  Y_pred_cross <- X_cross %*% true.coef[136:190]
  lin.prob <- (Ypred_1 - Ypred_2) + Y_pred_cross + 0.5
  lin.prob[lin.prob >= 0.99] <- 0.99
  lin.prob[lin.prob <= 0.01] <- 0.01
  outcome <- rbinom(n=length(lin.prob), size = 1, prob = lin.prob)

  return(outcome)
}

# #####################################
# Setup for Randomization
# #####################################
base.level <- lapply(model.frame(formula_full, dfList[[1]])[, -1], levels)
factor_focus <- c("gender") # the special case
factor_remain <- setdiff(featureName, factor_focus)
main_ind <- match(factor_focus, featureName)
control_ind <- match(factor_remain, featureName)
n_fac <- length(featureName)

sim.size <- 100
sim_tab  <- expand.grid("marginal" = c("Rep Pop", "Dem Pop"),
                        "sub" = c("All candidates", "Congressional Candidates", "Presidential Candidates"))

# ########################################
# Start Simulations
# #######################################
set.seed(1234)
point_list <- se_list <- c()
point_list_p <- se_list_p <- c()
for(t in 1:nrow(sim_tab)){

  ## set target distributions
  print(paste("Simulation Design: Type", t))
  prop_each <- marUse[[as.character(sim_tab[t, 1])]]
  prop_each <- prop_each[match(featureName, names(prop_each))]

  target_dist <- list()
  for(i in 1:length(prop_each)){
    target_dist[[i]]<-  prop_each[[i]]
    names(target_dist[[i]]) <-  base.level[[i]]
  }
  names(target_dist) <- names(prop_each)

  data_use  <- dfList[[as.character(sim_tab[t, 2])]]
  sample.size <- nrow(data_use)/2

  ## Set the  true coefficients
  fitAME_base <- model_pAMCE(formula = formula_u,
                             data = data_use,
                             pair = TRUE, pair_id = data_use$pair_id,
                             cross_int = TRUE,
                             target_dist = list(marUse[[1]]),
                             target_type = "marginal", reg = FALSE)
  true.coef.use <- fitAME_base$coef


  ame.Pop.est   <- matrix(NA, nrow = 1, ncol = sim.size)
  ame.Pop.se    <- matrix(NA, nrow = 1, ncol = sim.size)
  ame.Mixed.est <- ame.Mixed.se <- matrix(NA, ncol = sim.size, nrow = 1)


  ## Design 1: Randomize with Mixed Distribution
  cat("\nMixed Randomization: \n")
  for(i in 1:sim.size){
    paste(cat(i), ".", sep = "")
    seed.b <- i*1000
    set.seed(seed.b)

    ame.Mixed_est_b <- ame.Mixed_se_b <- c()

    ## Focus on Gender
    Mixed.r1 <- Mixed.r2 <- as.data.frame(matrix(NA, nrow = sample.size, ncol=n_fac))

    ## Mixed Randomize
    # Uniform for Main Factors
    for(j in main_ind){
      ind_sim_1 <- sample(x = seq(1:length(base.level[[j]])), size = sample.size, replace=TRUE)
      Mixed.r1[,j] <- factor(base.level[[j]][ind_sim_1], levels = base.level[[j]])

      ind_sim_2 <- sample(x = seq(1:length(base.level[[j]])), size = sample.size, replace=TRUE)
      Mixed.r2[,j] <- factor(base.level[[j]][ind_sim_2], levels = base.level[[j]])
    }
    # Pop for Control Factors
    for(j in control_ind){
      ind_sim_1 <- sample(x = seq(1:length(base.level[[j]])), size = sample.size, prob = prop_each[[j]], replace=TRUE)
      Mixed.r1[,j] <- factor(base.level[[j]][ind_sim_1], levels = base.level[[j]])

      ind_sim_2 <- sample(x = seq(1:length(base.level[[j]])), size = sample.size, prob = prop_each[[j]], replace=TRUE)
      Mixed.r2[,j] <- factor(base.level[[j]][ind_sim_2], levels = base.level[[j]])
    }
    colnames(Mixed.r1) <- colnames(Mixed.r2) <- featureName
    Mixed.r1$pair_id <- seq(1:nrow(Mixed.r1))
    Mixed.r2$pair_id <- seq(1:nrow(Mixed.r2))

    # Outcome
    outcome.Mixed <- TruePF_pair(Mixed.r1, Mixed.r2, true.coef.use)
    Mixed.r1$Y <- outcome.Mixed
    Mixed.r2$Y <- 1 - outcome.Mixed

    Mixed.all <- rbind(Mixed.r1, Mixed.r2)
    Mixed.all <- Mixed.all[order(Mixed.all$pair_id), ]
    Mixed.all$id <- data_use$id
    Mixed.all_rp <- rbind(Mixed.r2, Mixed.r1)

    ## Design-based Estimator
    out <- design_pAMCE(formula = formula_u, factor_name = "gender",
                        data = Mixed.all, pair_id = Mixed.all$pair_id,
                        cluster_id = Mixed.all$id,  cross_int = TRUE,
                        target_dist  = target_dist, target_type = "marginal")

    cat(paste("Mixed Num = ", i, "..\n", sep=""))
    ame.Mixed.est[1,i]  <- out$AMCE$gender[1, "estimate"]
    ame.Mixed.se[1,i]   <- out$AMCE$gender[1, "se"]
  }

  Mixed_point <- apply(ame.Mixed.est, 1, mean)
  Mixed_point_se <- apply(ame.Mixed.se, 1, mean)
  Mixed_se    <- apply(ame.Mixed.est, 1, sd)

  point_list[t] <- Mixed_point
  se_list[t]   <- Mixed_point_se


  cat("\nPopulation Randomization: \n")
  for(i in 1:sim.size){
    paste(cat(i), ".", sep = "")
    seed.b <- i*1000
    set.seed(seed.b)

    ame.Pop_est_b <- ame.Pop_se_b <- c()

    ## Focus on Gender
    Pop.r1 <- Pop.r2 <- as.data.frame(matrix(NA, nrow = sample.size, ncol=n_fac))

    # Pop for all Factors
    for(j in 1:n_fac){
      ind_sim_1 <- sample(x = seq(1:length(base.level[[j]])), size = sample.size, prob = prop_each[[j]], replace=TRUE)
      Pop.r1[,j] <- factor(base.level[[j]][ind_sim_1], levels = base.level[[j]])

      ind_sim_2 <- sample(x = seq(1:length(base.level[[j]])), size = sample.size, prob = prop_each[[j]], replace=TRUE)
      Pop.r2[,j] <- factor(base.level[[j]][ind_sim_2], levels = base.level[[j]])
    }
    colnames(Pop.r1) <- colnames(Pop.r2) <- featureName
    Pop.r1$pair_id <- seq(1:nrow(Pop.r1))
    Pop.r2$pair_id <- seq(1:nrow(Pop.r2))

    # Outcome
    outcome.Pop <- TruePF_pair(Pop.r1, Pop.r2, true.coef.use)
    Pop.r1$Y <- outcome.Pop
    Pop.r2$Y <- 1 - outcome.Pop

    Pop.all <- rbind(Pop.r1, Pop.r2)
    Pop.all <- Pop.all[order(Pop.all$pair_id), ]
    Pop.all$id <- data_use$id
    #  save(Pop.all, file = "../data/OnoBurden_data_pr.rdata")

    ## Design-based Estimator
    out_p <- design_pAMCE(formula = formula_u, factor_name = "gender",
                          data = Pop.all, pair_id = Pop.all$pair_id,
                          cluster_id = Pop.all$id, cross_int = TRUE,
                          target_dist  = target_dist, target_type = "marginal")

    # formula
    cat(paste("Pop Num = ", i, "..\n", sep=""))
    ame.Pop.est[1,i]  <- out_p$AMCE$gender[1, "estimate"]
    ame.Pop.se[1,i]   <- out_p$AMCE$gender[1, "se"]
  }

  Pop_point <- apply(ame.Pop.est, 1, mean)
  Pop_point_se <- apply(ame.Pop.se, 1, mean)
  Pop_se    <- apply(ame.Pop.est, 1, sd)

  point_list_p[t] <- Pop_point
  se_list_p[t]   <- Pop_point_se
}

ono_design_result <- list("sim_tab" = sim_tab,
                          "point_list" = point_list,
                          "se_list" = se_list,
                          "point_list_p" = point_list_p,
                          "se_list_p" = se_list_p)
# save(ono_design_result, file = "../result/ono_results_design_final.rdata")



