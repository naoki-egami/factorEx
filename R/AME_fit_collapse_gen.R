# (Internal Functions) Help functions for AME_estimate_collapse_genlasso.R
globalVariables("i")
AME.collapse.genlasso.crossfit.boot <- function(formula,
                                                data,
                                                pair = FALSE, cross_int,
                                                fac.level, ord.fac,
                                                cv.type = "cv.1Std",
                                                nfolds = 2,
                                                marginal_dist,
                                                marginal_type,
                                                joint_dist_u_list,
                                                formula_three_c,
                                                difference = FALSE,
                                                boot = 100,
                                                tableAME_base,
                                                coefAME_base_l,
                                                eps = 0.0001,
                                                numCores,
                                                seed){



  factor_l <- length(all.vars(formula)[-1])
  combMat <- combn(factor_l,2); intNames <- c()
  for(k in 1:ncol(combMat)){
    intNames[k] <- paste(all.vars(formula)[-1][combMat[1,k]], "*", all.vars(formula)[-1][combMat[2,k]], sep = "")
  }
  # formula_full <- as.formula(paste(all.vars(formula)[1], "~", paste(intNames, collapse = "+"), sep=""))
  formula_full0 <- paste(all.vars(formula)[1], "~", paste(intNames, collapse = "+"), sep="")
  if(is.null(formula_three_c) == TRUE){
    formula_full <- as.formula(formula_full0)
    combMat3 <- NULL

  }else{
    formula_full <- as.formula(paste(formula_full0, formula_three_c, sep = "+"))

    bterm <- terms(as.formula(paste("~", formula_three_c, sep = "")))
    bterm2 <- attr(bterm, "factors")[, attr(bterm, "order") == 3]
    row_num <- match(rownames(attr(bterm, "factors")), all.vars(formula)[-1])
    combMat3 <- row_num*bterm2
    if(is.matrix(combMat3)==TRUE) combMat3 <- apply(combMat3, 2, function(x) sort(x[x!=0]))
    else{ combMat3 <- matrix(sort(combMat3), ncol = 1, nrow = 3)}
  }

  # setup beta weights
  beta_weight <- makeWeight(formula = formula, data = data, pair = pair,
                            fac.level = fac.level, ord.fac = ord.fac)


  message("Cross-Validation: ", appendLF = FALSE)
  set.seed(seed)
  cv.lambda <- col.base.genlasso(formula = formula,
                                 data = data, pair = pair,
                                 fac.level = fac.level, ord.fac = ord.fac,
                                 beta_weight = beta_weight)

  cv.fit <- cv.genlasso(formula = formula,
                        data = data, pair = pair,
                        cv.lambda = cv.lambda,
                        fac.level = fac.level, ord.fac = ord.fac,
                        nfolds = nfolds,
                        cv.type = cv.type,
                        beta_weight = beta_weight,
                        seed = seed)

  lambda <- cv.fit

  # if(pair == TRUE){
  #   remove_id <- names(table(data$cluster))[table(data$cluster) %% 2 != 0]
  #   data <- data[!(data$cluster %in% remove_id), ]
  # }
  data <- data[order(data$cluster), ]

  message(paste("\nBootstrap (", boot, "):\n", sep = ""))
  fit.mat  <- c()
  coef.mat <- matrix(NA, nrow = boot, ncol = length(coefAME_base_l))
  all_eq <- all(table(data$cluster) == table(data$cluster)[1])


  if(Sys.info()[['sysname']] == 'Windows') {

    # ------------
    # Windows
    # ------------

    if (numCores == 1){
      fit_boot <- pblapply(1:boot, function(x)
        crossFitPar(x,
                    formula = formula,
                    formula_full = formula_full,
                    data = data,
                    pair = pair, cross_int = cross_int,
                    fac.level = fac.level, ord.fac = ord.fac,
                    combMat3 = combMat3,
                    lambda = lambda,
                    marginal_dist = marginal_dist,
                    marginal_type = marginal_type,
                    joint_dist_u_list =  joint_dist_u_list,
                    difference = difference,
                    tableAME_base = tableAME_base,
                    coefAME_base_l =  coefAME_base_l,
                    eps = eps,
                    beta_weight = beta_weight,
                    all_eq = all_eq,
                    seed = seed))
    }else {

      cl <- makeCluster(numCores)
      registerDoParallel(cl)
      # pb <- txtProgressBar(max = boot, style = 3)
      # progress <- function(n) setTxtProgressBar(pb, n)
      # opts <- list(progress = progress)
      #
      #     cl <- makeCluster(numCores)
      #     registerDoParallel(cl)
      #     on.exit(stopCluster(cl))

      fit_boot <- foreach(i = 1:boot,
                          .export = c("prepare_data"),
                          .packages = c("genlasso", "prodlim")) %dopar% {
                            crossFitPar(x = i,
                                        formula = formula,
                                        formula_full = formula_full,
                                        data = data,
                                        pair = pair, cross_int = cross_int,
                                        fac.level = fac.level, ord.fac = ord.fac,
                                        combMat3 = combMat3,
                                        lambda = lambda,
                                        marginal_dist = marginal_dist,
                                        marginal_type = marginal_type,
                                        joint_dist_u_list =  joint_dist_u_list,
                                        difference = difference,
                                        tableAME_base = tableAME_base,
                                        coefAME_base_l =  coefAME_base_l,
                                        eps = eps,
                                        beta_weight = beta_weight,
                                        all_eq = all_eq,
                                        seed = seed)
                          }
      #close(pb)
      stopCluster(cl)
    }
  }else{

    # ------------
    # Mac
    # ------------
    fit_boot <- pbmclapply(seq(1:boot), function(x) crossFitPar(x,
                                                                formula = formula,
                                                                formula_full = formula_full,
                                                                data = data,
                                                                pair = pair, cross_int = cross_int,
                                                                fac.level = fac.level, ord.fac = ord.fac,
                                                                combMat3 = combMat3,
                                                                lambda = lambda,
                                                                marginal_dist = marginal_dist,
                                                                marginal_type = marginal_type,
                                                                joint_dist_u_list =  joint_dist_u_list,
                                                                difference = difference,
                                                                tableAME_base = tableAME_base,
                                                                coefAME_base_l =  coefAME_base_l,
                                                                eps = eps,
                                                                beta_weight = beta_weight,
                                                                all_eq = all_eq,
                                                                seed = seed),
                           mc.cores = numCores)
  }


  for(b in 1:boot){
    coef.mat[b, 1:length(coefAME_base_l)] <- fit_boot[[b]]$coef
    fit.mat <- cbind(fit.mat, fit_boot[[b]]$tableAME_full[,4])
  }
  fit <- fit_boot[[1]]$tableAME_full
  estimate <- apply(fit.mat, 1, mean)
  se <- apply(fit.mat, 1, sd)
  low.95ci <- apply(fit.mat, 1, function(x) quantile(x, 0.025))
  high.95ci <- apply(fit.mat, 1, function(x) quantile(x, 0.975))

  fit[,4] <- estimate
  fit$se <- se
  fit$low.95ci <- low.95ci
  fit$high.95ci <- high.95ci

  out <- list("fit" = fit, "fit.mat" = fit.mat, "coef.mat" = coef.mat)
  return(out)
}

crossFitPar <- function(x,
                        formula,
                        formula_full,
                        data,
                        pair, cross_int,
                        fac.level, ord.fac,
                        combMat3,
                        lambda,
                        marginal_dist,
                        marginal_type,
                        joint_dist_u_list,
                        difference,
                        tableAME_base,
                        coefAME_base_l,
                        eps,
                        beta_weight,
                        all_eq,
                        seed){

  seed.b <- 1000*x + seed
  set.seed(seed.b)
  boot_id <- sample(unique(data$cluster), size = length(unique(data$cluster)), replace=TRUE)
  # create bootstap sample with sapply
  boot_which <- sapply(boot_id, function(x) which(data$cluster == x))
  if(all_eq == TRUE){new_boot_id <- rep(seq(1:length(boot_id)), each = table(data$cluster)[1])
  }else{new_boot_id <- rep(seq(1:length(boot_id)), times = unlist(lapply(boot_which, length)))}
  data_boot <- data[unlist(boot_which),]
  data_boot$cluster <- new_boot_id
  data_boot$pair_id <- paste0(data_boot$cluster, data_boot$pair_id)

  fit <- AME.collapse.gen.crossfit(formula = formula,
                                   formula_full = formula_full,
                                   data = data_boot,
                                   pair = pair, cross_int = cross_int,
                                   fac.level = fac.level, ord.fac = ord.fac,
                                   combMat3 = combMat3,
                                   lambda = lambda,
                                   marginal_dist = marginal_dist,
                                   marginal_type = marginal_type,
                                   joint_dist_u_list =  joint_dist_u_list,
                                   difference = difference,
                                   tableAME_base = tableAME_base,
                                   coefAME_base_l = coefAME_base_l,
                                   eps = eps,
                                   beta_weight = beta_weight,
                                   seed_use = seed.b)
  return(fit)
}

AME.collapse.gen.crossfit <- function(formula,
                                      formula_full,
                                      data,
                                      pair=FALSE, cross_int,
                                      combMat3 = NULL,
                                      lambda,
                                      fac.level,
                                      ord.fac,
                                      marginal_dist,
                                      marginal_type,
                                      joint_dist_u_list,
                                      difference = FALSE,
                                      tableAME_base,
                                      coefAME_base_l,
                                      eps,
                                      beta_weight,
                                      seed_use){

  set.seed(seed_use)
  data <- data[order(data$cluster), ]

  train_id <- sample(unique(data$cluster), size = floor(length(unique(data$cluster))/2), replace = FALSE)
  test_id  <- setdiff(unique(data$cluster), train_id)

  train_which <- unlist(sapply(train_id, function(x) which(data$cluster == x)))
  test_which  <- unlist(sapply(test_id, function(x) which(data$cluster == x)))

  data_train <- data[train_which, ]
  data_test  <- data[test_which, ]

  # Fit 1
  fit_col_1 <- collapse.fit.gen(formula = formula,
                                data = data_train, pair = pair,
                                lambda = lambda,
                                fac.level = fac.level, ord.fac = ord.fac,
                                eps = eps,
                                beta_weight = beta_weight)

  fitAME_1 <- fit.after.collapse.gen(formula_full = formula_full,
                                     newdata = data_test,
                                     collapse_level = fit_col_1,
                                     pair = pair, cross_int = cross_int,
                                     marginal_dist = marginal_dist,
                                     marginal_type = marginal_type,
                                     joint_dist_u_list =  joint_dist_u_list,
                                     tableAME_base = tableAME_base,
                                     coefAME_base_l = coefAME_base_l,
                                     difference = difference, combMat3 = combMat3)

  tableAME_1 <- fitAME_1$tableAME_new
  coefAME_1  <- fitAME_1$coef

  # Fit 2
  fit_col_2 <- collapse.fit.gen(formula = formula,
                                data = data_test, pair = pair,
                                lambda = lambda,
                                fac.level = fac.level, ord.fac = ord.fac,
                                eps = eps,
                                beta_weight = beta_weight)

  fitAME_2 <- fit.after.collapse.gen(formula_full = formula_full,
                                     newdata = data_train,
                                     collapse_level = fit_col_2,
                                     pair = pair, cross_int = cross_int,
                                     marginal_dist = marginal_dist,
                                     marginal_type = marginal_type,
                                     joint_dist_u_list =  joint_dist_u_list,
                                     tableAME_base = tableAME_base,
                                     coefAME_base_l =  coefAME_base_l,
                                     difference = difference, combMat3 = combMat3)

  tableAME_2 <- fitAME_2$tableAME_new
  coefAME_2  <- fitAME_2$coef

  coefAME_full <- (coefAME_1 + coefAME_2)/2

  # check
  # print(paste("Check:", all(tableAME_1[,3] == tableAME_2[,3]), sep = ""))
  tableAME_1$estimate <- (tableAME_1$estimate + tableAME_2$estimate)/2

  marginal_dist_u_list <- list()
  for(z in 1:length(marginal_dist)){
    marginal_dist_u_list[[z]] <- data.frame(matrix(NA, ncol=0, nrow=nrow(marginal_dist[[z]])))
    marginal_dist_u_list[[z]]$level <- paste(marginal_dist[[z]][,1], marginal_dist[[z]][,2],sep="")
    marginal_dist_u_list[[z]]$prop  <- marginal_dist[[z]][,3]
  }
  marginal_dist_u_base <- marginal_dist_u_list[[1]]

  ## Add STD
  table_STD <- AME.fit.STD.sep(formula = formula,
                               data = data,
                               pair = pair,
                               marginal_dist = marginal_dist,
                               marginal_dist_u_base = marginal_dist_u_base)

  # Insert into the main table
  uniq_fac <- unique(tableAME_1$factor)
  tableAME_full <- c()
  for(i in 1:length(uniq_fac)){
    tableAME_1_main <- tableAME_1[tableAME_1$factor == uniq_fac[i], ]
    uniq_level <- unique(tableAME_1_main$level)
    for(j in 1:length(uniq_level)){
      tableAME_1_main_m <- tableAME_1_main[tableAME_1_main$level == uniq_level[j], ]
      STD_base <- table_STD[table_STD$fac == uniq_fac[i] & table_STD$level == uniq_level[j], ]
      tableAME_1_main_m <- rbind(STD_base, tableAME_1_main_m)

      if(difference == TRUE){
        tableAME_add <- data.frame(matrix(NA, ncol = 0, nrow = length(marginal_type)))
        dif_est  <- tableAME_1_main_m[2:(1+length(marginal_type)), "estimate"] - STD_base$estimate
        tableAME_add$type <- paste(marginal_type,"-sample AMCE",sep="")
        tableAME_add$factor <- rep(uniq_fac[i], length(marginal_type))
        tableAME_add$level <- rep(uniq_level[j], length(marginal_type))
        tableAME_add$estimate <- dif_est
        tableAME_1_main_m <- rbind(tableAME_1_main_m, tableAME_add)
      }
      tableAME_full <- rbind(tableAME_full, tableAME_1_main_m)
    }
  }

  out <- list("tableAME_full" = tableAME_full, "coef" = coefAME_full)
  return(out)
}


# Collapsing (only main effects)
collapse.fit.gen <- function(formula,
                             data, pair = TRUE,
                             lambda,
                             fac.level, ord.fac,
                             eps = 0.0001,
                             beta_weight){

  beta <- col.genlasso(formula = formula,
                       data = data, pair = pair,
                       lambda = lambda,
                       fac.level = fac.level, ord.fac = ord.fac,
                       beta_weight = beta_weight)

  fac.name <- all.vars(formula)[-1]
  # collapsing
  collapse_level <-  Collapse.genlasso(beta = beta, fac.level = fac.level,
                                       ord.fac = ord.fac,
                                       fac.name = fac.name, eps = eps)$collapse.level

  # adjust (if collapse everything, make it binary)
  collapse_level <- lapply(collapse_level, FUN = function(x) if(length(unique(x)) == 1) c(1, rep(2, length(x) - 1)) else x)

  return(collapse_level)

}


fit.after.collapse.gen <- function(formula_full,
                                   newdata,
                                   collapse_level,
                                   pair=FALSE, cross_int,
                                   marginal_dist,
                                   marginal_type,
                                   joint_dist_u_list,
                                   tableAME_base,
                                   coefAME_base_l,
                                   difference = FALSE,
                                   combMat3 = NULL){

  original_level <- lapply(model.frame(formula_full, data = newdata)[,-1], levels)

  c_data_mar <- prepare_data(formula_full, data = newdata,
                             marginal_dist = marginal_dist,
                             original_level = original_level,
                             collapse_level = collapse_level)

  collapse_level_name <- lapply(model.frame(formula_full, c_data_mar$data_new)[,-1], levels)

  # Transform marginal_dist (for internal simplisity) ----------
  marginal_dist_u_list <- list()
  for(z in 1:length(marginal_dist)){
    marginal_dist_u_list[[z]] <- data.frame(matrix(NA, ncol=0, nrow=nrow(marginal_dist[[z]])))
    marginal_dist_u_list[[z]]$level <- paste(marginal_dist[[z]][,1], marginal_dist[[z]][,2],sep="")
    marginal_dist_u_list[[z]]$prop  <- marginal_dist[[z]][,3]
  }
  marginal_dist_u_base <- marginal_dist_u_list[[1]]

  three_way <- is.null(combMat3) == FALSE

  # Just Get Coefficients
  fitAME <- AME.fit(formula_full,
                    data = c_data_mar$data_new, pair = pair, cross_int = cross_int,
                    marginal_dist = NULL,
                    marginal_dist_u_list = NULL,
                    marginal_dist_u_base = NULL,
                    marginal_type = NULL,
                    difference = difference,
                    three_way = three_way,
                    joint_dist_u_list = NULL,
                    est_AME = FALSE)

  # tableAME <- fitAME$table_AME
  coefAME  <- fitAME$coef
  ind_b    <- fitAME$ind_b

  # Expand Coefficients
  n_fac <- length(all.vars(formula_full)) - 1
  # For main effects
  coefAME_main <- coefAME[1]
  for(z in 1:n_fac){
    coefAME_sub  <- coefAME[ind_b == z]
    collapse_level_b <- collapse_level[[z]][-1]
    coefAME_m0 <- c(0, coefAME_sub)[collapse_level_b]
    ## coefAME_m0 <- c(rep(0, times = sum(collapse_level[[z]] == 1) - 1), coefAME_sub[collapse_level[[z]] - 1]) (only for ordered collapsing)
    coefAME_main <- c(coefAME_main, coefAME_m0)
  }

  # For Two-way Interaction effects (within profiles)
  start <- n_fac
  combMat <- combn(n_fac, 2)
  coefAME_int <- c()
  for(z in 1:ncol(combMat)){
    coefAME_sub <- coefAME[ind_b == (z + start)]
    c_1 <- seq(from = 2, to = max(collapse_level[[combMat[1,z]]]))
    c_2 <- seq(from = 2, to = max(collapse_level[[combMat[2,z]]]))
    c_ind <- paste(rep(c_1, times = length(c_2)), rep(c_2, each = length(c_1)), sep = "_")

    l_ind <- paste(rep(collapse_level[[combMat[1,z]]][-1], times = length(collapse_level[[combMat[2,z]]]) - 1),
                   rep(collapse_level[[combMat[2,z]]][-1], each = length(collapse_level[[combMat[1,z]]]) - 1),
                   sep = "_")
    coefAME_i0 <- coefAME_sub[match(l_ind, c_ind)]
    coefAME_i0[is.na(coefAME_i0)] <- 0
    coefAME_int <- c(coefAME_int, coefAME_i0)
  }
  coefAME_long <- c(coefAME_main, coefAME_int)
  start <- start + ncol(combMat)

  # For Three-way Interaction effects (within profiles)
  if(three_way == TRUE){

    coefAME_int3 <- c()
    for(z in 1:ncol(combMat3)){
      coefAME_sub <- coefAME[ind_b == (z + start)]
      c_1 <- seq(from = 2, to = max(collapse_level[[combMat3[1,z]]]))
      c_2 <- seq(from = 2, to = max(collapse_level[[combMat3[2,z]]]))
      c_3 <- seq(from = 2, to = max(collapse_level[[combMat3[3,z]]]))
      c_ind <- paste(rep(c_1, times = length(c_2)), rep(c_2, each = length(c_1)), sep = "_")

      c_ind <- paste(rep(c_1, times = (length(c_2)*length(c_3))),
                     rep(rep(c_2, each = length(c_1)), times = length(c_3)),
                     rep(c_3, each = (length(c_1)*length(c_2))),
                     sep = "_")
      l_1 <- collapse_level[[combMat3[1,z]]][-1]
      l_2 <- collapse_level[[combMat3[2,z]]][-1]
      l_3 <- collapse_level[[combMat3[3,z]]][-1]
      l_ind <- paste(rep(l_1, times = (length(l_2)*length(l_3))),
                     rep(rep(l_2, each = length(l_1)), times = length(l_3)),
                     rep(l_3, each = (length(l_1)*length(l_2))),
                     sep = "_")

      coefAME_i0 <- coefAME_sub[match(l_ind, c_ind)]
      coefAME_i0[is.na(coefAME_i0)] <- 0
      coefAME_int3 <- c(coefAME_int3, coefAME_i0)
    }
    coefAME_long <- c(coefAME_long, coefAME_int3)
    start <- start + ncol(combMat3)
  }

  # For Interaction effects (within profiles)
  if(cross_int == TRUE){
    coefAME_cross_int <- c()
    for(z in 1:n_fac){
      coefAME_sub <- coefAME[ind_b == (z + start)]
      c_1 <- c_2 <- seq(from = 2, to = max(collapse_level[[z]]))
      c_ind <- paste(rep(c_1, times = length(c_2)), rep(c_2, each = length(c_1)), sep = "_")

      l_ind <- paste(rep(collapse_level[[z]][-1], times = length(collapse_level[[z]]) - 1),
                     rep(collapse_level[[z]][-1], each = length(collapse_level[[z]]) - 1),
                     sep = "_")
      coefAME_i0 <- coefAME_sub[match(l_ind, c_ind)]
      coefAME_i0[is.na(coefAME_i0)] <- 0
      coefAME_cross_int <- c(coefAME_cross_int, coefAME_i0)
    }
    coefAME_long <- c(coefAME_long, coefAME_cross_int)
  }

  names(coefAME_long) <- coefAME_base_l

  # Estimate after Expand
  tableAME_new <- coefIntAME(coefInt = coefAME_long, vcovInt = NULL, SE = FALSE,
                          marginal_dist = marginal_dist, marginal_dist_u_list = marginal_dist_u_list,
                          marginal_dist_u_base = marginal_dist_u_base, marginal_type = marginal_type,
                          difference = difference, cross_int = cross_int, three_way = three_way,
                          joint_dist_u_list = joint_dist_u_list)

  # Expand (tableAME)
  # type_l <- length(unique(tableAME_base$type))
  # fac_name <- unique(tableAME_base$factor)
  # tableAME_new <- matrix(NA, nrow = 0, ncol = 5)
  # for(z in 1:length(fac_name)){
  #   tableAME_sub  <- tableAME[tableAME$factor == fac_name[z],]
  #   tableAME_sub$level_num <- match(tableAME_sub$level, collapse_level_name[[fac_name[z]]])
  #
  #   tableAME_base_sub <- tableAME_base[tableAME_base$factor == fac_name[z],]
  #   tableAME_base_sub$level_num <- collapse_level[[fac_name[z]]][match(tableAME_base_sub$level, original_level[[fac_name[z]]])]
  #
  #   tableAME_use <- merge(tableAME_base_sub, tableAME_sub[, c("type", "level_num", "estimate")],
  #                         by = c("type", "level_num"), all.x = TRUE, all.y = FALSE)
  #   tableAME_use <- tableAME_use[row.match(tableAME_base_sub[, c("type", "level")], tableAME_use[, c("type", "level")]), ]
  #   tableAME_use$estimate[tableAME_use$level_num == 1] <- 0
  #   tableAME_new <- rbind(tableAME_new, tableAME_use)
  # }
  # tableAME_new$level_num <- NULL


  out <- list("tableAME_new" = tableAME_new, "coef" = coefAME_long)

  return(out)
}

# #####################################
# GenLasso Helpers
# #####################################

# GenLasso Collapsing
col.genlasso <- function(formula,
                         data, pair = TRUE,
                         lambda,
                         fac.level, ord.fac,
                         beta_weight){

  # Setup y and X

  if(pair == TRUE){
    data0 <- data[order(data$pair_id),]
    side <- rep(c(1,0), times=nrow(data0)/2)
    data1 <- data0[side==1,]
    data2 <- data0[side==0,]


    X1 <- model.matrix(formula, data=data1)[ ,-1]
    X2 <- model.matrix(formula, data=data2)[ ,-1]
    X <- cbind(1, X1 - X2)

    y <- model.frame(formula,data=data1)[ ,1]
  }else{
    X <- model.matrix(formula, data=data)
    y <- model.frame(formula, data=data)[,1]
  }

  # Setup D
  D <- makeD(fac.level = fac.level, ord.fac = ord.fac)

  # incorporate weights
  D_u <- as.numeric(beta_weight) * D

  # Fit Main
  fit   <- genlasso(X = X, y = y, D = D_u)

  beta_fit <- coef(fit, lambda = lambda)$beta

  return(beta_fit)
}

# (I just need  to have collapsing Results)

Collapse.genlasso <- function(beta, fac.level, ord.fac, fac.name, eps){

  n.fac <- length(fac.level)
  # beta: [1] intercept from [2]
  start_ind <- as.numeric(c(2, c(cumsum(fac.level -1) + 2)[-length(fac.level)]))
  end_ind   <- as.numeric(1 + cumsum(fac.level -1))

  ## Do Collapsing for each factor
  Collapse <- list()
  for(z in 1:length(fac.level)){
    ## First Order

    ## type.ind <- z
    MainD1 <- D1function.genlasso(nlevel = fac.level[z], ord = ord.fac[z])
    MainDif <- t(MainD1 %*% beta[start_ind[z]:end_ind[z]])

    Dif <- abs(MainDif)
    Collapse[[z]] <- apply(Dif <= eps, 2, all)
  }

  ## I got Collapsing Index, then decide which levels will be collapsed.
  collapse.level <- list()
  for(z in 1:n.fac){
    adj <- matrix(0, ncol=fac.level[z], nrow=fac.level[z])
    if(ord.fac[z]==TRUE){
      for(i in 1:length(Collapse[[z]])){
        adj[i, (i+1)] <- as.numeric(Collapse[[z]][i])
      }
      adj <- adj + t(adj)
    }else if(ord.fac[z]==FALSE){
      ref <- combn(seq(1:fac.level[z]),2)
      for(i in 1:length(Collapse[[z]])){
        adj[ref[1,i], ref[2,i]] <- as.numeric(Collapse[[z]][i])
      }
      adj <- adj + t(adj)
    }
    g <- graph_from_adjacency_matrix(adj, mode="undirected")
    collapse.level[[z]] <- components(g)$membership
  }
  names(collapse.level) <- fac.name

  ## Combine two weights.
  output <- list("Collapse.Index" = Collapse, "collapse.level" = collapse.level)
  return(output)
}


cv.genlasso <- function(formula,
                        data, pair = TRUE,
                        cv.lambda,
                        fac.level, ord.fac,
                        nfolds = 5,
                        cv.type = "cv.1Std",
                        beta_weight,
                        seed){

  data <- data[order(data$cluster), ]
  lambda_c <- c()

  set.seed(seed)
  for(i in 1:5){
    train_id <- sample(unique(data$cluster), size = floor(length(unique(data$cluster))/2), replace = FALSE)
    train_which <- unlist(sapply(train_id, function(x) which(data$cluster == x)))
    data_train <- data[train_which, ]

    cv.fit <- cv.genlasso.base(formula = formula,
                               data = data_train, pair = pair,
                               cv.lambda = cv.lambda,
                               fac.level = fac.level, ord.fac = ord.fac,
                               nfolds = nfolds,
                               cv.type = cv.type,
                               beta_weight = beta_weight,
                               seed = (seed + i))
    lambda_c[i] <- cv.fit
    message(paste(round(i*(100/5)),"%..",sep=""), appendLF = FALSE)
  }
  lambda <- mean(lambda_c)

  return(lambda)
}

cv.genlasso.base <- function(formula,
                             data, pair = TRUE,
                             cv.lambda,
                             fac.level, ord.fac,
                             nfolds = 5,
                             cv.type = "cv.1Std",
                             beta_weight,
                             seed){

  # Setup y and X

  if(pair == TRUE){
    data0 <- data[order(data$pair_id),]
    side <- rep(c(1,0), times=nrow(data0)/2)
    data1 <- data0[side==1,]
    data2 <- data0[side==0,]

    X1 <- model.matrix(formula, data=data1)[ ,-1]
    X2 <- model.matrix(formula, data=data2)[ ,-1]
    X <- cbind(1, X1 - X2)

    y <- model.frame(formula,data=data1)[ ,1]
  }else{
    X <- model.matrix(formula, data=data)
    y <- model.frame(formula, data=data)[,1]
  }

  # Setup D
  D <- makeD(fac.level = fac.level, ord.fac = ord.fac)

  # incorporate weights
  D_u <- as.numeric(beta_weight) * D

  # Fit Main
  fit   <- genlasso(X = X, y = y, D = D_u)

  # Cross Validation
  set.seed(seed)
  foldid <- sample(rep(seq(nfolds), length = length(y)))
  MSE <- matrix(0, nrow = nfolds, ncol = length(cv.lambda))
  for (i in seq(nfolds)) {
    which = foldid == i

    X.cv   <- X[!which, ]
    y.cv   <- y[!which]
    X.test <- X[which, ]
    y.test <- y[which]

    fit.cv   <- genlasso(X = X.cv, y = y.cv, D = D_u)
    beta.cv  <- coef(fit.cv, lambda = cv.lambda)$beta
    MSE[i, 1:length(cv.lambda)] <- apply(y.test - X.test%*% beta.cv, 2, function(x) mean(x^2))
  }
  cv.error <- apply(MSE, 2, mean)
  names(cv.error) <- cv.lambda

  # cv.min
  cv.min <- cv.lambda[which.min(cv.error)]

  #cv.sd1
  cv.sd.each <- apply(MSE, 2, sd)
  cv.sd1.value <- min(cv.error) + cv.sd.each[which.min(cv.error)]
  cv.sd1 <- max(cv.lambda[cv.error <= cv.sd1.value])

  if(cv.type == "cv.1Std") lambda_u <- cv.sd1
  if(cv.type == "cv.min")  lambda_u <- cv.min

  return(lambda_u)
}

col.base.genlasso <- function(formula,
                              data, pair = TRUE,
                              fac.level, ord.fac,
                              beta_weight){

  # train_id <- sample(unique(data$cluster), size = floor(length(unique(data$cluster))/2), replace = FALSE)
  # train_which <- unlist(sapply(train_id, function(x) which(data$cluster == x)))
  # data_h <- data[train_which, ]
  data_h <- data

  # Setup y and X
  if(pair == TRUE){
    data0 <- data_h[order(data_h$pair_id),]
    side <- rep(c(1,0), times=nrow(data0)/2)
    data1 <- data0[side==1,]
    data2 <- data0[side==0,]


    X1 <- model.matrix(formula, data=data1)[ ,-1]
    X2 <- model.matrix(formula, data=data2)[ ,-1]
    X <- cbind(1, X1 - X2)


    y <- model.frame(formula,data=data1)[ ,1]
  }else{
    X <- model.matrix(formula, data=data_h)
    y <- model.frame(formula, data=data_h)[,1]
  }

  # Setup D
  D <- makeD(fac.level = fac.level, ord.fac = ord.fac)

  # incorporate weights
  D_u <- as.numeric(beta_weight) * D

  # Fit Main
  fit   <- genlasso(X = X, y = y, D = D_u)

  cv.lambda <- fit$lambda
  return(cv.lambda)
}

center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

makeD <- function(fac.level, ord.fac){
  # colsize = 1 (intercept) + sum(fac.level)
  col.size <- 1 + sum(fac.level - 1)
  D <- matrix(0, nrow = 0, ncol = col.size)
  for(i in 1:length(fac.level)){
    uD <- D1function.genlasso(fac.level[i], ord.fac[i])
    uD_a <- augD(uD, ind = i, fac.level)
    D <- rbind(D, uD_a)
  }
  return(D)
}

augD <- function(mat, ind, fac.level){
  fac.level.u <- fac.level - 1

  if(ind == 1){right.D <- 0
  }else {right.D <- matrix(0, nrow = nrow(mat), ncol = (1 + sum(fac.level.u[1:(ind-1)])))}
  if(ind == length(fac.level)){left.D <- c()
  }else {left.D <- matrix(0, nrow = nrow(mat), ncol = sum(fac.level.u[(ind+1):length(fac.level.u)]))}
  aD <- cbind(right.D, mat, left.D); colnames(aD) <- NULL
  return(aD)
}

D1function.genlasso <- function(nlevel,ord = FALSE){
  if(ord==FALSE){
    qp = (nlevel)*(nlevel-1)/2
    D1.mat <- matrix(0, ncol=(nlevel), nrow=0)
    if (qp == 1) D1.mat <- matrix(c(-1,1), ncol=2, nrow=1)
    if (qp>1){
      for(i in 1 : (nlevel-1)){
        w.diag <- diag((nlevel-i))
        left.v <- c(rep(0,(i-1)),-1)
        left.mat <- matrix(rep(left.v,(nlevel-i)),nrow=(nlevel-i),byrow=TRUE)
        w.mat <- cbind(left.mat,w.diag)
        D1.mat <- rbind(D1.mat,w.mat)
      }
    }
  }else if(ord==TRUE){
    if(nlevel==2){
      D1.mat <- matrix(c(-1,1), ncol=2,nrow=1)
    }else{
      D1.mat <- cbind(0, diag(nlevel-1)) + cbind(-diag(nlevel-1),0)
    }

  }
  D1.final <- as.matrix(D1.mat[, -1])
  return(D1.final)
  ## nrow = differences between main effects
  ## ncol = mainlevel (coefficients for the main factor)
}

# multiply D by column
## number of samples, the number of factors,ordered
## usual estimates
## We only collapase with the first order model

makeWeight <- function(formula,
                       data, pair = TRUE,
                       fac.level, ord.fac){

  # Setup y and X
  if(pair == TRUE){
    data0 <- data[order(data$pair_id),]
    side <- rep(c(1,0), times=nrow(data0)/2)
    data1 <- data0[side==1,]
    data2 <- data0[side==0,]


    X1 <- model.matrix(formula, data=data1)[ ,-1]
    X2 <- model.matrix(formula, data=data2)[ ,-1]
    X <- cbind(1, X1 - X2)

    y <- model.frame(formula,data=data1)[ ,1]
  }else{
    X <- model.matrix(formula, data=data)
    y <- model.frame(formula, data=data)[,1]
  }

  # Fit the model ----------
  main_lm <- lm(y ~ X - 1)
  coefInt <- coef(main_lm)
  coefInt <- coefInt[is.na(coefInt) == FALSE]

  # Setup D
  D <- makeD(fac.level = fac.level, ord.fac = ord.fac)

  # original differences
  w_orig <- 1/abs(D %*% coefInt)

  # Standard weights
  s_size <- lapply(model.frame(formula, data=data)[, -1], function(x) table(x))
  w_add  <- c()
  for(i in 1:length(fac.level)){
    if(ord.fac[i] == TRUE){ # ordered
      s_l <- length(s_size[[i]])
      w_each <- sqrt(s_size[[i]][1:(s_l-1)] + s_size[[i]][2:s_l])
      w_add  <- c(w_add, w_each)
    }else { # unordered
      s_l <- length(s_size[[i]])
      combMat <- combn(s_l,2)
      w_each <- sqrt(s_size[[i]][combMat[1,]] + s_size[[i]][combMat[2,]])/(1 + s_l)*2
      w_add  <- c(w_add, w_each)
    }
  }
  # Adaptive weights
  w_final <- w_orig*w_add
  w_final <- w_final/mean(w_final)

  return(w_final)
}

