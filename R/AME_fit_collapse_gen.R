collapase.fit.gen <- function(formula,
                              data, pair = TRUE,
                              cv.lambda,
                              fac.level, ord.fac,
                              seed, nfolds,
                              cv.type, eps = 0.0001){

  beta <- col.genlasso(formula = formula,
                       data = data, pair = pair,
                       cv.lambda = cv.lambda,
                       fac.level = fac.level, ord.fac = ord.fac,
                       seed = seed, nfolds = nfolds,
                       cv.type = cv.type)

  fac.name <- all.vars(formula)[-1]
  # collapsing
  collapse_level <-  Collapse.genlasso(beta = beta, fac.level = fac.level,
                                       ord.fac = ord.fac,
                                       fac.name = fac.name, eps = eps)$collapse.level

  # adjust (if collapse everything, make it binary)
  collapse_level <- lapply(collapse_level, FUN = function(x) if(length(unique(x)) == 1) c(1, rep(2, length(x) - 1)) else x)

  return(collapse_level)

}


fit.after.collapse.gen <- function(formula,
                                   newdata,
                                   collapse_level,
                                   pair=FALSE,
                                   marginal_dist,
                                   marginal_type,
                                   tableAME_base,
                                   difference = FALSE){

  original_level <- lapply(model.frame(formula, data = newdata)[,-1], levels)

  c_data_mar <- prepare_data(formula, data = newdata,
                             marginal_dist = marginal_dist,
                             original_level = original_level,
                             collapse_level = collapse_level)

  collapse_level_name <- lapply(model.frame(formula, c_data_mar$data_new)[,-1], levels)

  # Transform marginal_dist (for internal simplisity) ----------
  marginal_dist_c <- c_data_mar$marginal_dist_new
  marginal_dist_u_list <- list()
  for(z in 1:length(marginal_dist_c)){
    marginal_dist_u_list[[z]] <- data.frame(matrix(NA, ncol=0, nrow=nrow(marginal_dist_c[[z]])))
    marginal_dist_u_list[[z]]$level <- paste(marginal_dist_c[[z]][,1], marginal_dist_c[[z]][,2],sep="")
    marginal_dist_u_list[[z]]$prop  <- marginal_dist_c[[z]][,3]
  }
  marginal_dist_u_base <- marginal_dist_u_list[[1]]


  tableAME <- AME.fit(formula,
                      data = c_data_mar$data_new, pair=pair,
                      marginal_dist = marginal_dist_c,
                      marginal_dist_u_list = marginal_dist_u_list,
                      marginal_dist_u_base = marginal_dist_u_base,
                      marginal_type = marginal_type,
                      difference = difference)

  # Expand
  type_l <- length(unique(tableAME_base$type))
  fac_name <- unique(tableAME_base$factor)
  tableAME_new <- matrix(NA, nrow = 0, ncol = 5)
  for(z in 1:length(fac_name)){
    tableAME_sub  <- tableAME[tableAME$factor == fac_name[z],]
    tableAME_sub$level_num <- match(tableAME_sub$level, collapse_level_name[[fac_name[z]]])

    tableAME_base_sub <- tableAME_base[tableAME_base$factor == fac_name[z],]
    tableAME_base_sub$level_num <- collapse_level[[fac_name[z]]][match(tableAME_base_sub$level, original_level[[fac_name[z]]])]

    tableAME_use <- merge(tableAME_base_sub, tableAME_sub[, c("type", "level_num", "estimate")],
                          by = c("type", "level_num"), all.x = TRUE, all.y = FALSE)
    tableAME_use <- tableAME_use[row.match(tableAME_base_sub[, c("type", "level")], tableAME_use[, c("type", "level")]), ]
    tableAME_use$estimate[tableAME_use$level_num == 1] <- 0
    tableAME_new <- rbind(tableAME_new, tableAME_use)
  }
  tableAME_new$level_num <- NULL

  return(tableAME_new)
}


AME.collapse.genlasso.crossfit.boot <- function(formula,
                                                data,
                                                pair = FALSE,
                                                fac.level, ord.fac,
                                                cv.type = "cv.1Std",
                                                seed = 1234, nfolds = 2,
                                                marginal_dist,
                                                marginal_type,
                                                difference = FALSE,
                                                boot = 100,
                                                tableAME_base,
                                                eps = 0.0001){



  factor_l <- length(all.vars(formula)[-1])
  combMat <- combn(factor_l,2); intNames <- c()
  for(k in 1:ncol(combMat)){
    intNames[k] <- paste(all.vars(formula)[-1][combMat[1,k]], "*", all.vars(formula)[-1][combMat[2,k]], sep = "")
  }
  formula_full <- as.formula(paste(all.vars(formula)[1], "~", paste(intNames, collapse = "+"), sep=""))

  cv.lambda <- col.base.genlasso(formula = formula,
                                 data = data, pair = pair,
                                 fac.level = fac.level, ord.fac = ord.fac)

  # if(pair == TRUE){
  #   remove_id <- names(table(data$cluster))[table(data$cluster) %% 2 != 0]
  #   data <- data[!(data$cluster %in% remove_id), ]
  # }
  data <- data[order(data$cluster), ]

  cat("Bootstrap:")
  fit.mat <- c()
  all_eq <- all(table(data$cluster) == table(data$cluster)[1])
  for(b in 1:boot){

    seed.b <- seed + 1000*b
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
                                     pair = pair,
                                     fac.level = fac.level, ord.fac = ord.fac,
                                     seed = seed.b, nfolds = nfolds,
                                     cv.type = cv.type,
                                     cv.lambda = cv.lambda,
                                     marginal_dist = marginal_dist,
                                     marginal_type = marginal_type,
                                     difference = difference,
                                     tableAME_base = tableAME_base,
                                     eps = eps)

    if(b == 1) fit_0 <- fit
    if(all(fit[,3] == fit_0[,3]) == FALSE) warning("check here")
    fit.mat <- cbind(fit.mat, fit[,4])

    if(b%%10 == 0) cat(paste(b, "..", sep=""))
  }
  estimate <- apply(fit.mat, 1, mean)
  se <- apply(fit.mat, 1, sd)
  low.ci <- apply(fit.mat, 1, function(x) quantile(x, 0.025))
  high.ci <- apply(fit.mat, 1, function(x) quantile(x, 0.975))

  fit[,4] <- estimate
  fit$se <- se
  fit$low.ci <- low.ci
  fit$high.ci <- high.ci
  return(fit)
}

AME.collapse.gen.crossfit <- function(formula,
                                      formula_full,
                                      data,
                                      pair=FALSE,
                                      cv.lambda,
                                      fac.level,
                                      ord.fac,
                                      seed, nfolds, cv.type,
                                      marginal_dist,
                                      marginal_type,
                                      difference = FALSE,
                                      tableAME_base,
                                      eps){

  data <- data[order(data$cluster), ]

  train_id <- sample(unique(data$cluster), size = floor(length(unique(data$cluster))/2), replace = FALSE)
  test_id  <- setdiff(unique(data$cluster), train_id)
  max_cl <- max(table(data$cluster))

  mat_train <- c(t(find.matches(train_id, data$cluster, maxmatch = max_cl)$matches))
  train_which <- mat_train[mat_train > 0]
  mat_test <- c(t(find.matches(test_id, data$cluster, maxmatch = max_cl)$matches))
  test_which <- mat_test[mat_test > 0]

  data_train <- data[train_which, ]
  data_test  <- data[test_which, ]

  # Fit 1
  fit_col_1 <- collapase.fit.gen(formula = formula,
                                 data = data_train, pair = pair,
                                 cv.lambda = cv.lambda,
                                 fac.level = fac.level, ord.fac = ord.fac,
                                 seed = seed, nfolds = nfolds, cv.type = cv.type, eps = eps)

  tableAME_1 <- fit.after.collapse.gen(formula = formula_full,
                                       newdata = data_test,
                                       collapse_level = fit_col_1,
                                       pair = pair,
                                       marginal_dist = marginal_dist,
                                       marginal_type = marginal_type,
                                       tableAME_base = tableAME_base,
                                       difference = difference)

  # Fit 2
  fit_col_2 <- collapase.fit.gen(formula = formula,
                                 data = data_test, pair = pair,
                                 cv.lambda = cv.lambda,
                                 fac.level = fac.level, ord.fac = ord.fac,
                                 seed = seed, nfolds = nfolds, cv.type = cv.type, eps = eps)

  tableAME_2 <- fit.after.collapse.gen(formula = formula_full,
                                       newdata = data_train,
                                       collapse_level = fit_col_2,
                                       pair = pair,
                                       marginal_dist = marginal_dist,
                                       marginal_type = marginal_type,
                                       tableAME_base = tableAME_base,
                                       difference = difference)

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
  table_STD <- AME.fit.STD(formula = formula,
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
        tableAME_add$type <- paste(marginal_type,"-STD",sep="")
        tableAME_add$factor <- rep(uniq_fac[i], length(marginal_type))
        tableAME_add$level <- rep(uniq_level[j], length(marginal_type))
        tableAME_add$estimate <- dif_est
        tableAME_1_main_m <- rbind(tableAME_1_main_m, tableAME_add)
      }
      tableAME_full <- rbind(tableAME_full, tableAME_1_main_m)
    }
  }
  return(tableAME_full)
}

AME.fit <- function(formula_full,
                    data,
                    pair=FALSE,
                    marginal_dist,
                    marginal_dist_u_list,
                    marginal_dist_u_base,
                    marginal_type,
                    difference = FALSE){
  # Differencing ----------
  if(pair==TRUE){
    data0 <- data[order(data$pair_id),]
    side <- rep(c(1,0),times=nrow(data0)/2)
    data1 <- data0[side==1,]
    data2 <- data0[side==0,]
    cluster_original <- data$cluster
    cluster <- data$cluster[side==1]
    X1 <- model.matrix(formula_full, data=data1)[ ,-1]
    X2 <- model.matrix(formula_full, data=data2)[ ,-1]
    X <- cbind(1, X1 - X2)
    y <- model.frame(formula_full,data=data1)[ ,1]
    # base_name <- c("(Intercept)", colnames(X1))
  }else{
    cluster_original <- data$cluster
    X <- model.matrix(formula_full, data=data)
    y <- model.frame(formula_full, data=data)[,1]
    # base_name <- colnames(X)
    side <- NULL
  }

  # Fit the model ----------
  main_lm <- lm(y ~ X - 1)
  coefInt <- coef(main_lm)
  coefInt <- coefInt[is.na(coefInt) == FALSE]
  base_name <- sub("X", "", names(coefInt))
  # base_name <- gsub(" ", "", base_name)
  names(coefInt) <- base_name
  # vcovInt <- vcovCR(main_lm, cluster = as.factor(cluster), type = "CR2")
  # colnames(vcovInt) <- rownames(vcovInt) <- base_name

  # Estimate AMEs ----------
  table_AME <- c()
  for(m in 1:nrow(marginal_dist_u_base)){
    coef_focus <- coefInt[grep(marginal_dist_u_base$level[m], names(coefInt), fixed = T)]
    # vcov_focus <- vcovInt[grep(marginal_dist_u_base$level[m], rownames(vcovInt), fixed = T),
    #                       grep(marginal_dist_u_base$level[m], colnames(vcovInt), fixed = T)]
    if(length(coef_focus) > 0){
      estNames <- gsub(paste(marginal_dist_u_base$level[m], ":", sep = ""), "", names(coef_focus), fixed = T)
      estNames <- gsub(paste(":", names(coef_focus)[1], sep = ""), "", estNames, fixed = T)
      table_AME_m <- c()
      # For each marginal distribution,
      for(z in 1:length(marginal_dist_u_list)){
        marginal_dist_u <- marginal_dist_u_list[[z]]
        # Find weights
        coef_prop <- c(1, as.numeric(as.character(marginal_dist_u[match(estNames, marginal_dist_u[, "level"]), "prop"]))[-1])
        # Compute AMEs
        coef_AME <- sum(coef_focus * coef_prop)
        # se_AME <- sqrt(coef_prop%*%vcov_focus%*%coef_prop)
        AME <- data.frame(matrix(NA, ncol = 0, nrow=1))
        AME$type <- marginal_type[z]
        AME$factor   <- marginal_dist[[z]][m,1]; AME$level <- marginal_dist[[z]][m,2]
        AME$estimate <- coef_AME;
        # AME$se <- se_AME
        table_AME_m <- rbind(table_AME_m, AME)
      }
      if(difference == TRUE){
        for(z in 2:length(marginal_dist_u_list)){
          marginal_dist_u <- marginal_dist_u_list[[z]]
          # Find weights
          coef_prop <- c(1, as.numeric(as.character(marginal_dist_u[match(estNames, marginal_dist_u[, "level"]), "prop"]))[-1])
          coef_prop0 <- c(1, as.numeric(as.character(marginal_dist_u_base[match(estNames, marginal_dist_u_base[, "level"]), "prop"]))[-1])
          # Compute AMEs
          coef_prop_d <- (coef_prop - coef_prop0)
          coef_AME_dif <- sum(coef_focus * coef_prop_d)
          # se_AME_dif <- sqrt(coef_prop_d%*%vcov_focus%*%coef_prop_d)
          AME_dif <- data.frame(matrix(NA, ncol = 0, nrow=1))
          AME_dif$type <- paste(marginal_type[z],"-",marginal_type[1],sep="")
          AME_dif$factor   <- marginal_dist[[z]][m,1]; AME_dif$level <- marginal_dist[[z]][m,2]
          AME_dif$estimate <- coef_AME_dif;
          # AME_dif$se <- se_AME_dif
          table_AME_m <- rbind(table_AME_m, AME_dif)
        }
      }
      table_AME <- rbind(table_AME, table_AME_m)
    }
  }
  colnames(table_AME) <- c("type", "factor", "level", "estimate")
  return(table_AME)
}

AME.fit.STD <- function(formula,
                        data,
                        pair=FALSE,
                        marginal_dist,
                        marginal_dist_u_base){
  # Differencing ----------
  if(pair==TRUE){
    data0 <- data[order(data$pair_id),]
    side <- rep(c(1,0),times=nrow(data0)/2)
    data1 <- data0[side==1,]
    data2 <- data0[side==0,]
    cluster_original <- data$cluster
    cluster <- data$cluster[side==1]
    X1 <- model.matrix(formula, data=data1)[ ,-1]
    X2 <- model.matrix(formula, data=data2)[ ,-1]
    X <- cbind(1, X1 - X2)
    y <- model.frame(formula,data=data1)[ ,1]
    # base_name <- c("(Intercept)", colnames(X1))
  }else{
    cluster_original <- data$cluster
    X <- model.matrix(formula, data=data)
    y <- model.frame(formula, data=data)[,1]
    # base_name <- colnames(X)
    side <- NULL
  }

  # Fit the model ----------
  main_lm <- lm(y ~ X - 1)
  coefInt <- coef(main_lm)
  coefInt <- coefInt[is.na(coefInt) == FALSE]
  base_name <- sub("X", "", names(coefInt))
  # base_name <- gsub(" ", "", base_name)
  names(coefInt) <- base_name
  # vcovInt <- vcovCR(main_lm, cluster = as.factor(cluster), type = "CR2")
  # colnames(vcovInt) <- rownames(vcovInt) <- base_name

  # Estimate AMEs ----------
  table_AME <- c()
  for(m in 1:nrow(marginal_dist_u_base)){
    coef_focus <- coefInt[grep(marginal_dist_u_base$level[m], names(coefInt), fixed = T)]
    # vcov_focus <- vcovInt[grep(marginal_dist_u_base$level[m], rownames(vcovInt), fixed = T),
    #                       grep(marginal_dist_u_base$level[m], colnames(vcovInt), fixed = T)]
    if(length(coef_focus) > 0){
      estNames <- gsub(paste(marginal_dist_u_base$level[m], ":", sep = ""), "", names(coef_focus), fixed = T)
      estNames <- gsub(paste(":", names(coef_focus)[1], sep = ""), "", estNames, fixed = T)

      table_AME_m <- data.frame(matrix(NA, ncol = 0, nrow=1))
      table_AME_m$type <- "STD"
      table_AME_m$factor   <- marginal_dist[[1]][m,1]; table_AME_m$level <- marginal_dist[[1]][m,2]
      table_AME_m$estimate <- coef_focus;
      table_AME <- rbind(table_AME, table_AME_m)
    }
  }
  colnames(table_AME) <- c("type", "factor", "level", "estimate")
  return(table_AME)
}


# GenLasso Collapsing

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

col.genlasso <- function(formula,
                         data, pair = TRUE,
                         cv.lambda,
                         fac.level, ord.fac,
                         seed = 1234, nfolds = 5,
                         cv.type = "cv.1Std"){

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

  # Fit Main
  fit   <- genlasso(X = X, y = y, D = D)

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

    fit.cv   <- genlasso(X = X.cv, y = y.cv, D = D)
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
  cv.sd1 <- max(cv.lambda[cv.error < cv.sd1.value])

  if(cv.type == "cv.1Std") lambda_u <- cv.sd1
  if(cv.type == "cv.min")  lambda_u <- cv.min

  beta_fit <- coef(fit, lambda = lambda_u)$beta

  return(beta_fit)
}

col.base.genlasso <- function(formula,
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

  # Setup D
  D <- makeD(fac.level = fac.level, ord.fac = ord.fac)

  # Fit Main
  fit   <- genlasso(X = X, y = y, D = D)

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

