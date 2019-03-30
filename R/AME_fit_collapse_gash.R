#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param pair whether we use the paired-conjoint design
#' @param pair_id id for paired-conjoint design. required when 'pair = TRUE'
#' @importFrom FindIt cv.CausalANOVA CausalANOVA
#' @export

collapse.fit <- function(formula,
                          family,
                          data, pair = TRUE,
                          nway = 1, collapse = TRUE, collapse.cost,
                          fac.level, ord.fac){

  fit <- CausalANOVA(formula = formula,
                     family = family,
                     screen = FALSE,
                     data = data, diff = pair, pair.id = data$pair_id,
                     nway = nway, collapse = TRUE, collapse.cost = collapse.cost,
                     ord.fac = ord.fac, verbose = FALSE)

  # collapsing
  fac.name <- all.vars(formula)[-1]
  collapse_level <-  Collapse.gash(fit, fac.level = fac.level,
                                   ord.fac = ord.fac,
                                   fac.name = fac.name, eps = 0.001)$collapse.level

  # adjust (if collapse everything, make it binary)
  collapse_level <- lapply(collapse_level, FUN = function(x) if(length(unique(x)) == 1) c(1, rep(2, length(x) - 1)) else x)

  return(collapse_level)

}


fit.after.collapse <- function(formula_full,
                               newdata,
                               collapse_level,
                               pair=FALSE,
                               marginal_dist,
                               marginal_type,
                               tableAME_base,
                               difference = FALSE,
                               original_level){

  c_data_mar <- prepare_data(formula_full, data = newdata,
                             marginal_dist = marginal_dist,
                             original_level = original_level,
                             collapse_level = collapse_level)

  collapse_level_name <- lapply(model.frame(formula_full, c_data_mar$data_new)[,-1], levels)

  # Transform marginal_dist (for internal simplisity) ----------
  marginal_dist_c <- c_data_mar$marginal_dist_new
  marginal_dist_u_list <- list()
  for(z in 1:length(marginal_dist_c)){
    marginal_dist_u_list[[z]] <- data.frame(matrix(NA, ncol=0, nrow=nrow(marginal_dist_c[[z]])))
    marginal_dist_u_list[[z]]$level <- paste(marginal_dist_c[[z]][,1], marginal_dist_c[[z]][,2],sep="")
    marginal_dist_u_list[[z]]$prop  <- marginal_dist_c[[z]][,3]
  }
  marginal_dist_u_base <- marginal_dist_u_list[[1]]


  fitAME <- AME.fit(formula_full,
                    data = c_data_mar$data_new, pair=pair,
                    marginal_dist = marginal_dist_c,
                    marginal_dist_u_list = marginal_dist_u_list,
                    marginal_dist_u_base = marginal_dist_u_base,
                    marginal_type = marginal_type,
                    difference = difference)

  tableAME <- fitAME$table_AME
  coefAME  <- fitAME$coef
  ind_b    <- fitAME$ind_b

  # Expand Coefficients
  n_fac <- length(all.vars(formula_full)) - 1
  # For main effects
  coefAME_main <- coefAME[1]
  for(z in 1:n_fac){
    coefAME_sub  <- coefAME[ind_b == z]
    coefAME_m0 <- c(rep(0, times = sum(collapse_level[[z]] == 1) - 1), coefAME_sub[collapse_level[[z]] - 1])
    coefAME_main <- c(coefAME_main, coefAME_m0)
  }
  # For Interaction effects
  combMat <- combn(n_fac, 2)
  coefAME_int <- c()
  for(z in 1:ncol(combMat)){
    coefAME_sub <- coefAME[ind_b == (z + n_fac)]
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

  out <- list("tableAME_new" = tableAME_new, "coef" = coefAME_long)

  return(out)
}

crossGashFitPar <- function(x,
                        formula,
                        formula_full,
                        data,
                        pair,
                        nway,
                        fac.level, ord.fac,
                        collapse.cost,
                        marginal_dist,
                        marginal_type,
                        difference,
                        family,
                        tableAME_base,
                        original_level,
                        all_eq){

  seed.b <- 1000*x
  set.seed(seed.b)
  boot_id <- sample(unique(data$cluster), size = length(unique(data$cluster)), replace=TRUE)
  # create bootstap sample with sapply
  boot_which <- sapply(boot_id, function(x) which(data$cluster == x))
  if(all_eq == TRUE){new_boot_id <- rep(seq(1:length(boot_id)), each = table(data$cluster)[1])
  }else{new_boot_id <- rep(seq(1:length(boot_id)), times = unlist(lapply(boot_which, length)))}
  data_boot <- data[unlist(boot_which),]
  data_boot$cluster <- new_boot_id
  data_boot$pair_id <- paste0(data_boot$cluster, data_boot$pair_id)

  fit <- AME.collapse.crossfit(formula = formula,
                               formula_full = formula_full,
                               data = data_boot,
                               pair = pair,
                               nway = nway,
                               fac.level = fac.level, ord.fac = ord.fac,
                               collapse.cost = collapse.cost,
                               marginal_dist = marginal_dist,
                               marginal_type = marginal_type,
                               difference = difference,
                               family = family,
                               tableAME_base = tableAME_base,
                               original_level = original_level)

  return(fit)
}


#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param pair whether we use the paired-conjoint design
#' @param pair_id id for paired-conjoint design. required when 'pair = TRUE'
#' @importFrom FindIt cv.CausalANOVA CausalANOVA
#' @import arm
#' @export

AME.collapse.crossfit.boot <- function(formula,
                                       data,
                                       pair = FALSE,
                                       fac.level, ord.fac,
                                       nway = 1,
                                       cv.collapse.cost = c(0.1, 0.3, 0.5, 0.7),
                                       cv.type = "cv.1Std",
                                       marginal_dist,
                                       marginal_type,
                                       family,
                                       difference = FALSE,
                                       boot = 100,
                                       tableAME_base,
                                       coefAME_base_l,
                                       numCores){



  factor_l <- length(all.vars(formula)[-1])
  combMat <- combn(factor_l,2); intNames <- c()
  for(k in 1:ncol(combMat)){
    intNames[k] <- paste(all.vars(formula)[-1][combMat[1,k]], "*", all.vars(formula)[-1][combMat[2,k]], sep = "")
  }
  formula_full <- as.formula(paste(all.vars(formula)[1], "~", paste(intNames, collapse = "+"), sep=""))

  cat("Cross-Validation:\n")
  cv.fit <- cv.CausalANOVA(formula = formula,
                           data = data, diff = pair,
                           pair.id = data$pair_id,
                           cluster = data$cluster,
                           ord.fac = ord.fac,
                           family = family,
                           cv.collapse.cost = cv.collapse.cost,
                           nway = nway, verbose = FALSE)

  if(cv.type == "cv.1Std") collapse.cost <- cv.fit$cv.1Std
  if(cv.type == "cv.min")  collapse.cost <- cv.fit$cv.min

  cat(paste("Selected cost parameter:", collapse.cost, sep = ""))

  data <- data[order(data$cluster), ]

  original_level <- lapply(model.frame(formula, data = data)[,-1], levels)

  cat("\nBootstrap:\n")
  fit.mat <- c()
  coef.mat <- matrix(NA, nrow = boot, ncol = coefAME_base_l)
  all_eq <- all(table(data$cluster) == table(data$cluster)[1])

  if(Sys.info()[['sysname']] == 'Windows') {

    if (numCores == 1){
      fit_boot <- pblapply(1:boot, function(x) crossGashFitPar(x,
                                                               formula = formula,
                                                               formula_full = formula_full,
                                                               data = data,
                                                               pair = pair,
                                                               nway = nway,
                                                               fac.level = fac.level, ord.fac = ord.fac,
                                                               collapse.cost = collapse.cost,
                                                               marginal_dist = marginal_dist,
                                                               marginal_type = marginal_type,
                                                               difference = difference,
                                                               family = family,
                                                               tableAME_base = tableAME_base,
                                                               original_level = original_level,
                                                               all_eq = all_eq))
    }else {

      cl <- makeCluster(numCores)
      registerDoSNOW(cl)
      pb <- txtProgressBar(max = boot, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      #
      #     cl <- makeCluster(numCores)
      #     registerDoParallel(cl)
      #     on.exit(stopCluster(cl))

      fit_boot <- foreach(i = 1:boot,
                          .export = c("prepare_data"),
                          .packages = c("FindIt", "prodlim"),
                          .options.snow = opts) %dopar% {
                            crossGashFitPar(i,
                                            formula = formula,
                                            formula_full = formula_full,
                                            data = data,
                                            pair = pair,
                                            nway = nway,
                                            fac.level = fac.level, ord.fac = ord.fac,
                                            collapse.cost = collapse.cost,
                                            marginal_dist = marginal_dist,
                                            marginal_type = marginal_type,
                                            difference = difference,
                                            family = family,
                                            tableAME_base = tableAME_base,
                                            original_level = original_level,
                                            all_eq = all_eq)
                          }
      close(pb)
      stopCluster(cl)
    }

  }else{

    # -----
    # Mac
    # -----

  fit_boot <- pbmclapply(seq(1:boot), function(x) crossGashFitPar(x,
                                                                formula = formula,
                                                                formula_full = formula_full,
                                                                data = data,
                                                                pair = pair,
                                                                nway = nway,
                                                                fac.level = fac.level, ord.fac = ord.fac,
                                                                collapse.cost = collapse.cost,
                                                                marginal_dist = marginal_dist,
                                                                marginal_type = marginal_type,
                                                                difference = difference,
                                                                family = family,
                                                                tableAME_base = tableAME_base,
                                                                original_level = original_level,
                                                                all_eq = all_eq),
                       mc.cores = numCores)
  }

  for(b in 1:boot){
    coef.mat[b, 1:coefAME_base_l] <- fit_boot[[b]]$coef
    fit.mat <- cbind(fit.mat, fit_boot[[b]]$tableAME_full[,4])
  }

  # for(b in 1:boot){
  #
  #   # seed.b <- seed + 1000*b
  #   # set.seed(seed.b)
  #   boot_id <- sample(unique(data$cluster), size = length(unique(data$cluster)), replace=TRUE)
  #   # create bootstap sample with sapply
  #   boot_which <- sapply(boot_id, function(x) which(data$cluster == x))
  #   if(all_eq == TRUE){new_boot_id <- rep(seq(1:length(boot_id)), each = table(data$cluster)[1])
  #   }else{new_boot_id <- rep(seq(1:length(boot_id)), times = unlist(lapply(boot_which, length)))}
  #   data_boot <- data[unlist(boot_which),]
  #   data_boot$cluster <- new_boot_id
  #   data_boot$pair_id <- paste0(data_boot$cluster, data_boot$pair_id)
  #
  #   fitC <- AME.collapse.crossfit(formula = formula,
  #                                 formula_full = formula_full,
  #                                 data = data_boot,
  #                                 pair = pair,
  #                                 nway = nway,
  #                                 fac.level = fac.level, ord.fac = ord.fac,
  #                                 collapse.cost = collapse.cost,
  #                                 marginal_dist = marginal_dist,
  #                                 marginal_type = marginal_type,
  #                                 difference = difference,
  #                                 family = family,
  #                                 tableAME_base = tableAME_base,
  #                                 original_level = original_level)
  #
  #   # Store coefficients
  #   coef.mat[b, 1:coefAME_base_l] <- fitC$coef
  #   fit <- fitC$tableAME_full
  #
  #   if(b == 1) fit_0 <- fit
  #   if(all(fit[,3] == fit_0[,3]) == FALSE) warning("check here")
  #   fit.mat <- cbind(fit.mat, fit[,4])
  #
  #   if(b%%10 == 0) cat(paste(b, "..", sep=""))
  # }

  fit <- fit_boot[[1]]$tableAME_full
  estimate <- apply(fit.mat, 1, mean)
  se <- apply(fit.mat, 1, sd)

  fit[,4] <- estimate
  fit$se <- se

  out <- list("fit" = fit, "fit.mat" = fit.mat, "coef.mat" = coef.mat)
  return(out)
}

AME.collapse.crossfit <- function(formula,
                                  formula_full,
                                  data,
                                  pair=FALSE,
                                  nway = 1,
                                  collapse.cost,
                                  fac.level, ord.fac,
                                  marginal_dist,
                                  marginal_type,
                                  difference = FALSE,
                                  family,
                                  tableAME_base,
                                  original_level){

  # training
  data <- data[order(data$cluster), ]

  train_id <- sample(unique(data$cluster), size = floor(length(unique(data$cluster))/2), replace = FALSE)
  test_id  <- setdiff(unique(data$cluster), train_id)

  train_which <- unlist(sapply(train_id, function(x) which(data$cluster == x)))
  test_which  <- unlist(sapply(test_id, function(x) which(data$cluster == x)))

  data_train <- data[train_which, ]
  data_test  <- data[test_which, ]

  # Fit 1
  fit_col_1 <- tryCatch({collapse.fit(formula = formula,
                                       family = family,
                                       data = data_train, pair = pair,
                                       nway = nway, collapse = TRUE, collapse.cost = collapse.cost,
                                       fac.level = fac.level, ord.fac = ord.fac)
  }, error=function(e){
    cat(paste("warning: no collapsing"))
    fit_col_1 <- lapply(original_level, function(x) seq(1:length(x)))
    return(fit_col_1)
  })
#
#   fit_col_1 <- collapse.fit(formula = formula,
#                              family = family,
#                              data = data_train, pair = pair,
#                              nway = nway, collapse = TRUE, collapse.cost = collapse.cost,
#                              fac.level = fac.level, ord.fac = ord.fac)

  fitAME_1 <- fit.after.collapse(formula_full = formula_full,
                                 newdata = data_test,
                                 collapse_level = fit_col_1,
                                 pair = pair,
                                 marginal_dist = marginal_dist,
                                 marginal_type = marginal_type,
                                 tableAME_base = tableAME_base,
                                 difference = difference,
                                 original_level = original_level)

  tableAME_1 <- fitAME_1$tableAME_new
  coefAME_1  <- fitAME_1$coef

  # Fit 2
  # fit_col_2 <- collapse.fit(formula = formula,
  #                            family = family,
  #                            data = data_test, pair = pair,
  #                            nway = nway, collapse = TRUE, collapse.cost = collapse.cost,
  #                            fac.level = fac.level, ord.fac = ord.fac)

  fit_col_2 <- tryCatch({collapse.fit(formula = formula,
                                       family = family,
                                       data = data_test, pair = pair,
                                       nway = nway, collapse = TRUE, collapse.cost = collapse.cost,
                                       fac.level = fac.level, ord.fac = ord.fac)
  }, error=function(e){
    cat(paste("warning: no collapsing"))
    fit_col_2 <- lapply(original_level, function(x) seq(1:length(x)))
    return(fit_col_2)
  })

  fitAME_2 <- fit.after.collapse(formula_full = formula_full,
                                 newdata = data_train,
                                 collapse_level = fit_col_2,
                                 pair = pair,
                                 marginal_dist = marginal_dist,
                                 marginal_type = marginal_type,
                                 tableAME_base = tableAME_base,
                                 difference = difference,
                                 original_level = original_level)

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
  out <- list("tableAME_full" = tableAME_full, "coef" = coefAME_full)
  return(out)
}

Collapse.gash <- function(fitGash, fac.level, ord.fac, fac.name, eps = 0.001){

  beta <- c(0.5, unlist(lapply(fitGash$AME, function(x) (x - x[1])[-1])))

  n.fac <- length(fac.level)
  # beta: [1] intercept from [2]
  start_ind <- as.numeric(c(2, c(cumsum(fac.level -1) + 2)[-length(fac.level)]))
  end_ind   <- as.numeric(1 + cumsum(fac.level -1))

  ## Do Collapsing for each factor
  Collapse <- list()
  for(z in 1:length(fac.level)){
    ## First Order

    ## type.ind <- z
    MainD1  <- D1function.gash(nlevel = fac.level[z], ord = ord.fac[z])
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

D1function.gash <- function(nlevel,ord = FALSE){
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
