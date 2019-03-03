#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param pair whether we use the paired-conjoint design
#' @param pair_id id for paired-conjoint design. required when 'pair = TRUE'
#' @export

collapase.fit <- function(formula,
                          family,
                          data, pair = TRUE,
                          nway = 1, collapse = TRUE, collapse.cost,
                          ord.fac){

  fit <- CausalANOVA(formula = formula,
                     family = family,
                     screen = FALSE,
                     data = data, diff = pair, pair.id = data$pair_id,
                     nway = nway, collapse = TRUE, collapse.cost = collapse.cost,
                     ord.fac = ord.fac, verbose = FALSE)

  # collapsing
  collapse_level <- Collapsing(fit)$collapse.level

  # adjust (if collapse everything, make it binary)
  collapse_level <- lapply(collapse_level, FUN = function(x) if(length(unique(x)) == 1) c(1, rep(2, length(x) - 1)) else x)

  return(collapse_level)

}


fit.after.collapse <- function(formula,
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


AME.collapse.crossfit.boot <- function(formula,
                                       data,
                                       pair = FALSE,
                                       ord.fac,
                                       nway = 1,
                                       cv.collapse.cost = c(0.1, 0.3, 0.5, 0.7),
                                       cv.type = "cv.1Std",
                                       marginal_dist,
                                       marginal_type,
                                       family,
                                       difference = FALSE,
                                       boot = 100,
                                       tableAME_base,
                                       seed){



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

  data <- data[order(data$cluster), ]

  cat("Bootstrap:")
  fit.mat <- c()
  max_cl <- max(table(data$cluster))

  for(b in 1:boot){

    seed.b <- seed + 1000*b
    set.seed(seed.b)
    boot_id   <- sample(unique(data$cluster), size = length(unique(data$cluster)), replace = TRUE)
    mat_base  <- find.matches(boot_id, data$cluster, maxmatch = max_cl)$matches
    mat_boot <- c(t(mat_base))
    boot_which <- mat_boot[mat_boot > 0]
    data_boot <- data[boot_which, ]
    data_boot$cluster <- rep(seq(1:nrow(mat_base)), times = apply(mat_base, 1, function(x) sum(x > 0)))
    data_boot$pair_id <- paste0(data_boot$cluster, data_boot$pair_id)

    fit <- AME.collapse.crossfit(formula = formula,
                                 formula_full = formula_full,
                                 data = data_boot,
                                 pair = pair,
                                 nway = nway, ord.fac = ord.fac,
                                 collapse.cost = collapse.cost,
                                 marginal_dist = marginal_dist,
                                 marginal_type = marginal_type,
                                 difference = difference,
                                 family = family,
                                 tableAME_base = tableAME_base)

    if(b == 1) fit_0 <- fit
    if(all(fit[,3] == fit_0[,3]) == FALSE) warning("check here")
    fit.mat <- cbind(fit.mat, fit[,4])

    if(b%%10 == 0) cat(paste(b, "..", sep=""))
  }
  estimate <- apply(fit.mat, 1, mean) # Added na.rm argument Brandon 2/9/2019
  se <- apply(fit.mat, 1, sd) # Added na.rm argument Brandon 2/9/2019

  fit[,4] <- estimate
  fit$se <- se
  return(fit)
}

AME.collapse.crossfit <- function(formula,
                                  formula_full,
                                  data,
                                  pair=FALSE,
                                  nway = 1,
                                  collapse.cost,
                                  ord.fac,
                                  marginal_dist,
                                  marginal_type,
                                  difference = FALSE,
                                  family,
                                  tableAME_base){

  # training
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
  fit_col_1 <- collapase.fit(formula = formula,
                             family = family,
                             data = data_train, pair = pair,
                             nway = nway, collapse = TRUE, collapse.cost = collapse.cost,
                             ord.fac = ord.fac)

  tableAME_1 <- fit.after.collapse(formula = formula_full,
                                   newdata = data_test,
                                   collapse_level = fit_col_1,
                                   pair = pair,
                                   marginal_dist = marginal_dist,
                                   marginal_type = marginal_type,
                                   tableAME_base = tableAME_base,
                                   difference = difference)

  # Fit 2
  fit_col_2 <- collapase.fit(formula = formula,
                             family = family,
                             data = data_test, pair = pair,
                             nway = nway, collapse = TRUE, collapse.cost = collapse.cost,
                             ord.fac = ord.fac)

  tableAME_2 <- fit.after.collapse(formula = formula_full,
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
