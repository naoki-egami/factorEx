#' Estimating the population AMCE using the model-based approach
#' @param formula Formula
#' @param formula_three Formula for three-way interactions
#' @param data Data
#' @param reg  TRUE (regularization) or FALSE (no regularization)
#' @param ord.fac Whether we assume each factor is ordered. When not specified, we assume all of them are ordered
#' @param pair Whether we use a paired-choice conjoint design
#' @param cross_int Include interactions across profiles
#' @param cluster unique identifiers for cluster
#' @param marginal_dist Marginal distributions of profiles to be used
#' @param marginal_type Names of marginal distributions
#' @param difference Whether we compute the difference between different pAMCEs
#' @param cv.type (when reg = TRUE)  `cv.1Std`` (stronger) or `cv.min` (weaker).
#' @param boot The number of bootstrap samples
#' @param seed Seed for bootstrap
#' @importFrom prodlim row.match
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom pbmcapply pbmclapply
#' @importFrom pbapply pblapply
#' @import parallel
#' @import arm
#' @import genlasso
#' @export

pAMCE <- function(formula,
                  formula_three = NULL,
                  data,
                  reg = TRUE,
                  ord.fac,
                  pair = FALSE, pair_id = NULL, cross_int = TRUE,
                  cluster = NULL,
                  target_dist,
                  target_type,
                  difference = FALSE,
                  cv.type = "cv.1Std", nfolds = 5,
                  boot = 100,
                  seed = 1234, numCores = NULL){

  ###########
  ## Check ##
  ###########

  if(missing(pair_id) == TRUE) pair_id <- NULL
  if(missing(cluster) == TRUE) cluster <- seq(1:nrow(data))
  if(missing(target_type) == TRUE){
    target_type <- rep("marginal", length(target_dist))
  }

  if(all(target_type %in% c("marginal", "joint", "target_data")) == FALSE){
    warning(" 'target_type' should be 'marginal', 'joint' or 'target_data' ")
  }

  if(length(target_type) != length(target_dist)){
    stop(" length of 'target_type' should be the same as length of 'target_dist' ")
  }

  if(is.null(names(target_dist)) == TRUE){
    target_name <- paste("dist_", seq(1:length(target_dist)), sep = "")
  }else{
    target_name <- names(target_dist)
  }

  if(is.null(numCores) == FALSE){
    if(numCores >= detectCores()) numCores <- detectCores() - 1
  }
  if(is.null(numCores)) numCores <- detectCores() - 1

  ###########
  ## Check ##
  ###########
  if(any(is.na(data))==TRUE){
    stop("Remove NA before using this function.")
  }
  if(pair==TRUE & is.null(pair_id)==TRUE){
    stop("When 'pair=TRUE', specify 'pair_id'.")
  }
  if(pair==TRUE & all(table(pair_id)==2)==FALSE){
    stop("When 'pair=TRUE', each of 'pair_id' should have two observations")
  }
  if(is.null(pair_id) == FALSE) pair <- TRUE

  if(pair == FALSE){
    cross_int <- FALSE
  }

  # if(is.null(pair_var) == FALSE){
  #   if(all(is.element(pair_var, all.vars(formula))) == FALSE){
  #     stop(" 'pair_var' should be variables listed in 'formula' ")}
  #   if(pair == FALSE){
  #     stop(" 'pair_var' is ignored when 'pair=FALSE' ")}
  # }

  if(difference==TRUE & length(target_dist) < 2){
    stop("if 'difference = TRUE', target_dist should contain more than one distribution.")
  }
  if(class(target_dist) != "list"){
    target_dist <- list(target_dist)
  }
  if(is.list(target_dist)==FALSE){
    stop("target_dist should be 'list'.")
  }

  ## Check each distribution
  for(i in 1:length(target_dist)){
    target_dist[[i]] <- checkDist(target_dist[[i]], type = target_type[i], formula = formula, data = data)
  }

  ## Create Marginals and 2d-Joints from target_dist
  marginal_dist <- list()
  for(i in 1:length(target_dist)){
    if(target_type[i]  == "marginal"){marginal_dist[[i]] <- target_dist[[i]]}
    if(target_type[i]  == "joint"){marginal_dist[[i]] <- Joint2Marginal(target_dist[[i]])}
    if(target_type[i]  == "target_data"){marginal_dist[[i]] <- createDist(formula = formula,
                                                                          target_data = target_dist[[i]],
                                                                          exp_data = data, type  = "marginal")}
    marginal_dist[[i]] <- checkDist(marginal_dist[[i]], type = "marginal", formula = formula, data = data)
  }
  if(is.null(formula_three) == TRUE){joint_dist <- NULL
  } else if(is.null(formula_three) == FALSE){
    joint_dist <- list()
    for(i in 1:length(target_dist)){
      if(target_type[i]  == "marginal"){joint_dist[[i]] <- Marginal2Joint(target_dist[[i]])}
      if(target_type[i]  == "joint"){joint_dist[[i]] <- target_dist[[i]]}
      if(target_type[i]  == "target_data"){joint_dist[[i]] <- createDist(formula = formula,
                                                                         target_data = target_dist[[i]],
                                                                         exp_data = data, type  = "joint")}
      joint_dist[[i]] <- checkDist(joint_dist[[i]], type = "joint", formula = formula, data = data)
    }
  }

  # make marginal_dist to data.frame
  marginal_dist_internal <- marginal_dist
  marginal_dist <- lapply(marginal_dist, list2data)
  names(marginal_dist) <- target_name

  # Check Marginal Distributions
  marginal_name_check <- lapply(marginal_dist, colnames)
  marginal_name_check_all <- all(unlist(lapply(marginal_name_check,
                                               function(x) all(x == c("factor", "levels", "prop")))))
  if(marginal_name_check_all == FALSE){
    stop(" 'colnames' of 'marginal_dist' should be c('factor', 'levels', 'prop') ")
  }
  if(length(target_type) > 1){
    for(z in 2:length(target_type)){
      if(all(marginal_dist[[1]][,1] == marginal_dist[[z]][,1]) == FALSE) stop("marginal_dist should have the same order for factor.")
      if(all(marginal_dist[[1]][,2] == marginal_dist[[z]][,2]) == FALSE) stop("marginal_dist should have the same order for levels.")
    }
  }
  # make them as character
  for(z in 1:length(marginal_dist)){
    marginal_dist[[z]]$factor <- as.character(marginal_dist[[z]]$factor)
    marginal_dist[[z]]$levels <- as.character(marginal_dist[[z]]$levels)
  }

  ## Check Baselines
  factor_l <- length(all.vars(formula)[-1])
  combMat <- combn(factor_l,2); intNames <- c()
  for(k in 1:ncol(combMat)){
    intNames[k] <- paste(all.vars(formula)[-1][combMat[1,k]], "*", all.vars(formula)[-1][combMat[2,k]], sep = "")
  }
  formula_full <- as.formula(paste(all.vars(formula)[1], "~", paste(intNames, collapse = "+"), sep=""))

  check_levels <- model.matrixBayes(formula_full, data = data)
  original_level <- lapply(model.frame(formula, data = data)[,-1], FUN = function(x) levels(x))
  mar_length <- length(unlist(original_level))
  rest_name <- colnames(check_levels)[(1 + mar_length):ncol(check_levels)][apply(check_levels[, (1 + mar_length):ncol(check_levels)],
                                                                                 2, function(x) mean(x)) == 0]

  rest_level <- unlist(strsplit(rest_name, ":"))
  baseline <- baseline_orig <- lapply(model.frame(formula, data=data)[,-1], FUN = function(x) levels(x)[1])
  basenames <- paste(all.vars(formula)[-1], unlist(baseline), sep = "")
  rest_base <- basenames[is.element(basenames, rest_level)]
  rest_fac  <- all.vars(formula)[-1][is.element(basenames, rest_level)]
  if(any(is.element(basenames, rest_level))){
    wa1 <- paste("The following ", length(rest_name), " pairs have no observation:", sep = "")
    wa2 <- paste(paste(paste(" (", seq(1:length(rest_name)), ") ", sep = ""), rest_name, sep = ""), collapse = ", ")
    wa3 <- "To incorporate the restriction, for each factor, select a baseline category that has no restriction."
    wa4 <- paste("Currently, ", length(rest_base), " baselines (",
                 paste(rest_base, collapse = ","), ") have restrictions.", sep = "")
    wa5 <- paste("Change baselines for ", paste(rest_fac, collapse = " and "), " using relevel().", sep = "")
    stop(paste("\n", wa1, wa2, wa3, wa4, wa5, sep = "\n"))
  }

  # ############
  # Renaming
  # ############
  {
    # Rename formula
    formula_orig <- formula
    fac_size <- length(all.vars(formula)[-1])
    fac_name_orig <- all.vars(formula)[-1]
    fac_name <- paste("f_", seq(1:fac_size), "_f", sep = "")
    formula  <- as.formula(paste("Y ~", paste(fac_name, collapse = "+"), sep = ""))
    rename_fac <- cbind(all.vars(formula_orig), c("Y", fac_name))
    colnames(rename_fac) <- c("original", "internal")

    if(is.null(formula_three) == FALSE){
      formula_three_c <- as.character(formula_three)[2]
      for(i in 2:nrow(rename_fac)){
        formula_three_c <- gsub(rename_fac[i,1], rename_fac[i,2], formula_three_c)
      }
    }else{
      formula_three_c <- NULL
    }

    # Rename levels
    original_level <- lapply(model.frame(formula_orig, data = data)[,-1], FUN = function(x) levels(x))
    internal_level <- list()
    for(i in 1:length(original_level)){
      internal_level[[i]] <- paste("x_", i, "_", seq(1:length(original_level[[i]])), "_x", sep = "")
    }
    names(internal_level) <- fac_name
    fac_name_0   <- rep(names(internal_level), unlist(lapply(internal_level, length)))
    rename_level <- cbind(fac_name_0, unlist(internal_level), unlist(original_level))

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

    # Rename joint_dist
    joint_dist_orig <- joint_dist
    if(is.null(joint_dist) == FALSE){
      for(z in 1:length(joint_dist)){
        for(j in 1:length(joint_dist[[z]])){
          joint_dist[[z]][[j]]$factor_1 <- rename_fac[match(joint_dist[[z]][[j]]$factor_1, rename_fac[, 1]), 2]
          joint_dist[[z]][[j]]$factor_2 <- rename_fac[match(joint_dist[[z]][[j]]$factor_2, rename_fac[, 1]), 2]
          rename_level_1 <- rename_level[rename_level[,1]%in% joint_dist[[z]][[j]]$factor_1, ]
          rename_level_2 <- rename_level[rename_level[,1]%in% joint_dist[[z]][[j]]$factor_2, ]
          joint_dist[[z]][[j]]$levels_1 <- rename_level_1[match(joint_dist[[z]][[j]]$levels_1, rename_level_1[, 3]), 2]
          joint_dist[[z]][[j]]$levels_2 <- rename_level_2[match(joint_dist[[z]][[j]]$levels_2, rename_level_2[, 3]), 2]
        }
      }
    }

  }

  if(reg == FALSE){
    out <-  AME_estimate(formula = formula,
                         data = data,
                         pair = pair, pair_id = pair_id, cross_int = cross_int,
                         cluster = cluster,
                         marginal_dist = marginal_dist,
                         marginal_type = target_name,
                         joint_dist = joint_dist,
                         boot = boot,
                         difference = difference, formula_three_c = formula_three_c)
  }else if(reg == TRUE){
    if(missing(ord.fac)) ord.fac <- rep(TRUE, (length(all.vars(formula)) - 1))
    out <- AME_estimate_collapse_genlasso(formula = formula,
                                          data = data,
                                          ord.fac = ord.fac,
                                          pair = pair, pair_id = pair_id, cross_int = cross_int,
                                          cluster = cluster,
                                          marginal_dist = marginal_dist,
                                          marginal_type = target_name,
                                          joint_dist = joint_dist,
                                          formula_three_c = formula_three_c,
                                          difference = difference,
                                          cv.type = cv.type,
                                          nfolds = nfolds,
                                          boot = boot,
                                          eps = 0.0001,
                                          numCores = numCores,
                                          seed = seed)
  }

  ## Approximate F-test
  if(is.null(formula_three_c) == FALSE){
    if(reg == TRUE) coef_f <- apply(out$boot_coef, 2, mean)
    if(reg == FALSE) coef_f <- out$coef
    data_u <- out$input$data

    Ftest <- Fthree(formula = formula,
                    formula_three_c = formula_three_c,
                    data = data_u,
                    pair = pair, cross_int = cross_int,
                    coef_f = coef_f)

    out$Ftest <- as.numeric(Ftest)
  }



  # ##############
  # Name back
  # ##############
  {
    out$input$formula <- formula_orig
    out$input$data <- data_orig
    out$input$marginal_dist <- marginal_dist_orig
    out$baseline <- baseline_orig

    marginal_dist_u_list <- list()
    for(z in 1:length(marginal_dist_orig)){
      marginal_dist_u_list[[z]] <- data.frame(matrix(NA, ncol=0, nrow=nrow(marginal_dist_orig[[z]])))
      marginal_dist_u_list[[z]]$level <- paste(marginal_dist_orig[[z]][,1], marginal_dist_orig[[z]][,2],sep="")
      marginal_dist_u_list[[z]]$prop  <- marginal_dist_orig[[z]][,3]
    }
    marginal_dist_u_base <- marginal_dist_u_list[[1]]
    out$input$marginal_dist_u_list <- marginal_dist_u_list
    out$input$marginal_dist_u_base <- marginal_dist_u_base

    # AME
    for(i in 1:length(out$AME)){
      match_level <- match(out$AME[[i]]$level, internal_level[[out$AME[[i]]$factor[1]]])
      out$AME[[i]]$factor <- rename_fac[,"original"][match(out$AME[[i]]$factor[1], rename_fac[,"internal"])]
      orignal_use <- original_level[[out$AME[[i]]$factor[1]]]
      out$AME[[i]]$level <- orignal_use[match_level]
    }
    names(out$AME) <- rename_fac[,"original"][match(names(out$AME), rename_fac[,"internal"])]
    # coefficients
    for(i in 1:fac_size){
      colnames(out$boot_coef) <- gsub(rename_fac[(i+1), "internal"], rename_fac[(i+1), "original"], colnames(out$boot_coef))
      for(j in 1:length(internal_level[[i]])){
        colnames(out$boot_coef) <- gsub(internal_level[[i]][j], original_level[[i]][j], colnames(out$boot_coef))
      }
    }

  }

  return(out)
}
