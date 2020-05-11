#' Estimating the population AMCE using a model-based approach
#' @param formula Formula
#' @param formula_three Formula for three-way interactions (optional)
#' @param data Data
#' @param reg  TRUE (regularization) or FALSE (no regularization). Default is TRUE
#' @param ord_fac Whether we assume each factor is ordered. When not specified, we assume all of them are ordered
#' @param pair Whether we use a paired-choice conjoint design
#' @param pair_id Unique identifiers for pairs in the paired-choice conjoint design  (optional)
#' @param cluster_id Unique identifiers for computing cluster standard errors (optional).
#' @param cross_int Include interactions across profiles. Default is FALSE.
#' @param target_dist Target profile distributions to be used. This argument should be `list`
#' @param target_type Types of target profile distributions. `marginal` or `target_data`. See Examples for details.
#' @param difference Whether we compute the differences between the multiple pAMCEs. Default is FALSE.
#' @param cv_type (optimal only when `reg = TRUE``)  `cv.1Std`` (stronger regularization; default) or `cv.min` (weaker regularization).
#' @param nfolds Number of cross validation folds. Default is 5.
#' @param boot The number of bootstrap samples.
#' @param seed Seed for bootstrap.
#' @param numCores Number of cores to be used for parallel computing. If not specified, detect the number of available cores internally.
#' @import parallel
#' @import stringr
#' @import arm
#' @import genlasso
#' @importFrom prodlim row.match
#' @importFrom igraph graph_from_adjacency_matrix components
#' @importFrom pbmcapply pbmclapply
#' @importFrom pbapply pblapply
#' @importFrom sandwich sandwich estfun
#' @importFrom mvtnorm rmvnorm
#' @importFrom grDevices dev.off palette
#' @importFrom graphics Axis abline legend par plot segments
#' @importFrom stats as.formula coef density lm model.frame model.matrix pf pnorm quantile sd terms
#' @importFrom utils combn setTxtProgressBar txtProgressBar globalVariables
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach "%dopar%" "%do%" foreach
#' @return \code{model_pAMCE} returns an object of \code{pAMCE} class.
#'  \itemize{
#'    \item \code{AMCE}: Estimates of the pAMCE for all factors.
#'    \item \code{boot_AMCE}: Estimates of the pAMCE for all factors in each bootstrap sample.
#'    \item \code{boot_coef}: Estimates of coefficients for the linear probability model in each bootstrap sample.
#'    \item \code{approach}: "model_based"
#'    \item \code{input}: Input into the function.
#'    \item \code{...}: Values for internal use.
#'  }
#' @examples
#'   # Small example
#'   target_dist_marginal <- OnoBurden$target_dist_marginal
#'   OnoBurden_data <- OnoBurden$OnoBurden_data
#'   OnoBurden_data_small <- OnoBurden_data[1:300, ]
#'   target_dist_marginal_small <- target_dist_marginal[c("gender", "race")]
#'
#'   # model-based estimation without regularization
#'   out_model_s <-
#'       model_pAMCE(formula = Y ~ gender + race,
#'            data = OnoBurden_data_small, reg = FALSE,
#'            pair_id = OnoBurden_data_small$pair_id,
#'            cluster_id = OnoBurden_data_small$id,
#'            target_dist  = target_dist_marginal_small,
#'            target_type = "marginal")
#'
#' \donttest{
#'   # Example
#'   data("OnoBurden")
#'   OnoBurden_data <- OnoBurden$OnoBurden_data
#'   OnoBurden_data_cong <- OnoBurden_data[OnoBurden_data$office == "Congress", ]
#'   target_dist_marginal <- OnoBurden$target_dist_marginal
#'
#'   # model-based estimation with regularization
#'   out_model <-
#'     model_pAMCE(formula = Y ~ gender + age + family + race + experience + party + pos_security,
#'                  data = OnoBurden_data_cong,
#'                  pair_id = OnoBurden_data_cong$pair_id,
#'                  cluster_id = OnoBurden_data_cong$id,
#'                  target_dist  = target_dist_marginal, target_type = "marginal")
#'  summary(out_model, factor_name = c("gender"))
#'
#'  # decompose the difference in the pAMCEs
#'  decompose_pAMCE(out_model, effect_name = c("gender", "Female"))
#' }
#' @description \code{model_pAMCE} implements the model-based approach to estimate the pAMCE. See de la Cuesta, Egami, and Imai (2019+) for details. More examples are available at the GitHub page of \code{factorEx}.
#' @references de la Cuesta, Egami, and Imai. (2019+). Improving the External Validity of Conjoint Analysis: The Essential Role of Profile Distribution. (Working Paper). Available at \url{https://scholar.princeton.edu/sites/default/files/negami/files/conjoint_profile.pdf}.
#' @references Egami and Imai. (2019). Causal Interaction in Factorial Experiments: Application to Conjoint Analysis. Journal of the American Statistical Association, Vol.114, No.526 (June), pp. 529–540. Available at \url{https://scholar.princeton.edu/sites/default/files/negami/files/causalint.pdf}.
#' @export


model_pAMCE <- function(formula,
                        formula_three = NULL,
                        data,
                        reg = TRUE,
                        ord_fac,
                        pair = FALSE, pair_id = NULL, cross_int = FALSE,
                        cluster_id = NULL,
                        target_dist,
                        target_type,
                        difference = FALSE,
                        cv_type = "cv.1Std", nfolds = 5,
                        boot = 100,
                        seed = 1234, numCores = NULL){

  ##################
  ## HouseKeeping ##
  ##################

  if(missing(pair_id) == TRUE) pair_id <- NULL
  if(missing(cluster_id) == TRUE) cluster_id <- seq(1:nrow(data))
  if(missing(target_type) == TRUE){
    target_type <- rep("marginal", length(target_dist))
  }

  if(all(target_type %in% c("marginal", "target_data")) == FALSE){
    warning(" 'target_type' should be 'marginal' or 'target_data' ")
  }

  if(length(target_type) == 1){
    if(target_type  == "marginal"){
      if(length(target_dist) != 1){
        target_dist <- list(target_dist)
      }
    }
  }

  if(length(target_type) != length(target_dist)){
    stop(" length of 'target_type' should be the same as length of 'target_dist' ")
  }

  if(class(target_dist) != "list"){
    target_dist <- list(target_dist)
  }
  if(is.list(target_dist)==FALSE){
    stop("target_dist should be 'list'.")
  }
  if(is.null(names(target_dist)) == TRUE){
    target_name <- paste("target_", seq(1:length(target_dist)), sep = "")
    names(target_dist) <- target_name
  }else{
    target_name <- names(target_dist)
  }

  if(is.null(numCores) == FALSE){
    if(numCores >= detectCores()) numCores <- detectCores() - 1
  }
  if(is.null(numCores)) numCores <- detectCores() - 1

  if(any(is.na(data))==TRUE){
    stop("Remove NA before using the function.")
  }
  if(pair==TRUE & is.null(pair_id)==TRUE){
    stop("When 'pair=TRUE', specify 'pair_id'.")
  }
  if(pair==TRUE & all(table(pair_id)==2)==FALSE){
    stop("When 'pair=TRUE', each of 'pair_id' should have two observations.")
  }
  if(is.null(pair_id) == FALSE) pair <- TRUE

  if(pair == FALSE){
    cross_int <- FALSE
  }

  if(difference == TRUE & length(target_dist) < 2){
    stop("if 'difference = TRUE', 'target_dist' should contain more than one distribution.")
  }
  if(boot < 500){
    message("Note: suggest 'boot' greater than 500 for final results.")
  }

  target_dist_orig <- target_dist
  # Add Marginal Distributions of Profiles used in Experiments
  Sample_Mar <- createDist(formula = formula, target_data = data,
                        exp_data = data, type  = "marginal")
  target_dist  <- c(list(Sample_Mar), target_dist)
  names(target_dist)[1] <- "sample"  # (marginal distributions used for randomization)
  target_name <- c("sample", target_name)
  target_type <- c("marginal", target_type)

  # ########################################
  # Check each target profile distribution #
  # ########################################
  for(i in 1:length(target_dist)){
    target_dist[[i]] <- checkDist(target_dist[[i]], type = target_type[i], formula = formula, data = data)
  }

  ## Create Marginals and 2d-Joints from "target_dist"
  marginal_dist <- list()
  for(i in 1:length(target_dist)){
    if(target_type[i]  == "marginal"){marginal_dist[[i]] <- target_dist[[i]]}
    # if(target_type[i]  == "joint"){marginal_dist[[i]] <- Joint2Marginal(target_dist[[i]])}
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
      # if(target_type[i]  == "joint"){joint_dist[[i]] <- target_dist[[i]]}
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

  # #################
  # Check Baselines #
  # #################
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

  # #############################
  # Renaming (for internal use)
  # #############################
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

  # Estimating pAMCE
  if(reg == FALSE){
    out <-  AME_estimate(formula = formula,
                         data = data,
                         pair = pair, pair_id = pair_id, cross_int = cross_int,
                         cluster = cluster_id,
                         marginal_dist = marginal_dist,
                         marginal_type = target_name,
                         joint_dist = joint_dist,
                         boot = boot, seed =  seed,
                         difference = difference, formula_three_c = formula_three_c)
  }else if(reg == TRUE){
    if(missing(ord_fac)) ord_fac <- rep(TRUE, (length(all.vars(formula)) - 1))
    out <- AME_estimate_collapse_genlasso(formula = formula,
                                          data = data,
                                          ord.fac = ord_fac,
                                          pair = pair, pair_id = pair_id, cross_int = cross_int,
                                          cluster = cluster_id,
                                          marginal_dist = marginal_dist,
                                          marginal_type = target_name,
                                          joint_dist = joint_dist,
                                          formula_three_c = formula_three_c,
                                          difference = difference,
                                          cv.type = cv_type,
                                          nfolds = nfolds,
                                          boot = boot,
                                          eps = 0.0001,
                                          numCores = numCores,
                                          seed = seed)
  }

  ## Approximate F-test
  if(is.null(formula_three_c) == FALSE){
    if(reg == TRUE)  coef_f <- apply(out$boot_coef, 2, mean)
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
    out$input$reg <- reg
    out$input$target_dist  <- target_dist_orig
    out$input$target_dist_name

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
    for(i in 1:length(out$AMCE)){
      match_level <- match(out$AMCE[[i]]$level, internal_level[[out$AMCE[[i]]$factor[1]]])
      out$AMCE[[i]]$factor <- rename_fac[,"original"][match(out$AMCE[[i]]$factor[1], rename_fac[,"internal"])]
      orignal_use <- original_level[[out$AMCE[[i]]$factor[1]]]
      out$AMCE[[i]]$level <- orignal_use[match_level]
    }
    names(out$AMCE) <- rename_fac[,"original"][match(names(out$AMCE), rename_fac[,"internal"])]

    # order of coefficients
    out$coef_order <- str_count(colnames(out$boot_coef), ":") + 1

    # coefficients
    for(i in 1:fac_size){
      colnames(out$boot_coef) <- gsub(rename_fac[(i+1), "internal"], rename_fac[(i+1), "original"], colnames(out$boot_coef))
      for(j in 1:length(internal_level[[i]])){
        colnames(out$boot_coef) <- gsub(internal_level[[i]][j], original_level[[i]][j], colnames(out$boot_coef))
      }
    }
    names(out$coef)  <- colnames(out$boot_coef)
  }
  out$approach <- "model_based"

  class(out)  <- c(class(out), "pAMCE")
  return(out)
}
