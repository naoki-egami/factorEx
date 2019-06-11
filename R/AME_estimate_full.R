#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param type "No-Reg", "gash-anova", or "genlasso"
#' @param ord.fac whether we assume each factor is ordered. When not specified, we assume all of them are ordered.
#' @param pair whether we use the paired-conjoint design
#' @param cross_int include interactions across profiles
#' @param cluster id for cluster
#' @param marginal_dist marginal distributions
#' @param marginal_type names of marginal distributions
#' @param difference whether we compute the difference between different estimators
#' @param family (only when 'type = gash-anova') when outcomes are binary, "binomial". when outcomes are continuous, "gaussian"
#' @param nway (only when `type = gash-anova`) Should be 1 almost always.
#' @param cv.collapse.cost (only when `type = gash-anova`) a grid for cross-validation in gash-anova
#' @param cv.type (when type = gash-anova or genlasso) `cv.1Std`` (stronger) or `cv.min` (weaker).
#' @param boot (when type = gash-anova or genlasso) the number of bootstrap
#' @param seed seed for bootstrap
#' @importFrom FindIt cv.CausalANOVA CausalANOVA
#' @importFrom prodlim row.match
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom pbmcapply pbmclapply
#' @importFrom pbapply pblapply
#' @import parallel
#' @import arm
#' @export

AME_estimate_full <- function(formula,
                              data,
                              type = "genlasso",
                              ord.fac,
                              pair = FALSE, pair_id = NULL, cross_int = TRUE,
                              cluster = NULL,
                              marginal_dist,
                              marginal_type,
                              difference = FALSE,
                              family = "binomial",
                              nway = 1,
                              cv.collapse.cost = c(0.1, 0.3, 0.5, 0.7),
                              cv.type = "cv.1Std", nfolds = 5,
                              boot = 100,
                              seed = 1234,
                              numCores = NULL){

  cat("Using version-conditional_effect:\n")

  ###########
  ## Check ##
  ###########
  if((type %in% c("No-Reg","gash-anova", "genlasso")) == FALSE){
    warning(" 'type' should be one of 'No-Reg', 'gash-anova' and 'genalsso' ")
  }

  if(missing(pair_id) == TRUE) pair_id <- NULL
  if(missing(cluster) == TRUE) cluster <- seq(1:nrow(data))

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
  if(pair == FALSE){
    cross_int <- FALSE
  }
  # if(is.null(pair_var) == FALSE){
  #   if(all(is.element(pair_var, all.vars(formula))) == FALSE){
  #     stop(" 'pair_var' should be variables listed in 'formula' ")}
  #   if(pair == FALSE){
  #     stop(" 'pair_var' is ignored when 'pair=FALSE' ")}
  # }
  if(difference==TRUE & length(marginal_type) < 2){
    stop("if 'difference = TRUE', marginal_dist should contain more than one distribution.")
  }
  if(class(marginal_dist) != "list"){
    marginal_dist <- list(marginal_dist)
  }
  if(is.list(marginal_dist)==FALSE){
    stop("marginal_dist should be 'list'.")
  }
  marginal_name_check <- lapply(marginal_dist, colnames)
  marginal_name_check_all <- all(unlist(lapply(marginal_name_check,
                                               function(x) all(x == c("factor", "levels", "prop")))))
  if(marginal_name_check_all == FALSE){
    stop(" 'colnames' of 'marginal_dist' should be c('factor', 'levels', 'prop') ")
  }
  if(length(marginal_type) > 1){
    for(z in 2:length(marginal_type)){
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

  if(type == "No-Reg"){
    out <-  AME_estimate(formula = formula,
                         data = data,
                         pair = pair, pair_id = pair_id, cross_int = cross_int,
                         cluster = cluster,
                         marginal_dist = marginal_dist,
                         marginal_type = marginal_type,
                         difference = difference)
  }else if(type == "gash-anova"){
    if(missing(ord.fac)) ord.fac <- rep(TRUE, (length(all.vars(formula)) - 1))
    out <- AME_estimate_collapse_gash(formula = formula,
                                      data = data,
                                      ord.fac = ord.fac,
                                      pair = pair, pair_id = pair_id,
                                      cluster = cluster,
                                      marginal_dist = marginal_dist,
                                      marginal_type = marginal_type,
                                      difference = difference,
                                      family = family,
                                      nway = nway,
                                      cv.collapse.cost = cv.collapse.cost,
                                      cv.type = cv.type,
                                      boot = boot,
                                      numCores = numCores,
                                      seed = seed)

  }else if(type == "genlasso"){
    if(missing(ord.fac)) ord.fac <- rep(TRUE, (length(all.vars(formula)) - 1))
    out <- AME_estimate_collapse_genlasso(formula = formula,
                                          data = data,
                                          ord.fac = ord.fac,
                                          pair = pair, pair_id = pair_id, cross_int = cross_int,
                                          cluster = cluster,
                                          marginal_dist = marginal_dist,
                                          marginal_type = marginal_type,
                                          difference = difference,
                                          cv.type = cv.type,
                                          nfolds = nfolds,
                                          boot = boot,
                                          eps = 0.0001,
                                          numCores = numCores,
                                          seed = seed)
  }

  # ##############
  # Name back
  # ##############
  {
    out$input$formula <- formula_orig
    out$input$data <- data_orig
    out$input$marginal_dist <- marginal_dist_orig
    out$baseline <- baseline_orig

    # AME
    for(i in 1:length(out$AME)){
      match_level <- match(out$AME[[i]]$level, internal_level[[out$AME[[i]]$factor[1]]])
      out$AME[[i]]$factor <- rename_fac[,"original"][match(out$AME[[i]]$factor[1], rename_fac[,"internal"])]
      orignal_use <- original_level[[out$AME[[i]]$factor[1]]]
      out$AME[[i]]$level <- orignal_use[match_level]
    }
    names(out$AME) <- rename_fac[,"original"][match(names(out$AME), rename_fac[,"internal"])]
    # coefficients
    if(type != "No-Reg"){
      for(i in 1:fac_size){
        colnames(out$boot_coef) <- gsub(rename_fac[(i+1), "internal"], rename_fac[(i+1), "original"], colnames(out$boot_coef))
        for(j in 1:length(internal_level[[i]])){
          colnames(out$boot_coef) <- gsub(internal_level[[i]][j], original_level[[i]][j], colnames(out$boot_coef))
        }
      }
    }
  }

  return(out)
}
