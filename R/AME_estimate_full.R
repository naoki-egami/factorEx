#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param type "No-Reg", "gash-anova", or "genlasso"
#' @param ord.fac whether we assume each factor is ordered. When not specified, we assume all of them are ordered.
#' @param pair whether we use the paired-conjoint design
#' @param pair_id id for paired-conjoint design. required when 'pair = TRUE'
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
                              pair = FALSE, pair_id = NULL,
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

  if((type %in% c("No-Reg","gash-anova", "genlasso")) == FALSE){
    warning(" 'type' should be one of 'No-Reg', 'gash-anova' and 'genalsso' ")
  }

  if(missing(pair_id) == TRUE) pair_id <- NULL
  if(missing(cluster) == TRUE) cluster <- seq(1:nrow(data))

  if(is.null(numCores) == FALSE){
    if(numCores >= detectCores()) numCores <- detectCores() - 1
  }
  if(is.null(numCores)) numCores <- detectCores() - 1

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
  baseline <- lapply(model.frame(formula, data=data)[,-1], FUN = function(x) levels(x)[1])
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

  if(class(marginal_dist) != "list"){
      marginal_dist <- list(marginal_dist)
  }


  if(type == "No-Reg"){
    out <-  AME_estimate(formula = formula,
                         data = data,
                         pair = pair, pair_id = pair_id,
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
                                          pair = pair, pair_id = pair_id,
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
  return(out)
}
