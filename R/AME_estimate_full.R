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
#' @import sandwich
#' @importFrom igraph graph_from_adjacency_matrix
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
                              seed = 1234){

  cat("Using Version 0.1.2\n")

  if((type %in% c("No-Reg","gash-anova", "genlasso")) == FALSE){
    warning(" 'type' should be one of 'No-Reg', 'gash-anova' and 'genalsso' ")
  }

  if(missing(pair_id) == TRUE) pair_id <- NULL
  if(missing(cluster) == TRUE) cluster <- seq(1:nrow(data))

  set.seed(seed)

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
                                      boot = boot)

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
                                          eps = 0.0001)
  }
  return(out)
}
