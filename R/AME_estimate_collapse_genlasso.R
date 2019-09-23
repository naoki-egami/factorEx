## (Internal Function) Estimating pAMCE with generalized lasso regularization
## @param formula formula
## @param data data
## @param ord.fac Whether we assume each factor is ordered. When not specified, we assume all of them are ordered
## @param pair Whether we use the paired-conjoint design
## @param pair_id Unique id for paired-conjoint design. Required when 'pair = TRUE'
## @param cross_int Include interactions across profiles. Default is FALSE
## @param cluster Unique identifiers for computing cluster standard errors
## @param marginal_dist Target profile marginal distributions to be used. This argument should be `list`
## @param marginal_type Names of target profile marginal distributions.
## @param joint_dist Target profile 2-dimensional joint distributions to be used. This argument should be `list`
## @param formula_three_c Formula for three-way interactions (optional)
## @param difference Whether we compute the differences between the multiple pAMCEs. Default is FALSE.
## @param cv.type (optimal only when `reg = TRUE``)  `cv.1Std`` (stronger regularization; default) or `cv.min` (weaker regularization).
## @param nfolds Number of cross validation folds. Default is 5.
## @param boot The number of bootstrap samples.
## @param seed Seed for bootstrap.
## @param numCores Number of cores to be used for parallel computing. If not specified, detect the number of available cores internally.
## @param eps (internal use) tolerance range for collapsing

AME_estimate_collapse_genlasso <- function(formula,
                                           data,
                                           ord.fac,
                                           pair = FALSE, pair_id = NULL, cross_int = TRUE,
                                           cluster = NULL,
                                           marginal_dist,
                                           marginal_type,
                                           joint_dist,
                                           formula_three_c,
                                           difference = FALSE,
                                           cv.type = "cv.1Std",
                                           nfolds = 5,
                                           boot = 100,
                                           eps = 0.0001,
                                           numCores,
                                           seed){

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
  if(difference==TRUE & length(marginal_type) < 2){
    stop("if 'difference = TRUE', marginal_dist should contain more than one distribution.")
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
  all.fac <- all(unlist(lapply(model.frame(formula, data=data)[,-1],
                               FUN=function(x) is.element("factor",class(x)))))
  if(all.fac == FALSE) stop("Design matrix should contain only factors."); rm(all.fac)

  ###########
  ## Setup ##
  ###########
  # foctorl number  ----------
  factor_l <- length(all.vars(formula)[-1])

  # Baseline ----------
  baseline <- lapply(model.frame(formula, data=data)[,-1],
                     FUN = function(x) levels(x)[1])

  fac.level <- unlist(lapply(model.frame(formula, data=data)[,-1],
                             FUN = function(x) length(levels(x))))

  original_level <- lapply(model.frame(formula, data=data)[,-1], FUN = function(x) levels(x))

  # Check ----------
  marginal_dist_u0 <- marginal_dist[[1]]
  marginal_dist_u <- data.frame(matrix(NA, ncol=0, nrow=nrow(marginal_dist_u0)))
  marginal_dist_u$level <- paste(marginal_dist_u0[,1], marginal_dist_u0[,2],sep="")
  marginal_dist_u$prop  <- marginal_dist_u0[,3]
  mar_name <- colnames(model.matrix(formula, data=data))[-1]
  if(all(is.element(mar_name, marginal_dist_u[,1])) == FALSE){
    stop("The first and second columns of `marginal.dist' should have correct names for 'factors' and 'levels'")
  }
  rm(marginal_dist_u)

  # Transform marginal_dist (for internal simplisity) ----------
  marginal_dist_u_list <- list()
  for(z in 1:length(marginal_dist)){
    marginal_dist_u_list[[z]] <- data.frame(matrix(NA, ncol=0, nrow=nrow(marginal_dist[[z]])))
    marginal_dist_u_list[[z]]$level <- paste(marginal_dist[[z]][,1], marginal_dist[[z]][,2],sep="")
    marginal_dist_u_list[[z]]$prop  <- marginal_dist[[z]][,3]
  }
  marginal_dist_u_base <- marginal_dist_u_list[[1]]

  if(is.null(joint_dist) == TRUE){
    joint_dist_u_list <- NULL
  }else{
    joint_dist_u_list <- list()
    for(z in 1:length(joint_dist)){
      temp <- do.call("rbind", joint_dist[[z]])
      joint_dist_u_list[[z]] <- data.frame(matrix(NA, ncol=0, nrow=nrow(temp)))
      joint_dist_u_list[[z]]$level_1 <- paste(temp[,1], temp[,3],sep="")
      joint_dist_u_list[[z]]$level_2 <- paste(temp[,2], temp[,4],sep="")
      joint_dist_u_list[[z]]$prop  <- temp[,5]
    }
  }

  # three-way
  three_way <- is.null(formula_three_c) == FALSE

  # base (just keep track of names)
  fitAME_base <- AME.fit(formula_full = formula,
                         data = data, pair = pair, cross_int = cross_int,
                         marginal_dist = marginal_dist,
                         marginal_dist_u_list = marginal_dist_u_list,
                         marginal_dist_u_base = marginal_dist_u_base,
                         marginal_type = marginal_type,
                         difference = difference,
                         three_way = three_way, est_AME = TRUE)

  tableAME_base <- fitAME_base$table_AME
  tableAME_base$estimate <- NULL
  coefAME_base  <- coefMake(original_level, cross_int = cross_int, formula_three_c = formula_three_c)

  # Collapsing
  if(pair == TRUE)  data$pair_id <- pair_id
  if(pair == FALSE) data$pair_id <- NA
  data$cluster <- cluster

  if(missing(ord.fac)) ord.fac <- rep(TRUE, factor_l)

  table_AME_f <- AME.collapse.genlasso.crossfit.boot(formula = formula,
                                                     data = data,
                                                     pair = pair, cross_int = cross_int,
                                                     fac.level = fac.level, ord.fac = ord.fac,
                                                     marginal_dist = marginal_dist,
                                                     marginal_type = marginal_type,
                                                     joint_dist_u_list =  joint_dist_u_list,
                                                     formula_three_c = formula_three_c,
                                                     difference = difference,
                                                     boot = boot,
                                                     cv.type = cv.type,
                                                     nfolds = nfolds,
                                                     tableAME_base = tableAME_base,
                                                     coefAME_base_l = coefAME_base,
                                                     eps = eps,
                                                     numCores = numCores,
                                                     seed = seed)

  table_AME <- table_AME_f$fit
  boot_AME  <- table_AME_f$fit.mat
  boot_coef <- table_AME_f$coef.mat
  colnames(boot_coef) <- coefAME_base

  ## For Each Factor
  AME <- list()
  for(g in 1:length(unique(table_AME$factor))){
    AME[[g]] <- table_AME[table_AME$factor == unique(table_AME$factor)[g], ]
  }
  names(AME) <- unique(table_AME$factor)
  type_all   <- unique(table_AME$type)
  type_difference   <- setdiff(unique(table_AME$type), marginal_type)

  coef_final <- apply(boot_coef, 2, mean)

  input  <- list("formula" = formula, "data" = data,
                 "pair" = pair, "pair_id" = pair_id, "cross_int" = cross_int,
                 "marginal_dist" = marginal_dist,
                 "marginal_dist_u_list" = marginal_dist_u_list,
                 "marginal_dist_u_base" = marginal_dist_u_base,
                 "marginal_type" = marginal_type, "difference" = difference)

  output <- list("AMCE" = AME, "baseline" = baseline, "coef" =  coef_final,
                 "type_all" = type_all, "type_difference" = type_difference,
                 "boot_AMCE" = boot_AME, "boot_coef" = boot_coef,
                 "input" = input)
  return(output)
}
