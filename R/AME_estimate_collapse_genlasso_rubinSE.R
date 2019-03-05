#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param pair whether we use the paired-conjoint design
#' @param pair_id id for paired-conjoint design. required when 'pair = TRUE'
#' @export

AME_estimate_collapse_genlasso_rubinSE <- function(formula,
                                                   data,
                                                   ord.fac,
                                                   pair=FALSE, pair_id = NULL,
                                                   cluster = NULL,
                                                   marginal_dist,
                                                   marginal_type,
                                                   difference = FALSE,
                                                   cv.type = "cv.1Std",
                                                   seed = 1234, nfolds = 2,
                                                   boot = 100,
                                                   eps = 0.0001){

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
  # Create two-way interaction formula ----------
  factor_l <- length(all.vars(formula)[-1])
  combMat <- combn(factor_l,2); intNames <- c()
  for(k in 1:ncol(combMat)){
    intNames[k] <- paste(all.vars(formula)[-1][combMat[1,k]], "*", all.vars(formula)[-1][combMat[2,k]], sep = "")
  }
  formula_full <- as.formula(paste(all.vars(formula)[1], "~", paste(intNames, collapse = "+"), sep=""))

  # Baseline ----------
  baseline <- lapply(model.frame(formula, data=data)[,-1],
                     FUN = function(x) levels(x)[1])

  fac.level <- unlist(lapply(model.frame(formula, data=data)[,-1],
                             FUN = function(x) length(levels(x))))

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

  # base
  tableAME_base <- AME.fit(formula,
                           data = data, pair = pair,
                           marginal_dist = marginal_dist,
                           marginal_dist_u_list = marginal_dist_u_list,
                           marginal_dist_u_base = marginal_dist_u_base,
                           marginal_type = marginal_type,
                           difference = difference)
  tableAME_base$estimate <- NULL




  # Collapsing
  if(pair == TRUE)  data$pair_id <- pair_id
  if(pair == FALSE) data$pair_id <- NA
  data$cluster <- cluster

  if(missing(ord.fac)) ord.fac <- rep(TRUE, factor_l)

  table_AME <- AME.collapse.genlasso.crossfit.boot.rubinSE(formula = formula,
                                                           data = data,
                                                           pair = pair,
                                                           fac.level = fac.level, ord.fac = ord.fac,
                                                           marginal_dist = marginal_dist,
                                                           marginal_type = marginal_type,
                                                           difference = difference,
                                                           boot = boot,
                                                           cv.type = cv.type,
                                                           nfolds = nfolds, seed = seed,
                                                           tableAME_base = tableAME_base,
                                                           eps = eps)


  ## For Each Factor
  AME <- list()
  for(g in 1:length(unique(table_AME$factor))){
    AME[[g]] <- table_AME[table_AME$factor == unique(table_AME$factor)[g], ]
  }
  names(AME) <- unique(table_AME$factor)
  type_all   <- unique(table_AME$type)
  type_difference   <- setdiff(unique(table_AME$type), marginal_type)

  input  <- list("formula" = formula, "data" = data,
                 "pair" = pair, "pair_id" = pair_id,
                 "marginal_dist" = marginal_dist,
                 "marginal_type" = marginal_type, "difference" = difference)

  output <- list("AME" = AME, "baseline" = baseline,
                 "type_all" = type_all, "type_difference" = type_difference,
                 "input" = input)
  return(output)
}





