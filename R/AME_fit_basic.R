#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param pair whether we use the paired-conjoint design
#' @param pair_id id for paired-conjoint design. required when 'pair = TRUE'
#' @export


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

AME.fit.STD.se <- function(formula,
                           data,
                           pair = FALSE,
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
    cluster <- data$cluster
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
  vcovInt <- cluster_se_glm(main_lm, cluster = as.factor(cluster))
  colnames(vcovInt) <- rownames(vcovInt) <- base_name

  # Estimate AMEs ----------
  table_AME <- c()
  for(m in 1:nrow(marginal_dist_u_base)){
    coef_focus <- coefInt[grep(marginal_dist_u_base$level[m], names(coefInt), fixed = T)]
    vcov_focus <- vcovInt[grep(marginal_dist_u_base$level[m], rownames(vcovInt), fixed = T),
                          grep(marginal_dist_u_base$level[m], colnames(vcovInt), fixed = T)]
    if(length(coef_focus) > 0){
      estNames <- gsub(paste(marginal_dist_u_base$level[m], ":", sep = ""), "", names(coef_focus), fixed = T)
      estNames <- gsub(paste(":", names(coef_focus)[1], sep = ""), "", estNames, fixed = T)

      table_AME_m <- data.frame(matrix(NA, ncol = 0, nrow=1))
      table_AME_m$type <- "STD"
      table_AME_m$factor   <- marginal_dist[[1]][m,1]; table_AME_m$level <- marginal_dist[[1]][m,2]
      table_AME_m$estimate <- coef_focus; table_AME_m$se <- sqrt(vcov_focus)
      table_AME <- rbind(table_AME, table_AME_m)
    }
  }
  colnames(table_AME) <- c("type", "factor", "level", "estimate", "se")
  return(table_AME)
}

cluster_se_glm <- function(model, cluster){

  #  Drop unused cluster indicators, if cluster var is a factor
  if (class(cluster) == "factor") {
    cluster <- droplevels(cluster)
  }

  if (nrow(model.matrix(model)) != length(cluster)) {
    stop("check your data: cluster variable has different N than model - you may have observations with missing data")
  }

  M <- length(unique(cluster))
  N <- length(cluster)
  K <- model$rank

  ## if(M<50) {
  ##     warning("Fewer than 50 clusters, variances may be unreliable (could try block bootstrap instead).")
  ## }

  dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
  uj  <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
  rcse.cov <- dfc * sandwich(model, meat. = crossprod(uj)/N)
  return(rcse.cov)
}

cluster_se_glm <- function(model, cluster){

  #  Drop unused cluster indicators, if cluster var is a factor
  if (class(cluster) == "factor") {
    cluster <- droplevels(cluster)
  }

  if (nrow(model.matrix(model)) != length(cluster)) {
    stop("check your data: cluster variable has different N than model - you may have observations with missing data")
  }

  M <- length(unique(cluster))
  N <- length(cluster)
  K <- model$rank

  ## if(M<50) {
  ##     warning("Fewer than 50 clusters, variances may be unreliable (could try block bootstrap instead).")
  ## }

  dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
  uj  <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
  rcse.cov <- dfc * sandwich(model, meat. = crossprod(uj)/N)
  return(rcse.cov)
}

