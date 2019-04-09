#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param pair whether we use the paired-conjoint design
#' @param pair_id id for paired-conjoint design. required when 'pair = TRUE'
#' @importFrom sandwich sandwich estfun
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
    X1 <- model.matrix(formula_full, data=data1)
    X2 <- model.matrix(formula_full, data=data2)
    X <- cbind(1, X1[,-1] - X2[,-1])
    y <- model.frame(formula_full, data=data1)[ ,1]
    # base_name <- c("(Intercept)", colnames(X1))
    ind_b <- attr(X1, "assign")
  }else{
    cluster_original <- data$cluster
    X <- model.matrix(formula_full, data=data)
    y <- model.frame(formula_full, data=data)[,1]
    # base_name <- colnames(X)
    side <- NULL
    ind_b <- attr(X, "assign")
  }

  # Fit the model ----------
  main_lm <- lm(y ~ X - 1)
  coefInt <- coef(main_lm)
  # coefInt <- coefInt[is.na(coefInt) == FALSE]
  coefInt[is.na(coefInt) == TRUE] <- 0 ## Insert 0 to non-existing pairs
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
  out <- list("table_AME" = table_AME, "coef" = coefInt, "ind_b" = ind_b)
  return(out)
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
  coefInt <- summary(main_lm)$coef[,1]
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
  coefInt <- summary(main_lm)$coef[,1]
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


#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param pair whether we use the paired-conjoint design
#' @param pair_id id for paired-conjoint design. required when 'pair = TRUE'
#' @export

AME_from_boot <- function(outAME, difference = FALSE){

  marginal_dist <- outAME$input$marginal_dist
  marginal_dist_u_list  <- outAME$input$marginal_dist_u_list
  marginal_dist_u_base  <- outAME$input$marginal_dist_u_base
  marginal_type  <- outAME$input$marginal_type

  boot_coef <- outAME$boot_coef

  # Estimate AMEs ----------
  table_AME <- c()
  for(m in 1:nrow(marginal_dist_u_base)){
    coef_focus <- boot_coef[, grep(marginal_dist_u_base$level[m], colnames(boot_coef), fixed = T)]
    # vcov_focus <- vcovInt[grep(marginal_dist_u_base$level[m], rownames(vcovInt), fixed = T),
    #                       grep(marginal_dist_u_base$level[m], colnames(vcovInt), fixed = T)]
    if(length(coef_focus) > 0){
      estNames <- gsub(paste(marginal_dist_u_base$level[m], ":", sep = ""), "", colnames(coef_focus), fixed = T)
      estNames <- gsub(paste(":", colnames(coef_focus)[1], sep = ""), "", estNames, fixed = T)
      table_AME_m <- c()
      # For each marginal distribution,
      for(z in 1:length(marginal_dist_u_list)){
        marginal_dist_u <- marginal_dist_u_list[[z]]
        # Find weights
        coef_prop <- c(1, as.numeric(as.character(marginal_dist_u[match(estNames, marginal_dist_u[, "level"]), "prop"]))[-1])
        # Compute AMEs
        coef_AME <- mean(coef_focus %*% coef_prop)
        # coef_AME <- sum(coef_focus * coef_prop)
        se_AME   <- sd(coef_focus %*% coef_prop)
        # se_AME <- sqrt(coef_prop%*%vcov_focus%*%coef_prop)
        AME <- data.frame(matrix(NA, ncol = 0, nrow=1))
        AME$type <- marginal_type[z]
        AME$factor   <- marginal_dist[[z]][m,1]; AME$level <- marginal_dist[[z]][m,2]
        AME$estimate <- coef_AME;
        AME$se <- se_AME
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
          coef_AME_dif <- mean(coef_focus %*% coef_prop_d)
          # coef_AME_dif <- sum(coef_focus * coef_prop_d)
          se_AME_dif <- sd(coef_focus %*% coef_prop_d)
          # se_AME_dif <- sqrt(coef_prop_d%*%vcov_focus%*%coef_prop_d)
          AME_dif <- data.frame(matrix(NA, ncol = 0, nrow=1))
          AME_dif$type <- paste(marginal_type[z],"-",marginal_type[1],sep="")
          AME_dif$factor   <- marginal_dist[[z]][m,1]; AME_dif$level <- marginal_dist[[z]][m,2]
          AME_dif$estimate <- coef_AME_dif;
          AME_dif$se <- se_AME_dif
          table_AME_m <- rbind(table_AME_m, AME_dif)
        }
      }
      table_AME <- rbind(table_AME, table_AME_m)
    }
  }
  colnames(table_AME) <- c("type", "factor", "level", "estimate", "se")
  out <- list("table_AME" = table_AME)
  return(out)
}

#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param pair whether we use the paired-conjoint design
#' @param pair_id id for paired-conjoint design. required when 'pair = TRUE'
#' @export

cAME_from_boot <- function(outAME, factor_name, level_name, difference = TRUE){

  marginal_dist <- outAME$input$marginal_dist
  marginal_dist_u_list  <- outAME$input$marginal_dist_u_list
  marginal_dist_u_base  <- outAME$input$marginal_dist_u_base
  marginal_type  <- outAME$input$marginal_type
  formula <- outAME$input$formula
  data    <- outAME$input$data

  leveluse <- lapply(model.frame(formula, data=data)[,-1], FUN = function(x) levels(x)[-1])
  leveluse_all <- lapply(model.frame(formula, data=data)[,-1], FUN = function(x) levels(x))

  boot_coef <- outAME$boot_coef

  # Estimate conditional AMEs ----------
  effect_name <- paste(factor_name, level_name, sep="")
  coef_focus  <- boot_coef[, grep(effect_name, colnames(boot_coef), fixed = T)]
  # vcov_focus <- vcovInt[grep(effect_name, rownames(vcovInt), fixed = T),
  #                       grep(effect_name, colnames(vcovInt), fixed = T)]

  if(length(coef_focus) == 0){
    stop(" 'level_name' cannot take the baseline level or undefined levels of the specified factor.")
  }

  factor_all <- all.vars(formula)[-1]
  part_pos <- lapply(leveluse[!(names(leveluse) %in% factor_name)], length)
  end_pos   <- cumsum(part_pos) + 1
  start_pos <- c(2, end_pos[-length(end_pos)] + 1)

  estNames <- gsub(paste(effect_name, ":", sep = ""), "", colnames(coef_focus), fixed = T)
  estNames <- gsub(paste(":", effect_name, sep = ""), "", estNames, fixed = T)

  # Estimate Conditional AMEs
  table_cAME <- table_cProp <- c()
  for(m in 1:(length(factor_all) - 1)){
    ind_focus_mat <- matrix(1, ncol = (as.numeric(part_pos[m][1]) + 1), nrow = ncol(coef_focus))
    leveluse_focus <- leveluse_all[!(names(leveluse_all) %in% factor_name)][m]
    leveluse_name  <- paste(names(leveluse_focus), unlist(leveluse_focus), sep="")
    table_cAME_m <- table_cProp_m <- c()
    # For each marginal distribution,
    for(z in 1:length(marginal_dist)){
      marginal_dist_u <- marginal_dist_u_list[[z]]
      # Find weights
      coef_prop <- c(1, as.numeric(as.character(marginal_dist_u[match(estNames, marginal_dist_u[, "level"]), "prop"]))[-1])
      coef_prop_use <- as.numeric(as.character(marginal_dist_u[match(leveluse_name, marginal_dist_u[, "level"]), "prop"]))
      prop_focus_mat_b <- ind_focus_mat * coef_prop
      prop_focus_mat_b[start_pos[m]:end_pos[m], 1:ncol(prop_focus_mat_b)] <- cbind(0, diag(1, (end_pos[m] - start_pos[m] + 1)))
      # prop_focus_mat <- t(prop_focus_mat_b) * coef_prop_use
      prop_focus_mat <- t(prop_focus_mat_b)

      # Compute AMEs
      coef_cAME_b <- prop_focus_mat %*% t(coef_focus)
      coef_cAME   <- apply(coef_cAME_b, 1, mean)
      se_cAME     <- apply(coef_cAME_b, 1, sd)
      cAME <- data.frame(matrix(NA, ncol = 0, nrow = length(coef_cAME)))
      cAME$type <- marginal_type[z]
      cAME$factor   <- names(leveluse_focus); cAME$level <- unlist(leveluse_focus)
      cAME$estimate <- coef_cAME; cAME$se <- se_cAME
      table_cAME_m <- rbind(table_cAME_m, cAME)

      # Store Distribution
      cProp <- data.frame(matrix(NA, ncol = 0, nrow = length(coef_cAME)))
      cProp$type <- marginal_type[z]
      cProp$factor   <- names(leveluse_focus); cProp$level <- unlist(leveluse_focus)
      cProp$prop <- coef_prop_use
      table_cProp_m <- rbind(table_cProp_m, cProp)
    }
    if(difference == TRUE){
      for(z in 2:length(marginal_dist)){
        marginal_dist_u <- marginal_dist_u_list[[z]]
        # Find weights
        coef_prop <- c(1, as.numeric(as.character(marginal_dist_u[match(estNames, marginal_dist_u[, "level"]), "prop"]))[-1])
        coef_prop_use <- as.numeric(as.character(marginal_dist_u[match(leveluse_name, marginal_dist_u[, "level"]), "prop"]))
        coef_prop0 <- c(1, as.numeric(as.character(marginal_dist_u_base[match(estNames, marginal_dist_u_base[, "level"]), "prop"]))[-1])
        coef_prop_use0 <- as.numeric(as.character(marginal_dist_u_base[match(leveluse_name, marginal_dist_u_base[, "level"]), "prop"]))

        prop_focus_mat_b <- ind_focus_mat * coef_prop
        prop_focus_mat_b[start_pos[m]:end_pos[m], 1:ncol(prop_focus_mat_b)] <- cbind(0, diag(1, (end_pos[m] - start_pos[m] + 1)))
        # prop_focus_mat <- t(prop_focus_mat_b) * coef_prop_use
        prop_focus_mat <- t(prop_focus_mat_b)

        prop_focus_mat_b0 <- ind_focus_mat * coef_prop0
        prop_focus_mat_b0[start_pos[m]:end_pos[m], 1:ncol(prop_focus_mat_b0)] <- cbind(0, diag(1, (end_pos[m] - start_pos[m] + 1)))
        # prop_focus_mat0 <- t(prop_focus_mat_b0) * coef_prop_use0
        prop_focus_mat0 <- t(prop_focus_mat_b0)

        # Compute difference in AMEs
        prop_focus_mat_d <- (prop_focus_mat - prop_focus_mat0)
        coef_cAME_dif_b  <- prop_focus_mat_d %*% t(coef_focus)
        coef_cAME_dif    <- apply(coef_cAME_dif_b, 1, mean)
        se_cAME_dif <- apply(coef_cAME_dif_b, 1, sd)
        cAME_dif <- data.frame(matrix(NA, ncol = 0, nrow=length(coef_cAME_dif)))
        cAME_dif$type <- paste(marginal_type[z],"-",marginal_type[1],sep="")
        cAME_dif$factor   <- names(leveluse_focus); cAME_dif$level <- unlist(leveluse_focus)
        cAME_dif$estimate <- coef_cAME_dif; cAME_dif$se <- se_cAME_dif
        table_cAME_m <- rbind(table_cAME_m, cAME_dif)

        # Store Distribution
        cProp_dif <- data.frame(matrix(NA, ncol = 0, nrow = length(coef_cAME_dif)))
        cProp_dif$type <- paste(marginal_type[z],"-",marginal_type[1],sep="")
        cProp_dif$factor   <- names(leveluse_focus); cProp_dif$level <- unlist(leveluse_focus)
        cProp_dif$prop <- coef_prop_use - coef_prop_use0
        table_cProp_m <- rbind(table_cProp_m, cProp_dif)
      }
    }
    table_cAME_m <- table_cAME_m[order(table_cAME_m$level),]
    table_cAME <- rbind(table_cAME, table_cAME_m)

    table_cProp_m <- table_cProp_m[order(table_cProp_m$level),]
    table_cProp <- rbind(table_cProp, table_cProp_m)
  }
  colnames(table_cAME)  <- c("type", "factor", "level", "estimate", "se")
  colnames(table_cProp) <- c("type", "factor", "level", "prop")

  ## For Each Factor
  cAME <- cProp <- list()
  for(g in 1:length(unique(table_cAME$factor))){
    cAME[[g]]  <- table_cAME[table_cAME$factor == unique(table_cAME$factor)[g], ]
    cProp[[g]] <- table_cProp[table_cProp$factor == unique(table_cProp$factor)[g], ]
  }
  names(cAME) <- names(cProp) <- unique(table_cAME$factor)
  type_all   <- unique(table_cAME$type)
  type_difference   <- setdiff(unique(table_cAME$type), marginal_type)

  ## Reorder Levels
  marginal_use <- marginal_dist[[1]]
  level_order <- tapply(marginal_use$levels, marginal_use$factor, unique)
  level_order <- level_order[match(names(cAME), names(level_order))]

  for(i in 1:length(cAME)){
    cAME[[i]]$level <- factor(cAME[[i]]$level, levels = level_order[[i]])
    cAME[[i]]  <- cAME[[i]][order(cAME[[i]]$level), ]
    cAME[[i]]$level <- as.character(cAME[[i]]$level)
    cProp[[i]]$level <- factor(cProp[[i]]$level, levels = level_order[[i]])
    cProp[[i]] <- cProp[[i]][order(cProp[[i]]$level), ]
    cProp[[i]]$level <- as.character(cProp[[i]]$level)
  }

  output <- list("cAME" = cAME, "cProp" = cProp,
                 "difference" = difference,
                 "type_all" = type_all,
                 "type_difference" = type_difference,
                 "warnings" = warnings)
  return(output)
}


#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param pair whether we use the paired-conjoint design
#' @param pair_id id for paired-conjoint design. required when 'pair = TRUE'
#' @export

plot_cAME <- function(came.out,
                      factor_name,
                      main = c("Conditional AMEs"),
                      col,
                      pch, mar,
                      legend_pos = "topright",
                      plot_all = TRUE,
                      plot_difference = "add",
                      marginal_prop, marginal_effect,
                      cex = 1){

  type_all <- came.out$type_all
  type_difference <- came.out$type_difference

  if(all(is.element(factor_name, names(came.out$cAME))) == FALSE){
    stop(" 'factor_name' can only take a subset of factors estimated in 'cAME.estimate'")
  }
  if(missing(col)){
    if(plot_difference == "add")  col <- rep("black", length(type_all))
    if(plot_difference == "no")   col <- rep("black", length(type_all) - length(type_difference))
    if(plot_difference == "only") col <- rep("black", length(type_difference))
  }
  if(missing(pch)){
    if(plot_difference == "add")  pch <- rep(19, length(type_all))
    if(plot_difference == "no")   pch <- rep(19, length(type_all) - length(type_difference))
    if(plot_difference == "only") pch <- rep(19, length(type_difference))
  }
  if(missing(mar)){
    mar <- 6
  }
  if(came.out$difference == FALSE & plot_difference == "only"){
    stop(" 'cAME.estiamte' should take `difference = TRUE` to plot their estimates")
  }

  ## type_plot
  if(plot_difference == "add")  type_plot <- type_all
  if(plot_difference == "no")   type_plot <- setdiff(type_all, type_difference)
  if(plot_difference == "only") type_plot <- type_difference

  if(missing(marginal_prop))   marginal_prop <- type_plot
  if(missing(marginal_effect)) marginal_effect <- type_plot

  ## col_plot
  col_prop <- col_effect <- col
  col_prop[is.element(type_plot, marginal_prop) == FALSE] <- NA
  col_effect[is.element(type_plot, marginal_effect) == FALSE] <- NA

  ## Correct all esimates ----------
  p_coef_full <- p_high_full <- p_low_full <- c()
  p_name_full <- p_name_f_full <- p_col_full_effect <- p_col_full_prop <- p_cProp_full <- p_pch_full <- c()
  for(g in 1:length(factor_name)){
    p_cAME  <- came.out$cAME[[factor_name[g]]]

    if(plot_difference == "add"){p_cAME  <- p_cAME}
    if(plot_difference == "no"){p_cAME   <- p_cAME[(p_cAME$type %in% type_difference) == FALSE, ]}
    if(plot_difference == "only"){p_cAME  <- p_cAME[(p_cAME$type %in% type_difference) == TRUE, ]}

    p_coef <- c(NA, p_cAME$estimate)
    p_se   <- c(NA, p_cAME$se)
    p_high <- p_coef + 1.96*p_se
    p_low  <- p_coef - 1.96*p_se
    p_name_t <- paste(factor_name[g], ":   ", sep="")
    p_name_b   <- p_cAME$level
    p_name_b[p_cAME$type != p_cAME$type[1]] <- ""
    p_name_f <- c(p_name_t, rep("", length(p_name_b)))
    p_name   <- c("", p_name_b)
    p_col_prop  <- c(NA, rep(col_prop, length(unique(p_cAME$level))))
    p_col_effect  <- c(NA, rep(col_effect, length(unique(p_cAME$level))))
    p_pch  <- c(NA, rep(pch, length(unique(p_cAME$level))))

    ## Store values
    p_coef_full <- c(p_coef_full, p_coef)
    p_high_full <- c(p_high_full, p_high)
    p_low_full  <- c(p_low_full, p_low)
    p_name_f_full <- c(p_name_f_full, p_name_f)
    p_name_full <- c(p_name_full, p_name)
    p_col_full_prop  <- c(p_col_full_prop, p_col_prop)
    p_col_full_effect  <- c(p_col_full_effect, p_col_effect)
    p_pch_full  <- c(p_pch_full, p_pch)

    ## Store Distributions
    p_cProp <- came.out$cProp[[factor_name[g]]]

    if(plot_difference == "add"){p_cProp  <- p_cProp}
    if(plot_difference == "no"){p_cProp   <- p_cProp[(p_cProp$type %in% type_difference) == FALSE, ]}
    if(plot_difference == "only"){p_cProp  <- p_cProp[(p_cProp$type %in% type_difference) == TRUE, ]}

    p_cProp_full     <- c(p_cProp_full, NA, p_cProp$prop)
  }

  ## Add AME
  # AME_coef <- c(came.out$AME$estimate)
  # AME_high  <- c(came.out$AME$estimate + 1.96*came.out$AME$se)
  # AME_low   <- c(came.out$AME$estimate - 1.96*came.out$AME$se)
  # p_high_full <- c(NA, AME.high, p_high_full)
  # p_low_full <- c(NA, AME.low, p_low_full)
  p_name_f_full_prop <- p_name_f_full
  p_name_full_prop <- p_name_full
  # p_col_full_prop <- p_col_prop_full
  p_pch_full_prop  <- p_pch_full
  # p_name_f_AME <- c(rep("", length(came.out$AME$estimate)))
  # p_name_AME <- unique(came.out$AME$type)
  # p_col_AME <- c(col)
  # p_pch_AME  <- c(pch)

  ## Plot Setup for conditional AMEs ----------
  p_type <- unique(p_cAME$type)
  p_x <- c(min(p_low_full, na.rm=TRUE), max(p_high_full, na.rm = TRUE))

  ## Plot Setup for Proportion----------
  p_type_prop <- unique(p_cProp$type)
  if(came.out$difference == TRUE) p_x_prop <- c(min(p_cProp_full, na.rm=TRUE), max(p_cProp_full, na.rm = TRUE))
  if(came.out$difference == FALSE) p_x_prop <- c(0, max(p_cProp_full, na.rm = TRUE))

  ## Plot ----------
  if(plot_all == FALSE){
    par(oma = c(1, mar,1,1), mar=c(4,2,4,1))
    plot(rev(p_coef_full), seq(1:length(p_coef_full)), pch=rev(p_pch_full),
         yaxt="n", ylab = "", main = main, xlim = p_x,
         xlab = "Estimated Effects", col = rev(p_col_full_effect))
    segments(rev(p_low_full), seq(1:length(p_coef_full)),
             rev(p_high_full), seq(1:length(p_coef_full)),
             col = rev(p_col_full_effect))
    Axis(side=2, at = seq(1:length(p_name_full)),
         labels=rev(p_name_full), las=1, font = 1, tick=F, cex.axis = cex)
    Axis(side=2, at = seq(1:length(p_name_f_full)),
         labels=rev(p_name_f_full), las=1, font = 2, tick=F, cex.axis = cex)
    abline(v=0, lty=2)
    if(is.character(legend_pos[1])==TRUE) legend(legend_pos, p_type, col= col, pch = pch)
    if(is.character(legend_pos[1])==FALSE) legend(x=legend_pos[1], y=legend_pos[2],
                                                  p_type, col= col, pch = pch)

  }else if(plot_all == TRUE){
    par(mfrow=c(1,2), oma = c(1, mar,1,1), mar=c(4,2,4,1))
    ## Plot ----------
    plot(rev(p_cProp_full), seq(1:length(p_cProp_full)), pch=rev(p_pch_full_prop),
         yaxt="n", ylab = "", main = main[2], xlim = p_x_prop,
         xlab = "Proportions", col = rev(p_col_full_prop))
    segments(rev(p_cProp_full), seq(1:length(p_cProp_full)),
             0, seq(1:length(p_cProp_full)),
             col = rev(p_col_full_prop))
    abline(v=0, lty=2)
    Axis(side=2, at = seq(1:length(p_name_full)),
         labels=rev(p_name_full), las=1, font = 1, tick=F, cex.axis = cex)
    Axis(side=2, at = seq(1:length(p_name_f_full)),
         labels=rev(p_name_f_full), las=1, font = 2, tick=F, cex.axis = cex)
    if(is.character(legend_pos[1])==TRUE) legend(legend_pos, p_type_prop, col= col, pch = pch)
    if(is.character(legend_pos[1])==FALSE) legend(x=legend_pos[1], y=legend_pos[2],
                                                  p_type_prop, col= col, pch = pch)

    plot(rev(p_coef_full), seq(1:length(p_coef_full)), pch=rev(p_pch_full),
         yaxt="n", ylab = "", main = main[1], xlim = p_x,
         xlab = "Estimated Effects", col = rev(p_col_full_effect))
    segments(rev(p_low_full), seq(1:length(p_coef_full)),
             rev(p_high_full), seq(1:length(p_coef_full)),
             col = rev(p_col_full_effect))
    abline(v=0, lty=2)
    if(is.character(legend_pos[1])==TRUE) legend(legend_pos, p_type, col= col, pch = pch)
    if(is.character(legend_pos[1])==FALSE) legend(x=legend_pos[1], y=legend_pos[2],
                                                  p_type, col= col, pch = pch)
  }
  # {
  #   p_x_AME <- c(min(AME_low, na.rm=TRUE), max(AME_high, na.rm = TRUE))
  #
  #
  #   layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  #   par(oma = c(1,6,1,1), mar=c(4,2,4,1))
  #   ## AMEs
  #   plot(rev(AME_coef), seq(1:length(AME_coef)), pch=rev(p_pch_AME),
  #        yaxt="n", ylab = "", main = "AMEs", xlim = p_x_AME,
  #        xlab = "Estimated Effects", col = rev(p_col_AME))
  #   segments(rev(AME_low), seq(1:length(AME_coef)),
  #            rev(AME_high), seq(1:length(AME_coef)),
  #            col = rev(p_col_AME))
  #   Axis(side=2, at = seq(1:length(p_name_AME)), labels=rev(p_name_AME), las=1, font = 1)
  #   Axis(side=2, at = seq(1:length(p_name_f_AME)), labels=rev(p_name_f_AME), las=1, font = 2)
  #   abline(v=0, lty=2)
  #
  #   plot(rev(p_coef_full), seq(1:length(p_coef_full)), pch=rev(p_pch_full),
  #        yaxt="n", ylab = "", main = main[1], xlim = p_x,
  #        xlab = "Estimated Effects", col = rev(p_col_full))
  #   segments(rev(p_low_full), seq(1:length(p_coef_full)),
  #            rev(p_high_full), seq(1:length(p_coef_full)),
  #            col = rev(p_col_full))
  #   Axis(side=2, at = seq(1:length(p_name_full)), labels=rev(p_name_full), las=1, font = 1, tick=F)
  #   Axis(side=2, at = seq(1:length(p_name_f_full)), labels=rev(p_name_f_full), las=1, font = 2, tick=F)
  #   abline(v=0, lty=2)
  #   if(is.character(legend_pos[1])==TRUE) legend(legend_pos, p_type, col= col, pch = pch)
  #   if(is.character(legend_pos[1])==FALSE) legend(x=legend_pos[1], y=legend_pos[2],
  #                                                 p_type, col= col, pch = pch)
  #   ## Plot ----------
  #   plot(rev(p_cProp_full), seq(1:length(p_cProp_full)), pch=rev(p_pch_full_prop),
  #        yaxt="n", ylab = "", main = main[2], xlim = p_x_prop,
  #        xlab = "Proportions", col = rev(p_col_full_prop))
  #   segments(rev(p_cProp_full), seq(1:length(p_cProp_full)),
  #            0, seq(1:length(p_cProp_full)),
  #            col = rev(p_col_full_prop))
  #   abline(v=0, lty=2)
  #   if(is.character(legend_pos[1])==TRUE) legend(legend_pos, p_type_prop, col= col, pch = pch)
  #   if(is.character(legend_pos[1])==FALSE) legend(x=legend_pos[1], y=legend_pos[2],
  #                                                 p_type_prop, col= col, pch = pch)
  # }

}



#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param pair whether we use the paired-conjoint design
#' @param pair_id id for paired-conjoint design. required when 'pair = TRUE'
#' @export

decompose_AME <- function(outAME, factor_name, level_name, marginal_diff = c("Pop", "Exp")){

  marginal_dist <- outAME$input$marginal_dist
  marginal_dist_u_list  <- outAME$input$marginal_dist_u_list
  marginal_dist_u_base  <- outAME$input$marginal_dist_u_base
  marginal_type  <- outAME$input$marginal_type

  if(all(is.element(marginal_diff, marginal_type)) == FALSE){
    stop("`marginal_diff` can only take names in `marginal_type` specified in `AME_estimate_full.`")
  }

  boot_coef <- outAME$boot_coef
  use_name <- paste(factor_name, level_name, sep = "")

  # Estimate Marginal Contributions to AMEs ----------
  coef_focus <- boot_coef[, grep(use_name, colnames(boot_coef), fixed = T)]
  # vcov_focus <- vcovInt[grep(marginal_dist_u_base$level[m], rownames(vcovInt), fixed = T),
  #                       grep(marginal_dist_u_base$level[m], colnames(vcovInt), fixed = T)]
  if((length(coef_focus) > 0) == FALSE){
    stop(" 'level_name' cannot take the baseline level or undefined levels of the specified factor")
  }
  estNames <- gsub(paste(use_name, ":", sep = ""), "", colnames(coef_focus), fixed = T)
  estNames <- gsub(paste(":", colnames(coef_focus)[1], sep = ""), "", estNames, fixed = T)
  table_AME_diff <- c()
  # For each marginal distribution,
  ind1 <- which(marginal_type == marginal_diff[1])
  ind2 <- which(marginal_type == marginal_diff[2])
  marginal_dist_u1   <- marginal_dist_u_list[[ind1]]
  marginal_dist_u1$fac <- marginal_dist[[ind1]]$fac
  marginal_dist_u2   <- marginal_dist_u_list[[ind2]]
  marginal_dist_u2$fac <- marginal_dist[[ind2]]$fac

  # all(marginal_dist_u1$level == marginal_dist_u2$level) ## Check
  # all(marginal_dist_u1$fac == marginal_dist_u2$fac) ## Check

  marginal_factor <- setdiff(unique(marginal_dist_u1$fac), factor_name)
  base_mar1 <- marginal_dist_u1[match(estNames, marginal_dist_u1[, "level"]), ]
  base_mar2 <- marginal_dist_u2[match(estNames, marginal_dist_u2[, "level"]), ]

  for(m in 1:length(marginal_factor)){
    # Find weights

    use_mar1 <- base_mar1
    use_mar1[use_mar1$fac != marginal_factor[m], "prop"] <- 0
    use_mar2 <- base_mar2
    use_mar2[use_mar2$fac != marginal_factor[m], "prop"] <- 0

    coef_prop1 <- c(0, as.numeric(as.character(use_mar1[match(estNames, use_mar1[, "level"]), "prop"]))[-1])
    coef_prop2 <- c(0, as.numeric(as.character(use_mar2[match(estNames, use_mar2[, "level"]), "prop"]))[-1])
    # Compute AMEs
    coef_AME_diff <- mean(coef_focus %*% (coef_prop1 - coef_prop2))
    # coef_AME <- sum(coef_focus * coef_prop)
    se_AME_diff   <- sd(coef_focus %*% (coef_prop1 - coef_prop2))
    low_AME_diff    <- quantile(coef_focus %*% (coef_prop1 - coef_prop2), probs = c(0.025))
    high_AME_diff   <- quantile(coef_focus %*% (coef_prop1 - coef_prop2), probs = c(0.975))
    # se_AME <- sqrt(coef_prop%*%vcov_focus%*%coef_prop)
    AME_diff <- data.frame(matrix(NA, ncol = 0, nrow=1))
    AME_diff$type <- paste(marginal_diff[1], marginal_diff[2], sep = "-")
    AME_diff$factor   <- marginal_factor[m]
    AME_diff$estimate <- coef_AME_diff;
    AME_diff$se <- se_AME_diff
    AME_diff$low.95ci  <- low_AME_diff
    AME_diff$high.95ci <- high_AME_diff
    table_AME_diff <- rbind(table_AME_diff, AME_diff)
  }
  return(table_AME_diff)
}

coefMake <- function(original_level){
  # Expand Coefficients
  n_fac <- length(original_level)
  original_name <- list()

  # For main effects
  main_name <- ""
  for(z in 1:length(original_level)){
    original_name[[z]]  <- paste(names(original_level)[z], original_level[[z]], sep = "")
    main_name <- c(main_name, original_name[[z]][-1])
  }

  # For Interaction effects
  combMat <- combn(n_fac, 2)
  int_name <- c()
  for(z in 1:ncol(combMat)){
    L1 <- original_name[[combMat[1,z]]]
    L2 <- original_name[[combMat[2,z]]]
    c_1 <- seq(from = 2, to = length(L1))
    c_2 <- seq(from = 2, to = length(L2))
    int0 <- paste(L1[rep(c_1, times = length(c_2))],
                  L2[rep(c_2, each = length(c_1))], sep = ":")
    int_name <- c(int_name, int0)
  }
  coef_name <- c(main_name, int_name)
  return(coef_name)
}






