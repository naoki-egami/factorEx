# (Internal Functions) Help functions
AME.fit <- function(formula_full,
                    data,
                    pair = FALSE, cross_int = TRUE,
                    marginal_dist = NULL,
                    marginal_dist_u_list = NULL,
                    marginal_dist_u_base = NULL,
                    joint_dist_u_list = NULL,
                    marginal_type,
                    difference = FALSE,
                    three_way = FALSE,
                    est_AME = TRUE){
  # Differencing ----------
  if(pair==TRUE){
    data0 <- data[order(data$pair_id),]
    side <- rep(c(1,0),times=nrow(data0)/2)
    data1 <- data0[side==1,]
    data2 <- data0[side==0,]
    cluster_original <- data$cluster
    cluster <- data$cluster[side==1]

    # incorporate cross-profile interactions
    X1 <- model.matrix(formula_full, data=data1)
    ind_b <- attr(X1, "assign")
    X1  <- X1[,-1]
    X2 <- model.matrix(formula_full, data=data2)[ ,-1]
    X <- cbind(1, X1 - X2)

    if(cross_int == TRUE){

      base_fac <- all.vars(formula_full)[-1]
      data2_u <- data2[,base_fac]; colnames(data2_u) <- paste(base_fac,"_rp",sep="")
      data_c <- cbind(data1[,base_fac], data2_u)
      for_cross <- paste("~", paste(paste(base_fac, paste(base_fac,"_rp",sep=""), sep="*"),
                                    collapse = "+"))
      for_cross0 <- paste("~", paste(base_fac, collapse = "+"))
      data_cross0 <- model.matrix(as.formula(for_cross0), data = data_c)[,-1]
      sing <- 2*ncol(data_cross0)
      data_cross <- model.matrix(as.formula(for_cross), data = data_c)
      ind_b_c <- attr(data_cross, "assign")[-1]
      data_cross <- data_cross[,-1]
      X_cross <- data_cross[, c((sing + 1): ncol(data_cross))] # don't need to rename
      ind_b_c <- ind_b_c[c((sing + 1): ncol(data_cross))]

      # modify X and ind_b
      X <- cbind(X, X_cross)
      ind_b <- c(ind_b, ind_b_c - min(ind_b_c) + 1 + max(ind_b))
    }
    y <- model.frame(formula_full, data=data1)[ ,1]
    # base_name <- c("(Intercept)", colnames(X1))

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
  coefInt <- coef(main_lm)  ## this includes intercept
  coefInt[is.na(coefInt) == TRUE] <- 0 ## Insert 0 to non-existing pairs
  base_name <- sub("X", "", names(coefInt))
  names(coefInt) <- base_name

  # Estimate AMEs ----------
  # Estimeate AME from two-ways (or three-ways)
  if(est_AME ==  TRUE){
    table_AME <- coefIntAME(coefInt = coefInt, vcovInt = NULL, SE = FALSE,
                            marginal_dist = marginal_dist, marginal_dist_u_list = marginal_dist_u_list,
                            marginal_dist_u_base = marginal_dist_u_base, marginal_type = marginal_type,
                            difference = difference, cross_int = cross_int, three_way = three_way,
                            joint_dist_u_list = joint_dist_u_list)
    out <- list("table_AME" = table_AME, "coef" = coefInt, "ind_b" = ind_b)
  }else if(est_AME ==  FALSE){
    out <- list("coef" = coefInt, "ind_b" = ind_b)
  }
  return(out)
}

coefIntAME <- function(coefInt, vcovInt, SE = FALSE, marginal_dist, marginal_dist_u_list, marginal_dist_u_base,
                       marginal_type, difference, cross_int, three_way = FALSE, joint_dist_u_list = NULL){

  # Estimate AMEs ----------
  table_AME <- c()
  for(m in 1:nrow(marginal_dist_u_base)){
    coef_focus <- coefInt[grep(marginal_dist_u_base$level[m], names(coefInt), fixed = T)]
    if(SE == TRUE){
      vcov_focus <- vcovInt[grep(marginal_dist_u_base$level[m], rownames(vcovInt), fixed = T),
                            grep(marginal_dist_u_base$level[m], colnames(vcovInt), fixed = T)]
    }
    if(length(coef_focus) > 0){
      estNames <- gsub(paste(marginal_dist_u_base$level[m], ":", sep = ""), "", names(coef_focus), fixed = T)
      estNames <- gsub(paste(":", names(coef_focus)[1], sep = ""), "", estNames, fixed = T)

      if(cross_int == TRUE){
        estNames <- sub(paste(marginal_dist[[1]]$factor[m],"_rp", sep = ""),
                        marginal_dist[[1]]$factor[m], estNames)
      }
      # allow for three-ways
      if(three_way == TRUE){
        estL <- strsplit(estNames, ":")
        estL <- do.call("rbind", lapply(estL, function(x) if(length(x)==1){c(x[1], NA)}else{x}))
      }

      table_AME_m <- c()
      # For each marginal distribution,
      for(z in 1:length(marginal_dist_u_list)){
        marginal_dist_u <- marginal_dist_u_list[[z]]
        joint_dist_u <- joint_dist_u_list[[z]]
        # Find weights
        if(three_way == FALSE){
          coef_prop <- c(1, as.numeric(as.character(marginal_dist_u[match(estNames, marginal_dist_u[, "level"]), "prop"]))[-1])
        }else if(three_way == TRUE){
          # allow for three-ways
          ind_two  <- which(apply(estL, 1, function(x) any(is.na(x))))[-1]
          prop_two <- marginal_dist_u[match(estL[ind_two,1], marginal_dist_u[, "level"]), "prop"]
          ind_three  <- which(apply(estL, 1, function(x) all(is.na(x) == FALSE)))
          prop_three <- joint_dist_u[row.match(as.data.frame(matrix(estL[ind_three, ], ncol = 2)),
                                               as.data.frame(joint_dist_u[, 1:2])), "prop"]
          coef_prop  <- rep(NA, nrow(estL))
          coef_prop[1] <- 1
          coef_prop[ind_two]   <-  prop_two
          coef_prop[ind_three] <- prop_three
        }

        # Compute AMEs
        coef_AME <- sum(coef_focus * coef_prop)
        if(SE == TRUE) se_AME <- sqrt(coef_prop%*%vcov_focus%*%coef_prop)
        AME <- data.frame(matrix(NA, ncol = 0, nrow=1))
        AME$type <- marginal_type[z]
        AME$factor   <- marginal_dist[[z]][m,1]; AME$level <- marginal_dist[[z]][m,2]
        AME$estimate <- coef_AME;
        if(SE == TRUE) AME$se <- se_AME
        table_AME_m <- rbind(table_AME_m, AME)
      }
      if(difference == TRUE){
        for(z in 2:length(marginal_dist_u_list)){
          marginal_dist_u <- marginal_dist_u_list[[z]]
          joint_dist_u <- joint_dist_u_list[[z]]
          joint_dist_u_base <- joint_dist_u_list[[1]]
          # Find weights
          if(three_way == FALSE){
            coef_prop <- c(1, as.numeric(as.character(marginal_dist_u[match(estNames, marginal_dist_u[, "level"]), "prop"]))[-1])
            coef_prop0 <- c(1, as.numeric(as.character(marginal_dist_u_base[match(estNames, marginal_dist_u_base[, "level"]), "prop"]))[-1])
          }else if(three_way == TRUE){

            ind_two  <- which(apply(estL, 1, function(x) any(is.na(x))))[-1]
            ind_three  <- which(apply(estL, 1, function(x) all(is.na(x)==FALSE)))

            prop_two <- marginal_dist_u[match(estL[ind_two,1], marginal_dist_u[, "level"]), "prop"]
            prop_three <- joint_dist_u[row.match(as.data.frame(estL[ind_three, ]), as.data.frame(joint_dist_u[, 1:2])), "prop"]
            coef_prop  <- rep(NA, nrow(estL))
            coef_prop[1] <- 1
            coef_prop[ind_two]   <-  prop_two
            coef_prop[ind_three] <- prop_three

            prop_two0 <- marginal_dist_u_base[match(estL[ind_two,1], marginal_dist_u_base[, "level"]), "prop"]
            prop_three0 <- joint_dist_u_base[row.match(as.data.frame(estL[ind_three, ]), as.data.frame(joint_dist_u_base[, 1:2])), "prop"]
            coef_prop0  <- rep(NA, nrow(estL))
            coef_prop0[1] <- 1
            coef_prop0[ind_two]   <-  prop_two0
            coef_prop0[ind_three] <- prop_three0
          }

          # Compute AMEs
          coef_prop_d <- (coef_prop - coef_prop0)
          coef_AME_dif <- sum(coef_focus * coef_prop_d)
          if(SE == TRUE) se_AME_dif <- sqrt(coef_prop_d%*%vcov_focus%*%coef_prop_d)
          AME_dif <- data.frame(matrix(NA, ncol = 0, nrow=1))
          AME_dif$type <- paste(marginal_type[z],"-",marginal_type[1],sep="")
          AME_dif$factor   <- marginal_dist[[z]][m,1]; AME_dif$level <- marginal_dist[[z]][m,2]
          AME_dif$estimate <- coef_AME_dif;
          if(SE == TRUE) AME_dif$se <- se_AME_dif
          table_AME_m <- rbind(table_AME_m, AME_dif)
        }
      }
      table_AME <- rbind(table_AME, table_AME_m)
    }
  }
  colnames(table_AME) <- c("type", "factor", "level", "estimate")
  if(SE == TRUE) colnames(table_AME) <- c("type", "factor", "level", "estimate", "se")
  return(table_AME)
}

AME.fit.STD.sep <- function(formula,
                            data,
                            pair = FALSE,
                            marginal_dist,
                            marginal_dist_u_base){
  # Differencing ----------
  if(pair==TRUE){
    data <- data[order(data$pair_id),]
  }
  cluster_original <- data$cluster
  y <- model.frame(formula, data=data)[,1]
  fac_use <- all.vars(formula)[-1]
  coefInt <- 0

  # Fit the model separately for each factor
  for(z in 1:length(fac_use)){
    for_use <- as.formula(paste("~ ", fac_use[z], sep = ""))
    X <- model.matrix(for_use, data = data)
    lm_use <- lm(y ~ X)
    coefSub <- summary(lm_use)$coef[-1,1]
    coefSub <- coefSub[is.na(coefSub) == FALSE]
    base_name <- colnames(X)[-1]
    names(coefSub) <- base_name
    coefInt <- c(coefInt, coefSub)
  }

  # Estimate AMEs ----------
  table_AME <- c()
  for(m in 1:nrow(marginal_dist_u_base)){
    coef_focus <- coefInt[grep(marginal_dist_u_base$level[m], names(coefInt), fixed = T)]
    if(length(coef_focus) > 0){
      estNames <- gsub(paste(marginal_dist_u_base$level[m], ":", sep = ""), "", names(coef_focus), fixed = T)
      estNames <- gsub(paste(":", names(coef_focus)[1], sep = ""), "", estNames, fixed = T)

      table_AME_m <- data.frame(matrix(NA, ncol = 0, nrow=1))
      table_AME_m$type <- "sample AMCE"
      table_AME_m$factor   <- marginal_dist[[1]][m,1]; table_AME_m$level <- marginal_dist[[1]][m,2]
      table_AME_m$estimate <- coef_focus;
      table_AME <- rbind(table_AME, table_AME_m)
    }
  }
  colnames(table_AME) <- c("type", "factor", "level", "estimate")
  return(table_AME)
}

AME.fit.STD.se.sep <- function(formula,
                           data,
                           pair = FALSE,
                           marginal_dist,
                           marginal_dist_u_base){
  # Differencing ----------
  if(pair==TRUE){
    data <- data[order(data$pair_id),]
  }
  cluster_original <- data$cluster
  y <- model.frame(formula, data=data)[,1]
  fac_use <- all.vars(formula)[-1]
  coefInt <- 0
  vcovInt <- 0

  # Fit the model separately for each factor
  for(z in 1:length(fac_use)){
    for_use <- as.formula(paste("~ ", fac_use[z], sep = ""))
    X <- model.matrix(for_use, data = data)
    lm_use <- lm(y ~ X)
    coefSub <- summary(lm_use)$coef[-1,1]
    coefSub <- coefSub[is.na(coefSub) == FALSE]
    base_name <- colnames(X)[-1]
    names(coefSub) <- base_name
    coefInt <- c(coefInt, coefSub)
    vcovSub <- diag(cluster_se_glm(lm_use, cluster = as.factor(data$cluster)))[-1]
    names(vcovSub) <- base_name
    vcovInt <- c(vcovInt, vcovSub)
  }

  # Estimate AMEs ----------
  table_AME <- c()
  for(m in 1:nrow(marginal_dist_u_base)){
    coef_focus <- coefInt[grep(marginal_dist_u_base$level[m], names(coefInt), fixed = T)]
    vcov_focus <- vcovInt[grep(marginal_dist_u_base$level[m], names(coefInt), fixed = T)]
    if(length(coef_focus) > 0){
      estNames <- gsub(paste(marginal_dist_u_base$level[m], ":", sep = ""), "", names(coef_focus), fixed = T)
      estNames <- gsub(paste(":", names(coef_focus)[1], sep = ""), "", estNames, fixed = T)

      table_AME_m <- data.frame(matrix(NA, ncol = 0, nrow=1))
      table_AME_m$type <- "sample AMCE"
      table_AME_m$factor   <- marginal_dist[[1]][m,1]; table_AME_m$level <- marginal_dist[[1]][m,2]
      table_AME_m$estimate <- coef_focus; table_AME_m$se <- sqrt(vcov_focus)
      table_AME <- rbind(table_AME, table_AME_m)
    }
  }
  colnames(table_AME) <- c("type", "factor", "level", "estimate", "se")
  return(table_AME)
}

Fthree <- function(formula,
                   formula_three_c,
                   data,
                   pair = FALSE, cross_int = TRUE,
                   coef_f){


  factor_l <- length(all.vars(formula)[-1])
  combMat <- combn(factor_l,2); intNames <- c()
  for(k in 1:ncol(combMat)){
    intNames[k] <- paste(all.vars(formula)[-1][combMat[1,k]], "*", all.vars(formula)[-1][combMat[2,k]], sep = "")
  }
  formula_full0 <- paste(all.vars(formula)[1], "~", paste(intNames, collapse = "+"), sep="")
  if(is.null(formula_three_c) == TRUE){
    formula_full <- as.formula(formula_full0)
  }else{
    formula_full <- as.formula(paste(formula_full0, formula_three_c, sep = "+"))
  }

  # Differencing ----------
  if(pair==TRUE){
    data0 <- data[order(data$pair_id),]
    side <- rep(c(1,0),times=nrow(data0)/2)
    data1 <- data0[side==1,]
    data2 <- data0[side==0,]

    # incorporate cross-profile interactions
    X1 <- model.matrix(formula_full, data=data1)[,-1]
    X2 <- model.matrix(formula_full, data=data2)[ ,-1]
    X <- cbind(1, X1 - X2)

    if(cross_int == TRUE){

      base_fac <- all.vars(formula_full)[-1]
      data2_u <- data2[,base_fac]; colnames(data2_u) <- paste(base_fac,"_rp",sep="")
      data_c <- cbind(data1[,base_fac], data2_u)
      for_cross <- paste("~", paste(paste(base_fac, paste(base_fac,"_rp",sep=""), sep="*"),
                                    collapse = "+"))
      for_cross0 <- paste("~", paste(base_fac, collapse = "+"))
      data_cross0 <- model.matrix(as.formula(for_cross0), data = data_c)[,-1]
      sing <- 2*ncol(data_cross0)
      data_cross <- model.matrix(as.formula(for_cross), data = data_c)
      data_cross <- data_cross[,-1]
      X_cross <- data_cross[, c((sing + 1): ncol(data_cross))] # don't need to rename

      # modify X and ind_b
      X <- cbind(X, X_cross)
    }
    y <- model.frame(formula_full, data=data1)[ ,1]
    # base_name <- c("(Intercept)", colnames(X1))

  }else{
    X <- model.matrix(formula_full, data=data)
    y <- model.frame(formula_full, data=data)[,1]
    # base_name <- colnames(X)
    side <- NULL
  }

  # Approximate F-test
  coef_order <- unlist(lapply(strsplit(names(coef_f), ":"), length))
  ind_three  <- which(coef_order == 3)
  R <-matrix(0, ncol = length(coef_f), nrow = length(ind_three))
  R[seq(1:length(ind_three)), ind_three] <- diag(length(ind_three))
  sigma2 <- sum((y - X%*%coef_f)^2)/(nrow(X) - ncol(X))
  Ft_n <- t(R%*%coef_f)%*%solve(R%*%solve(t(X)%*%X)%*%t(R))%*%(R%*%coef_f)
  Ft_d <- sigma2*nrow(X)
  Ft <- Ft_n/Ft_d
  pvalue <- pf(Ft, df1 = nrow(R), df2 = nrow(X) - ncol(X), lower.tail = FALSE)

  return(pvalue)
}


# Cluster standard errors
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

coefMake <- function(original_level, cross_int, formula_three_c){
  # Expand Coefficients
  n_fac <- length(original_level)
  original_name <- list()

  # For main effects
  main_name <- ""
  for(z in 1:length(original_level)){
    original_name[[z]]  <- paste(names(original_level)[z], original_level[[z]], sep = "")
    main_name <- c(main_name, original_name[[z]][-1])
  }

  # For Interaction effects (within profiles)
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

  # For three-way interactions
  if(is.null(formula_three_c) == FALSE){
      bterm <- terms(as.formula(paste("~", formula_three_c, sep = "")))
      bterm2 <- attr(bterm, "factors")[, attr(bterm, "order") == 3]
      row_num <- match(rownames(attr(bterm, "factors")), names(original_level))
      combMat <- row_num*bterm2
      if(is.matrix(combMat)==TRUE) combMat <- apply(combMat, 2, function(x) sort(x[x!=0]))
      else{ combMat <- matrix(sort(combMat), ncol = 1, nrow = 3)}

      int3_name <- c()
      for(z in 1:ncol(combMat)){
        L1 <- original_name[[combMat[1,z]]]
        L2 <- original_name[[combMat[2,z]]]
        L3 <- original_name[[combMat[3,z]]]
        c_1 <- seq(from = 2, to = length(L1))
        c_2 <- seq(from = 2, to = length(L2))
        c_3 <- seq(from = 2, to = length(L3))
        int0 <- paste(L1[rep(c_1, times = (length(c_2)*length(c_3)))],
                      rep(L2[rep(c_2, each = length(c_1))], times=length(c_3)),
                      L3[rep(c_3, each = (length(c_1)*length(c_2)))],
                      sep = ":")
        int3_name <- c(int3_name, int0)
      }
    coef_name <- c(coef_name, int3_name)
  }

  # For Interaction effects (across profiles)
  if(cross_int == TRUE){
    cross_int_name <- c()
    for(z in 1:length(original_level)){
      L1 <- original_name[[z]]
      L2  <- paste(paste(names(original_level)[z],"_rp",sep=""),
                                   original_level[[z]], sep = "")
      c_1 <- c_2 <- seq(from = 2, to = length(L1))
      int0 <- paste(L1[rep(c_1, times = length(c_2))],
                    L2[rep(c_2, each = length(c_1))], sep = ":")
      cross_int_name <- c(cross_int_name, int0)
    }
    coef_name <- c(coef_name, cross_int_name)
  }
  return(coef_name)
}

createDistBase <- function(formula, data, marginal){
    dataX <- model.frame(formula, data)[,-1]

    all_levels <- lapply(dataX, levels)

    if(marginal == TRUE){
      factor_u <-  rep(names(all_levels), unlist(lapply(all_levels, length)))
      levels_u <- c(unlist(all_levels))
      prop_u   <- unlist(lapply(dataX, function(x) prop.table(table(x))))

      target <- data.frame(matrix(NA, ncol = 0, nrow = length(factor_u)))
      target$factor   <- factor_u
      target$levels <- levels_u
      target$prop   <- prop_u
      attributes(target)$dist_type <- "marginal"

    }else if(marginal == FALSE){
      fac_name <- names(all_levels)
      combMat  <- combn(length(fac_name), 2)
      target   <- list()
      for(i in 1:ncol(combMat)){
        levels_u <- expand.grid(all_levels[[combMat[1, i]]], all_levels[[combMat[2, i]]])
        factor_u <- cbind(rep(fac_name[combMat[1, i]], nrow(levels_u)), rep(fac_name[combMat[2, i]], nrow(levels_u)))
        target_0 <- as.data.frame(cbind(factor_u, levels_u))
        b_tab <- table(dataX[, c(fac_name[combMat[1, i]], fac_name[combMat[2, i]])])
        b_tab <- b_tab[match(all_levels[[combMat[1, i]]], rownames(b_tab)),]
        b_tab <- b_tab[, match(all_levels[[combMat[2, i]]], colnames(b_tab))]
        prop_u   <- c(b_tab)/sum(b_tab)
        target_0$prop   <- prop_u
        colnames(target_0) <- c("factor_1", "factor_2", "levels_1", "levels_2", "prop")
        target[[i]] <- target_0
        names(target)[[i]] <- paste(fac_name[combMat[1, i]], ":", fac_name[combMat[2, i]], sep = "")
      }
      attributes(target)$dist_type <- "joint"
    }
    return(target)
}

createDist0 <- function(formula, data, exp_data, marginal){

  formula_use <- formula(paste("~",  as.character(formula)[3], sep  = ""))

  col_in <- is.element(all.vars(formula_use), colnames(data))

  if(all(col_in) == FALSE){
    formula_use <- formula(paste("~", paste(all.vars(formula_use)[col_in == TRUE], collapse = "+"), sep = ""))
  }
  dataX <- model.frame(formula_use, data)

  all_levels <- lapply(dataX, levels)

  if(marginal == TRUE){
    target   <- lapply(dataX, function(x) prop.table(c(table(x))))

    attributes(target)$dist_type <- "marginal"
    attributes(target)$factor <- names(target)
    attributes(target)$level  <- unlist(all_levels)

  }else if(marginal == FALSE){
    fac_name <- names(all_levels)
    combMat  <- combn(length(fac_name), 2)
    target   <- list()
    for(i in 1:ncol(combMat)){
      levels_u <- expand.grid(all_levels[[combMat[1, i]]], all_levels[[combMat[2, i]]])
      factor_u <- cbind(rep(fac_name[combMat[1, i]], nrow(levels_u)), rep(fac_name[combMat[2, i]], nrow(levels_u)))
      target_0 <- as.data.frame(cbind(factor_u, levels_u))
      b_tab <- table(dataX[, c(fac_name[combMat[1, i]], fac_name[combMat[2, i]])])
      b_tab <- b_tab[match(all_levels[[combMat[1, i]]], rownames(b_tab)),]
      b_tab <- b_tab[, match(all_levels[[combMat[2, i]]], colnames(b_tab))]
      prop_u   <- c(b_tab)/sum(b_tab)
      target_0$prop   <- prop_u
      colnames(target_0) <- c("factor_1", "factor_2", "levels_1", "levels_2", "prop")
      target[[i]] <- target_0
      names(target)[[i]] <- paste(fac_name[combMat[1, i]], ":", fac_name[combMat[2, i]], sep = "")
    }
    attributes(target)$dist_type <- "joint"
    attributes(target)$factor <- fac_name
    attributes(target)$level  <- unlist(all_levels)
  }
  return(target)
}

createDist <- function(formula, target_data, exp_data, type = "marginal"){
  # Create Distribution based target_data and augment with Experimental Data

  if(type == "marginal"){
    dist <- createDist0(formula = formula, data = target_data, exp_data = exp_data, marginal = TRUE)

    exp_marginal <- createDist0(formula, data = exp_data, exp_data = exp_data, marginal = TRUE)
    if(all(is.element(names(exp_marginal), names(dist))) ==  TRUE ){
      return(dist)
    }else{
      dist_new <- exp_marginal
      for(i in 1:length(exp_marginal)){
        if(names(exp_marginal)[i] %in% names(dist)){
          dist_new[[i]] <-  dist[[match(names(exp_marginal)[i], names(dist))]]
        }
      }
      return(dist_new)
    }
  }


  if(type == "joint"){
    target_marginal <- createDist0(formula = formula, data = target_data, exp_data = exp_data, marginal = TRUE)
    dist <- createDist0(formula = formula, data = target_data, exp_data = exp_data, marginal = FALSE)
    exp_joint <- createDist0(formula, data = exp_data, exp_data  = exp_data, marginal = FALSE)

    if(all(is.element(names(exp_joint), names(dist))) ==  TRUE){
      return(dist)
    }else{
      exp_fac_name  <- attributes(exp_joint)$factor
      dist_fac_name <- attributes(dist)$factor
      missing_fac <- which(is.element(exp_fac_name, dist_fac_name) == FALSE)
      combMat <- combn(seq(1:length(attributes(exp_joint)$factor)), 2)
      combMat_m <- matrix(combMat %in% missing_fac, nrow = 2)
      combMat_m_ind <- apply(combMat_m, 2, function(x) any(x == TRUE))
      dist_new <- exp_joint
      exp_marginal <- createDist0(formula, data = exp_data, exp_data  = exp_data, marginal = TRUE)
      for(i in 1:length(exp_joint)){
        if(combMat_m_ind[i] == FALSE){ # none is missing
          dist_new[[i]] <-  dist[[match(names(exp_joint)[i], names(dist))]]
        }else{
          # if(all(combMat_m[1:2, i] == TRUE))   # both are missing, just use exp_joint
          if(all(combMat_m[1:2, i]) ==  FALSE){  # when only one is missing
            if(combMat_m[1, i] == TRUE){ # the first is missing
              prop_1 <-  exp_marginal[[combMat[1, i]]]
              prop_2 <-  target_marginal[[match(exp_fac_name[combMat[2, i]],dist_fac_name)]]
            }else if(combMat_m[2, i] == TRUE){ # the second is missing
              prop_2 <-  exp_marginal[[combMat[2, i]]]
              prop_1 <-  target_marginal[[match(exp_fac_name[combMat[1, i]],dist_fac_name)]]
            }
            prop_new <- rep(prop_1, times = length(prop_2))*rep(prop_2, each = length(prop_1))
            dist_new[[i]]$prop <- prop_new
          }
        }
      }
      return(dist_new)
    }
  }
}


Joint2Marginal <- function(joint_dist){

  prop_u_1 <- prop_u_0  <- list()
  fac_name_1 <- fac_name_0 <- c()
  for(i in 1:length(joint_dist)){
    prop_u_0[[i]] <- tapply(joint_dist[[i]]$prop, joint_dist[[i]]$levels_1, sum)
    prop_u_1[[i]] <- tapply(joint_dist[[i]]$prop, joint_dist[[i]]$levels_2, sum)
    fac_name_0[i] <- as.character(joint_dist[[i]]$factor_1[1])
    fac_name_1[i] <- as.character(joint_dist[[i]]$factor_2[1])
  }
  fac_name <- c(fac_name_0, fac_name_1)
  prop_u   <- append(prop_u_0, prop_u_1)
  prop_u2  <- prop_u[match(unique(fac_name), fac_name)]
  fac_name2  <- fac_name[match(unique(fac_name), fac_name)]

  names(prop_u2) <- fac_name2
  target <- prop_u2
  attributes(target)$dist_type <- "marginal"

  return(target)
}

Marginal2Joint <- function(marginal_dist){
  fac_name  <- names(marginal_dist)
  all_levels <- marginal_dist
  level_names <- lapply(marginal_dist, names)

  for(i in 1:length(all_levels)){
    names(all_levels[[i]]) <- level_names[[i]]
  }
  combMat  <- combn(length(fac_name), 2)
  target   <- list()
  for(i in 1:ncol(combMat)){
    levels_u <- expand.grid(level_names[[combMat[1, i]]], level_names[[combMat[2, i]]])
    factor_u <- cbind(rep(fac_name[combMat[1, i]], nrow(levels_u)), rep(fac_name[combMat[2, i]], nrow(levels_u)))
    target_0 <- as.data.frame(cbind(factor_u, levels_u))
    prop_u_0 <- apply(expand.grid(all_levels[[combMat[1, i]]], all_levels[[combMat[2, i]]]), 1, prod)
    prop_u   <- apply(expand.grid(all_levels[[combMat[1, i]]], all_levels[[combMat[2, i]]]), 1, prod)
    target_0$prop   <- prop_u
    colnames(target_0) <- c("factor_1", "factor_2", "levels_1", "levels_2", "prop")
    target[[i]] <- target_0
    names(target)[[i]] <- paste(fac_name[combMat[1, i]], ":", fac_name[combMat[2, i]], sep = "")
  }
  attributes(target)$dist_type <- "joint"

  return(target)
}


checkDist <- function(check_dist, type = "marginal", formula, data){

  if(type == "marginal"){

    exp_marginal <- createDist(formula = formula, target_data = data, exp_data = data, type = "marginal")

    if(setequal(names(check_dist), names(exp_marginal)) == FALSE){
      stop(" 'factor' of 'target_dist' should one of factors in 'data' ")
    }
    if(setequal(names(unlist(check_dist)), names(unlist(exp_marginal))) == FALSE){
      stop(" 'levels' of 'target_dist' should one of levels in 'data' ")
    }

    # reorder factor and levels
    check_dist <- check_dist[match(names(exp_marginal), names(check_dist))]
    for(i in 1:length(check_dist)){
      check_dist[[i]] <- check_dist[[i]][match(names(exp_marginal[[i]]), names(check_dist[[i]]))]
    }
  }

  if(type == "joint"){
    # check factor and level names
    exp_marginal <- createDist(formula = formula, target_data = data, exp_data = data, type = "marginal")
    exp_joint <- createDist(formula = formula, target_data = data, exp_data = data, type = "joint")
    check_dist_name <- list()
    exp_dist_name <- list()

    if(length(check_dist) != length(exp_joint)){
      stop(" When target_type = joint, please specify the two dimensional joint distributions for all pairs between factors")
    }

    for(i in 1:length(check_dist)){
      if(all(colnames(check_dist[[i]]) == c("factor_1", "factor_2", "levels_1", "levels_2", "prop")) == FALSE){
        stop(" When target_type = joint, 'colnames' in each element of target_dist should be c('factor_1', 'factor_2', 'levels_1', `levels_2`, 'prop')")
      }

      if(all(is.element(c(as.character(check_dist[[i]][, c("factor_1")]),
                          as.character(check_dist[[i]][, c("factor_2")])), attributes(exp_marginal)$factor)) == FALSE){
        stop(" 'factor_1' and  'factor_2' in each element of target_dist should one of factors in 'data' ")
      }
      if(all(is.element(c(as.character(check_dist[[i]][, c("levels_1")]),
                          as.character(check_dist[[i]][, c("levels_2")])), attributes(exp_marginal)$level)) == FALSE){
        stop(" 'levels_1' and  'levels_2' in each element of target_dist should one of levels in 'data' ")
      }
      check_dist_name[[i]] <- c(unique(as.character(check_dist[[i]][, c("factor_1")])), unique(as.character(check_dist[[i]][, c("factor_2")])))
      exp_dist_name[[i]] <- c(unique(as.character(exp_joint[[i]][, c("factor_1")])), unique(as.character(exp_joint[[i]][, c("factor_2")])))
    }
    # reorder factor-interaction
    check_reorder <- c()
    for(i in 1:length(check_dist)){
      check_reorder[i]  <- which(unlist(lapply(check_dist_name, function(x) setequal(x, exp_dist_name[[i]]))))
    }
    check_dist  <- check_dist[check_reorder]
    names(check_dist) <-  names(exp_joint)

    # Finally reorder factors and levels within each element
    for(i in 1:length(check_dist)){
      exp_col <- c(as.character(exp_joint[[i]][1, c("factor_1")]), as.character(exp_joint[[i]][1, c("factor_2")]))
      check_col <- c(as.character(check_dist[[i]][1, c("factor_1")]), as.character(check_dist[[i]][1, c("factor_2")]))
      col_ind <- match(exp_col, check_col)
      check_dist[[i]]  <- cbind(check_dist[[i]][, c(1, 2)[col_ind]], check_dist[[i]][, c(3, 4)[col_ind]], check_dist[[i]][, 5])
      colnames(check_dist[[i]]) <- c("factor_1", "factor_2", "levels_1", "levels_2",  "prop")

      # Finally reorder factor and levels
      check_dist[[i]] <- check_dist[[i]][row.match(exp_joint[[i]][, c("factor_1", "factor_2", "levels_1", "levels_2")],
                                                   check_dist[[i]][, c("factor_1", "factor_2", "levels_1", "levels_2")]),]

      check_dist[[i]]$factor_1 <- as.character(check_dist[[i]]$factor_1)
      check_dist[[i]]$factor_2 <- as.character(check_dist[[i]]$factor_2)
      check_dist[[i]]$levels_1 <- as.character(check_dist[[i]]$levels_1)
      check_dist[[i]]$levels_2 <- as.character(check_dist[[i]]$levels_2)
    }
  }

  if(type == "target_data"){
    if(all(is.element(colnames(check_dist), all.vars(formula))) == FALSE){
      stop(" When target_type = target_data, 'colnames' in each element of target_dist should be one of the factors in formula")
    }
    new_marginal <- createDist(formula = formula, target_data = check_dist, exp_data = data, type = "marginal")
    exp_marginal <- createDist(formula = formula, target_data = data, exp_data = data, type = "marginal")

    if(all(is.element(attributes(new_marginal)$level, attributes(exp_marginal)$level)) == FALSE){
      stop(" 'levels' used in 'target_dist' should one of levels in 'data' ")
    }
  }

  return(check_dist)
}

list2data <- function(marginal_dist){
  levels_u <- unlist(lapply(marginal_dist, names))
  factor_u <-  rep(names(marginal_dist), unlist(lapply(marginal_dist, length)))

  target <- data.frame(matrix(NA, ncol = 0, nrow = length(factor_u)))
  target$factor   <- factor_u
  target$levels <- levels_u
  target$prop   <- unlist(marginal_dist)
  return(target)
}

prepare_data <- function(formula, data, marginal_dist,
                         original_level, collapse_level){
  # Task 1: Prepare Data

  # Housekeeping
  factor_use <- all.vars(formula)[-1]
  n_fac <- length(factor_use)

  # New Level Names
  collapse_u <- list()
  for(i in 1:n_fac){
    original_level_use <- original_level[[i]]
    new_name <- c()
    for(j in 1:length(unique(collapse_level[[i]]))){
      new_name[j] <- paste(original_level_use[collapse_level[[i]]==j], collapse="/")
    }
    collapse_u[[i]] <- new_name[collapse_level[[i]]]
  }

  # Prepare New Data
  data_new <- data

  # Releveling Data
  for(i in 1:n_fac){
    levels(data_new[, factor_use[i]]) <- collapse_u[[i]]
    data_new[,factor_use[i]] <- factor(data_new[,factor_use[i]], ordered=FALSE)
  }

  # Task 2: Prepare Marginal Distributions
  ## subsetting
  marginal_dist_new <- list()
  for(z in 1:length(marginal_dist)){
    marginal_dist_use <- marginal_dist[[z]]

    # Reorder based on the original levels
    new_mar <- c()
    for(i in 1:n_fac){
      base <- marginal_dist_use[marginal_dist_use$factor == factor_use[i], ]
      base <- base[match(original_level[[i]], base$levels), ]
      prop <- tapply(base$prop, collapse_level[[i]], sum)
      levels <- unique(collapse_u[[i]])
      factor <- rep(unique(base$factor), length(levels))
      base2 <- as.data.frame(cbind(factor, levels))
      base2$factor <- as.character(base2$factor)
      base2$levels <- as.character(base2$levels)
      base2$prop <- prop
      new_mar <- rbind(new_mar, base2)
    }
    marginal_dist_new[[z]] <- rbind(new_mar,
                                    marginal_dist_use[marginal_dist_use$factor %in% factor_use == FALSE,])
  }
  names(marginal_dist_new) <- names(marginal_dist)

  # output
  out <- list("data_new" = data_new, "marginal_dist_new" = marginal_dist_new)
  return(out)
}

