#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param pair whether we use the paired-conjoint design
#' @param pair_id id for paired-conjoint design. required when 'pair = TRUE'
#' @importFrom sandwich sandwich estfun
#' @export

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
  # coefInt <- coefInt[is.na(coefInt) == FALSE]
  coefInt[is.na(coefInt) == TRUE] <- 0 ## Insert 0 to non-existing pairs
  base_name <- sub("X", "", names(coefInt))
  # base_name <- gsub(" ", "", base_name)
  names(coefInt) <- base_name
  # vcovInt <- vcovCR(main_lm, cluster = as.factor(cluster), type = "CR2")
  # colnames(vcovInt) <- rownames(vcovInt) <- base_name

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

  # table_AME <- c()
  # for(m in 1:nrow(marginal_dist_u_base)){
  #   coef_focus <- coefInt[grep(marginal_dist_u_base$level[m], names(coefInt), fixed = T)]
  #   # vcov_focus <- vcovInt[grep(marginal_dist_u_base$level[m], rownames(vcovInt), fixed = T),
  #   #                       grep(marginal_dist_u_base$level[m], colnames(vcovInt), fixed = T)]
  #   if(length(coef_focus) > 0){
  #     estNames <- gsub(paste(marginal_dist_u_base$level[m], ":", sep = ""), "", names(coef_focus), fixed = T)
  #     estNames <- gsub(paste(":", names(coef_focus)[1], sep = ""), "", estNames, fixed = T)
  #     if(cross_int == TRUE){
  #       estNames <- sub(paste(marginal_dist[[1]]$factor[m],"_rp", sep = ""),
  #                       marginal_dist[[1]]$factor[m], estNames)
  #     }
  #     table_AME_m <- c()
  #     # For each marginal distribution,
  #     for(z in 1:length(marginal_dist_u_list)){
  #       marginal_dist_u <- marginal_dist_u_list[[z]]
  #       # Find weights
  #       coef_prop <- c(1, as.numeric(as.character(marginal_dist_u[match(estNames, marginal_dist_u[, "level"]), "prop"]))[-1])
  #       # Compute AMEs
  #       coef_AME <- sum(coef_focus * coef_prop)
  #       # se_AME <- sqrt(coef_prop%*%vcov_focus%*%coef_prop)
  #       AME <- data.frame(matrix(NA, ncol = 0, nrow=1))
  #       AME$type <- marginal_type[z]
  #       AME$factor   <- marginal_dist[[z]][m,1]; AME$level <- marginal_dist[[z]][m,2]
  #       AME$estimate <- coef_AME;
  #       # AME$se <- se_AME
  #       table_AME_m <- rbind(table_AME_m, AME)
  #     }
  #     if(difference == TRUE){
  #       for(z in 2:length(marginal_dist_u_list)){
  #         marginal_dist_u <- marginal_dist_u_list[[z]]
  #         # Find weights
  #         coef_prop <- c(1, as.numeric(as.character(marginal_dist_u[match(estNames, marginal_dist_u[, "level"]), "prop"]))[-1])
  #         coef_prop0 <- c(1, as.numeric(as.character(marginal_dist_u_base[match(estNames, marginal_dist_u_base[, "level"]), "prop"]))[-1])
  #         # Compute AMEs
  #         coef_prop_d <- (coef_prop - coef_prop0)
  #         coef_AME_dif <- sum(coef_focus * coef_prop_d)
  #         # se_AME_dif <- sqrt(coef_prop_d%*%vcov_focus%*%coef_prop_d)
  #         AME_dif <- data.frame(matrix(NA, ncol = 0, nrow=1))
  #         AME_dif$type <- paste(marginal_type[z],"-",marginal_type[1],sep="")
  #         AME_dif$factor   <- marginal_dist[[z]][m,1]; AME_dif$level <- marginal_dist[[z]][m,2]
  #         AME_dif$estimate <- coef_AME_dif;
  #         # AME_dif$se <- se_AME_dif
  #         table_AME_m <- rbind(table_AME_m, AME_dif)
  #       }
  #     }
  #     table_AME <- rbind(table_AME, table_AME_m)
  #   }
  # }
  # colnames(table_AME) <- c("type", "factor", "level", "estimate")
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


# Functions don't explicitly model two-way interactions
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

    # if(is.null(pair_var) == TRUE){
    X1 <- model.matrix(formula, data=data1)[ ,-1]
    X2 <- model.matrix(formula, data=data2)[ ,-1]
    X <- cbind(1, X1 - X2)
    # }else{
    #   Xf <- model.frame(formula, data=data1)
    #   Xf_o <- as.numeric(attr(terms(Xf), "order") == 1) + 1
    #   Xf <- attr(terms(Xf), "factors")
    #   Xf_ind <- which(apply(t(t(matrix(Xf[pair_var, ], ncol = ncol(Xf))*Xf_o)), 2, sum)  == 2)
    #   X01 <- model.matrix(formula, data=data1)
    #   X02 <- model.matrix(formula, data=data2)
    #   pair_var_ind <- is.element(attr(X01, "assign"), Xf_ind)
    #   X0 <- X01 - X02
    #   X0[, pair_var_ind == TRUE] <- X01[, pair_var_ind == TRUE]
    #   X <- cbind(1, X0[,-1])
    # }

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

    # if(is.null(pair_var) == TRUE){
    #   X1 <- model.matrix(formula, data=data1)[ ,-1]
    #   X2 <- model.matrix(formula, data=data2)[ ,-1]
    #   X <- cbind(1, X1 - X2)
    # }else{
    #   Xf <- model.frame(formula, data=data1)
    #   Xf_o <- as.numeric(attr(terms(Xf), "order") == 1) + 1
    #   Xf <- attr(terms(Xf), "factors")
    #   Xf_ind <- which(apply(t(t(matrix(Xf[pair_var, ], ncol = ncol(Xf))*Xf_o)), 2, sum)  == 2)
    #   X01 <- model.matrix(formula, data=data1)
    #   X02 <- model.matrix(formula, data=data2)
    #   pair_var_ind <- is.element(attr(X01, "assign"), Xf_ind)
    #   X0 <- X01 - X02
    #   X0[, pair_var_ind == TRUE] <- X01[, pair_var_ind == TRUE]
    #   X <- cbind(1, X0[,-1])
    # }

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


# #######################
# Help functions
# #######################

# 1. Cluster standard errors
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

  cross_int <- outAME$cross_int
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

      if(cross_int == TRUE){
        estNames <- sub(paste(marginal_dist[[1]]$factor[m],"_rp", sep = ""), marginal_dist[[1]]$factor[m], estNames)
      }

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

  cross_int <- outAME$input$cross_int

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
  if(cross_int == TRUE){
    part_pos <- append(lapply(leveluse[!(names(leveluse) %in% factor_name)], length), length(leveluse[[factor_name]]))
    names(part_pos)[length(part_pos)] <- factor_name
    leveluse_all_b <- append(leveluse_all[!(names(leveluse_all) %in% factor_name)], leveluse_all[factor_name])
  }else{
    part_pos <- lapply(leveluse[!(names(leveluse) %in% factor_name)], length)
    leveluse_all_b <- leveluse_all[!(names(leveluse_all) %in% factor_name)]
  }
  end_pos   <- cumsum(part_pos) + 1
  start_pos <- c(2, end_pos[-length(end_pos)] + 1)

  estNames <- gsub(paste(effect_name, ":", sep = ""), "", colnames(coef_focus), fixed = T)
  estNames <- gsub(paste(":", effect_name, sep = ""), "", estNames, fixed = T)

  if(cross_int == TRUE){
    estNames <- sub(paste(factor_name,"_rp", sep = ""), factor_name, estNames)
  }


  # Estimate Conditional AMEs
  table_cAME <- table_cProp <- c()
  for(m in 1:length(part_pos)){
    ind_focus_mat <- matrix(1, ncol = (as.numeric(part_pos[m][1]) + 1), nrow = ncol(coef_focus))
    leveluse_focus <- leveluse_all_b[m]
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
                      main = c("Proportions", "Conditional AMEs"),
                      col,
                      pch, mar,
                      legend_pos = "topright",
                      plot_all = TRUE,
                      plot_difference = "no",
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
         yaxt="n", ylab = "", main = "Conditional AMCEs", xlim = p_x,
         xlab = "Estimated Effects", col = rev(p_col_full_effect))
    segments(rev(p_low_full), seq(1:length(p_coef_full)),
             rev(p_high_full), seq(1:length(p_coef_full)),
             col = rev(p_col_full_effect))
    Axis(side=2, at = seq(1:length(p_name_full)),
         labels=rev(p_name_full), las=1, font = 1, tick=F, cex.axis = cex)
    Axis(side=2, at = seq(1:length(p_name_f_full)),
         labels=rev(p_name_f_full), las=1, font = 2, tick=F, cex.axis = cex)
    abline(v=0, lty=2)
    if(is.character(legend_pos[1])==TRUE)  legend(legend_pos, marginal_effect, col= p_col_full_effect, pch = p_pch_full)
    if(is.character(legend_pos[1])==FALSE) legend(x=legend_pos[1], y=legend_pos[2],
                                                  marginal_effect, col= p_col_full_effect, pch = p_pch_full)

  }else if(plot_all == TRUE){
    par(mfrow=c(1,2), oma = c(1, mar,1,1), mar=c(4,2,4,1))
    ## Plot ----------
    plot(rev(p_cProp_full), seq(1:length(p_cProp_full)), pch=rev(p_pch_full_prop),
         yaxt="n", ylab = "", main = main[1], xlim = p_x_prop,
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
         yaxt="n", ylab = "", main = main[2], xlim = p_x,
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
  cross_int <- outAME$input$cross_int

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

  if(cross_int == TRUE){
    estNames <- sub(paste(factor_name,"_rp", sep = ""), factor_name, estNames)
  }

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
  if(cross_int == TRUE) marginal_factor <- c(marginal_factor, factor_name)
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
    AME_diff$type <- paste(marginal_diff[1], marginal_diff[2], sep = " - ")
    AME_diff$factor   <- marginal_factor[m]
    AME_diff$estimate <- coef_AME_diff;
    AME_diff$se <- se_AME_diff
    AME_diff$low.95ci  <- low_AME_diff
    AME_diff$high.95ci <- high_AME_diff
    table_AME_diff <- rbind(table_AME_diff, AME_diff)
  }
  return(table_AME_diff)
}

#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param pair whether we use the paired-conjoint design
#' @param pair_id id for paired-conjoint design. required when 'pair = TRUE'
#' @export

plot_decompose <- function(outAME, factor_name, level_name, type = "decompose",
                     marginal_diff, mar = 8){

  if(is.element(type, c("decompose", "diagnostics")) == FALSE){
    stop("type should be `decompose` or `diagnostics.`")
  }

  cross_int <- outAME$input$cross_int


  if(type == "decompose"){

    if(missing(marginal_diff)){
      if(length(outAME$input$marginal_type) == 1){
        stop("Cannot compare two marginal distributions because `AME_full_estimate`` only includes one marginal distribution.")
      }else{
      marginal_diff <- outAME$input$marginal_type[c(2, 1)]
      }
    }

    dec_tab <- decompose_AME(outAME = outAME, factor_name = factor_name,
                             level_name = level_name,
                             marginal_diff = marginal_diff)

    point <- rev(dec_tab[,3])
    low  <- rev(dec_tab[,5])
    high   <- rev(dec_tab[,6])
    fac_name_p <- rev(dec_tab[,2])

    xmin <- min(low); xmax <- max(high)

    par(mar = c(4, mar, 6, 4))
    plot(point, seq(1:nrow(dec_tab)), pch = 19, ylim = c(0.5, nrow(dec_tab)+0.5), yaxt = "n",
        xlim = c(xmin, xmax), ylab =  "", xlab = "Change in Popuation AMCE",
        main = paste("Decompose Change in Population AMCE:\n", unique(dec_tab$type), sep =""))
    segments(low, seq(1:nrow(dec_tab)), high, seq(1:nrow(dec_tab)), lwd = 2)
    abline(v = 0, lty = 2)
    Axis(side = 2, at = seq(1:nrow(dec_tab)), labels = fac_name_p, las = 1)
  }
}


#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param pair whether we use the paired-conjoint design
#' @param pair_id id for paired-conjoint design. required when 'pair = TRUE'
#' @export

plot_AME <- function(ame.out,
                     factor_name,
                     col, pch,
                     legend_pos = "topright",
                     plot_difference = "both",
                     main = "AMEs", xlim,
                     mar, cex = 1,
                     plot_type,
                     plot_name){

  if(all(is.element(factor_name, names(ame.out$AME))) == FALSE){
    stop(" 'factor_name' can only take a subset of factors estimated in 'AME.estimate'")
  }
  if(missing(col)){
    col <- rep("black",length(unique(ame.out$AME[[1]]$type)))
  }
  if(missing(pch)){
    pch <- rep(19, length(col))
  }
  if(missing(mar)){
    mar <- 12
  }
  if(ame.out$input$difference == FALSE & plot_difference == "only"){
    stop(" 'AME.estiamte' should take `difference = TRUE` to plot their estimates")
  }
  if(missing(plot_name)) plot_name <- plot_type
  type_difference <- ame.out$type_difference

  ## Correct all esimates ----------
  p_coef_full <- p_high_full <- p_low_full <- c()
  p_name_full <- p_name_f_full <- p_col_full <- p_pch_full <- c()
  for(g in 1:length(factor_name)){
    p_AME  <- ame.out$AME[[factor_name[g]]]
    p_AME <- p_AME[order(factor(p_AME$level, levels = unique(p_AME$level))),] # updated on 12/27 (Naoki)

    if(plot_difference == "both"){p_AME  <- p_AME}
    if(plot_difference == "no"){p_AME   <- p_AME[(p_AME$type %in% type_difference) == FALSE, ]}
    if(plot_difference == "only"){p_AME  <- p_AME[(p_AME$type %in% type_difference) == TRUE, ]}

    if(missing(plot_type) == FALSE){
      p_AME  <- p_AME[(p_AME$type %in% plot_type) == TRUE, ]
      p_AME  <- p_AME[order(factor(p_AME$level, levels = unique(p_AME$level)),
                            factor(p_AME$type, levels = plot_type)),]
    }

    p_coef <- c(NA, p_AME$estimate)
    p_se   <- c(NA, p_AME$se)
    p_high <- p_coef + 1.96*p_se
    p_low  <- p_coef - 1.96*p_se
    p_name_t <- paste(factor_name[g], " (", ame.out$baseline[factor_name[g]], "):   ", sep="")
    p_name_b   <- p_AME$level
    p_name_b[p_AME$type != p_AME$type[1]] <- ""
    p_name_f <- c(p_name_t, rep("", length(p_name_b)))
    p_name   <- c("", p_name_b)
    p_col  <- c(NA, rep(col, length(unique(p_AME$level))))
    p_pch  <- c(NA, rep(pch, length(unique(p_AME$level))))

    ## Store values
    p_coef_full <- c(p_coef_full, p_coef)
    p_high_full <- c(p_high_full, p_high)
    p_low_full  <- c(p_low_full, p_low)
    p_name_f_full <- c(p_name_f_full, p_name_f)
    p_name_full <- c(p_name_full, p_name)
    p_col_full  <- c(p_col_full, p_col)
    p_pch_full  <- c(p_pch_full, p_pch)
  }

  ## Plot Setup ----------
  p_type <- unique(p_AME$type)
  if(missing(plot_type) == FALSE) p_type <- plot_type
  if(missing(xlim)){
    p_x <- c(min(p_low_full, na.rm=TRUE), max(p_high_full, na.rm = TRUE))
  }else{
    p_x <- xlim
  }

  ## Plot ----------
  par(mar=c(4,mar,4,2))
  plot(rev(p_coef_full), seq(1:length(p_coef_full)), pch=rev(p_pch_full),
       yaxt="n", ylab = "", main = main, xlim = p_x,
       xlab = "Estimated Effects", col = rev(p_col_full))
  segments(rev(p_low_full), seq(1:length(p_coef_full)),
           rev(p_high_full), seq(1:length(p_coef_full)),
           col = rev(p_col_full))
  Axis(side=2, at = seq(1:length(p_name_full)), labels=rev(p_name_full),
       las=1, font = 1, tick=F, cex.axis = cex)
  Axis(side=2, at = seq(1:length(p_name_full)), labels=rev(p_name_f_full),
       las=1, font = 2, tick=F, cex.axis = cex)
  abline(v=0, lty=2)
  if(is.character(legend_pos[1])==TRUE) legend(legend_pos, plot_name, col= col, pch = pch)
  if(is.character(legend_pos[1])==FALSE) legend(x=legend_pos[1], y=legend_pos[2],
                                                plot_name, col= col, pch = pch)
}

#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param pair whether we use the paired-conjoint design
#' @param pair_id id for paired-conjoint design. required when 'pair = TRUE'
#' @export

plot_dist <- function(factor_name,
                      marginal_dist, marginal_type,
                      main = "Marginal Distibution",
                      col, pch, mar,
                      legend_pos = "topright",
                      cex){

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

  if(all(is.element(factor_name, unique(marginal_dist[[1]]$factor))) == FALSE){
    stop(" 'factor_name' can only take a subset of factors estimated in 'cAME.estimate'")
  }
  if(missing(col)){
    col <- rep("black",length(unique(marginal_dist[[1]]$levels)))
  }
  if(missing(pch)){
    pch <- rep(19,length(unique(marginal_dist[[1]]$levels)))
  }
  if(missing(mar)){
    mar <- 6
  }

  ## Collect all Distributions ----------
  p_mar_coef <- list()
  for(z in 1: length(marginal_dist)){
    p_coef_full <- c()
    p_name_full <- p_name_f_full <- p_col_full <- p_pch_full <- c()
    for(g in 1:length(factor_name)){
      p_mar  <- marginal_dist[[z]][marginal_dist[[z]]$factor == factor_name[g], ]

      p_coef <- c(NA, p_mar$prop)
      p_name_t <- paste(factor_name[g], ":   ", sep="")
      p_name_b   <- as.character(p_mar$levels)
      p_name_b[p_mar$type != p_mar$type[1]] <- ""
      p_name_f <- c(p_name_t, rep("", length(p_name_b)))
      p_name   <- c("", p_name_b)
      p_col  <- c(NA, rep(col[z], length(unique(p_mar$levels))))
      p_pch  <- c(NA, rep(pch[z], length(unique(p_mar$levels))))

      ## Store values
      p_coef_full  <- c(p_coef_full, p_coef)
      p_name_f_full <- c(p_name_f_full, p_name_f)
      p_name_full <- c(p_name_full, p_name)
      p_col_full  <- c(p_col_full, p_col)
      p_pch_full  <- c(p_pch_full, p_pch)
    }
    p_mar_coef$p_coef_full[[z]] <- p_coef_full
    p_mar_coef$p_name_full[[z]] <- p_name_full
    p_mar_coef$p_name_f_full[[z]] <- p_name_f_full
    p_mar_coef$p_col_full[[z]] <- p_col_full
    p_mar_coef$p_pch_full[[z]] <- p_pch_full
  }


  ## Plot Setup for Proportion----------
  p_type_prop <- marginal_type
  p_x_prop <- c(0, max(unlist(p_mar_coef$p_coef_full), na.rm = TRUE))

  ## Plot ----------

  par(oma = c(1, mar,1,1), mar=c(4,2,4,1))
  ## Plot ----------
  plot(rev(p_mar_coef$p_coef_full[[1]]), seq(1:length(p_mar_coef$p_coef_full[[1]])),
       pch=pch[1],
       yaxt="n", ylab = "", main = main, xlim = p_x_prop, type = "o",
       ylim = c(1, length(p_name_f_full)),
       xlab = "Proportions", col = col[1])
  if(length(marginal_dist) > 1){
    for(z in 2:length(marginal_dist)){
      points(rev(p_mar_coef$p_coef_full[[z]]),
             seq(1:length(p_mar_coef$p_coef_full[[z]])) + 0.05*z,
             pch=pch[z], type = "o", col = col[z])
    }
  }
  Axis(side=2, at = 0.25 + seq(1:length(p_mar_coef$p_name_full[[1]])),
       labels = rev(p_mar_coef$p_name_full[[1]]), las=1, font = 1, tick=F,
       cex.axis = cex)
  Axis(side=2, at = 0.25 + seq(1:length(p_mar_coef$p_name_f_full[[1]])),
       labels = rev(p_mar_coef$p_name_f_full[[1]]), las=1, font = 2, tick=F,
       cex.axis = cex)
  abline(v=0, lty=2)
  if(is.character(legend_pos[1])==TRUE)  legend(legend_pos, p_type_prop, col= col, pch = pch)
  if(is.character(legend_pos[1])==FALSE) legend(x=legend_pos[1], y=legend_pos[2],
                                                p_type_prop, col= col, pch = pch)

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
    # factor_u <-  rep(names(all_levels), unlist(lapply(all_levels, length)))
    # levels_u <- c(unlist(all_levels))
    target   <- lapply(dataX, function(x) prop.table(c(table(x))))

    # target <- data.frame(matrix(NA, ncol = 0, nrow = length(factor_u)))
    # target$factor   <- factor_u
    # target$levels <- levels_u
    # target$prop   <- prop_u
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


#' Create Profile Distributions from data.frame
#' @param formula formula
#' @param data data.frame
#' @param marginal Whether we create marginal distributions (TRUE) or the two-dimensional joint distributions
#' @export

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

  # levels_u <- names(unlist(prop_u2))
  # factor_u <-  rep(fac_name2, unlist(lapply(prop_u2, length)))
  #
  # target <- data.frame(matrix(NA, ncol = 0, nrow = length(factor_u)))
  # target$factor   <- factor_u
  # target$levels <- levels_u
  # target$prop   <- unlist(prop_u2)

  attributes(target)$dist_type <- "marginal"

  return(target)
}

Marginal2Joint <- function(marginal_dist){

  # fac_name   <- unique(marginal_dist$factor)
  #
  # all_levels <- tapply(marginal_dist$prop, marginal_dist$factor, function(x) x)
  # all_levels <- all_levels[match(fac_name, names(all_levels))]
  #
  # level_names <- tapply(marginal_dist$levels, marginal_dist$factor, function(x) x)
  # level_names <- level_names[match(fac_name, names(level_names))]

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
    # if(all(colnames(check_dist) == c("factor", "levels", "prop")) == FALSE){
    #   stop(" When target_type = marginal, 'colnames' of 'target_dist' should be c('factor', 'levels', 'prop')")
    # }

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


