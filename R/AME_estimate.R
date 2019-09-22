## (Internal Function): Estimating pAMCE without regularization

AME_estimate <- function(formula,
                         data,
                         pair = FALSE, pair_id = NULL, cross_int = TRUE,
                         cluster,
                         marginal_dist,
                         marginal_type,
                         joint_dist = NULL,
                         boot,
                         difference = FALSE, formula_three_c){

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
  # Create two-way interaction formula (allow for three-ways) ----------
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

  # Baseline ----------
  baseline <- lapply(model.frame(formula, data=data)[,-1],
                     FUN = function(x) levels(x)[1])

  # pair_id and cluster
  if(missing(pair_id) == TRUE) pair_id <- NULL
  if(missing(cluster) == TRUE) cluster <- seq(1:nrow(data))
  data$pair_id <- pair_id
  data$cluster <- cluster
  # Differencing ----------
  if(pair==TRUE){
    data0 <- data[order(pair_id),]
    side <- rep(c(1,0), times=nrow(data0)/2)
    data1 <- data0[side==1,]
    data2 <- data0[side==0,]
    cluster_original <- cluster
    cluster <- cluster[side==1]

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
      X_cross <- data_cross[, c((sing + 1): ncol(data_cross))]

      # modify X and ind_b
      X <- cbind(X, X_cross)
    }

    y <- model.frame(formula_full,data=data1)[ ,1]
    # base_name <- c("(Intercept)", colnames(X1))
  }else{
    cluster_original <- cluster
    X <- model.matrix(formula_full, data=data)
    y <- model.frame(formula_full, data=data)[,1]
    # base_name <- colnames(X)
    side <- NULL
  }

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

  # Fit the model ----------
  main_lm <- lm(y ~ X - 1)
  coefInt <- coef(main_lm)
  coefInt <- coefInt[is.na(coefInt) == FALSE]
  base_name <- sub("X", "", names(coefInt))
  names(coefInt) <- base_name
  vcovInt <- cluster_se_glm(main_lm, cluster = as.factor(cluster))
  colnames(vcovInt) <- rownames(vcovInt) <- base_name

  # Transform marginal_dist (for internal simplisity) ----------
  marginal_dist_u_list <- list()
  for(z in 1:length(marginal_dist)){
    marginal_dist_u_list[[z]] <- data.frame(matrix(NA, ncol=0, nrow=nrow(marginal_dist[[z]])))
    marginal_dist_u_list[[z]]$level <- paste(marginal_dist[[z]][,1], marginal_dist[[z]][,2],sep="")
    marginal_dist_u_list[[z]]$prop  <- marginal_dist[[z]][,3]
  }
  marginal_dist_u_base <- marginal_dist_u_list[[1]]

  # Incorporate Joint Distribution
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
    joint_dist_u_base <- joint_dist_u_list[[1]]
  }

  # Estimate AMEs ----------
  ## Estimate AMEs from two-ways
  if(is.null(formula_three_c) == TRUE){
    table_AME <- coefIntAME(coefInt = coefInt, vcovInt = vcovInt, SE = TRUE,
                            marginal_dist = marginal_dist, marginal_dist_u_list = marginal_dist_u_list,
                            marginal_dist_u_base = marginal_dist_u_base, marginal_type = marginal_type,
                            difference = difference, cross_int = cross_int)
  }else if(is.null(formula_three_c) == FALSE){
    table_AME <- coefIntAME(coefInt = coefInt, vcovInt = vcovInt, SE = TRUE,
                            marginal_dist = marginal_dist, marginal_dist_u_list = marginal_dist_u_list,
                            marginal_dist_u_base = marginal_dist_u_base, marginal_type = marginal_type,
                            difference = difference, cross_int = cross_int, three_way = TRUE,
                            joint_dist_u_list = joint_dist_u_list)
  }

  # table_AME <- c()
  # for(m in 1:nrow(marginal_dist_u_base)){
  #   coef_focus <- coefInt[grep(marginal_dist_u_base$level[m], names(coefInt), fixed = T)]
  #   vcov_focus <- vcovInt[grep(marginal_dist_u_base$level[m], rownames(vcovInt), fixed = T),
  #                         grep(marginal_dist_u_base$level[m], colnames(vcovInt), fixed = T)]
  #   if(length(coef_focus) > 0){
  #     estNames <- gsub(paste(marginal_dist_u_base$level[m], ":", sep = ""), "", names(coef_focus), fixed = T)
  #     estNames <- gsub(paste(":", names(coef_focus)[1], sep = ""), "", estNames, fixed = T)
  #
  #     if(cross_int == TRUE){
  #       estNames <- sub(paste(marginal_dist[[1]]$factor[m],"_rp", sep = ""),
  #                       marginal_dist[[1]]$factor[m], estNames)
  #     }
  #
  #     table_AME_m <- c()
  #     # For each marginal distribution,
  #     for(z in 1:length(marginal_dist)){
  #       marginal_dist_u <- marginal_dist_u_list[[z]]
  #       # Find weights
  #       coef_prop <- c(1, as.numeric(as.character(marginal_dist_u[match(estNames, marginal_dist_u[, "level"]), "prop"]))[-1])
  #       # Compute AMEs
  #       coef_AME <- sum(coef_focus * coef_prop)
  #       se_AME <- sqrt(coef_prop%*%vcov_focus%*%coef_prop)
  #       AME <- data.frame(matrix(NA, ncol = 0, nrow=1))
  #       AME$type <- marginal_type[z]
  #       AME$factor   <- marginal_dist[[z]][m,1]; AME$level <- marginal_dist[[z]][m,2]
  #       AME$estimate <- coef_AME; AME$se <- se_AME
  #       table_AME_m <- rbind(table_AME_m, AME)
  #     }
  #     if(difference == TRUE){
  #       for(z in 2:length(marginal_dist)){
  #         marginal_dist_u <- marginal_dist_u_list[[z]]
  #         # Find weights
  #         coef_prop <- c(1, as.numeric(as.character(marginal_dist_u[match(estNames, marginal_dist_u[, "level"]), "prop"]))[-1])
  #         coef_prop0 <- c(1, as.numeric(as.character(marginal_dist_u_base[match(estNames, marginal_dist_u_base[, "level"]), "prop"]))[-1])
  #         # Compute AMEs
  #         coef_prop_d <- (coef_prop - coef_prop0)
  #         coef_AME_dif <- sum(coef_focus * coef_prop_d)
  #         se_AME_dif <- sqrt(coef_prop_d%*%vcov_focus%*%coef_prop_d)
  #         AME_dif <- data.frame(matrix(NA, ncol = 0, nrow=1))
  #         AME_dif$type <- paste(marginal_type[z],"-",marginal_type[1],sep="")
  #         AME_dif$factor   <- marginal_dist[[z]][m,1]; AME_dif$level <- marginal_dist[[z]][m,2]
  #         AME_dif$estimate <- coef_AME_dif; AME_dif$se <- se_AME_dif
  #         table_AME_m <- rbind(table_AME_m, AME_dif)
  #       }
  #     }
  #     table_AME <- rbind(table_AME, table_AME_m)
  #   }
  # }
  # colnames(table_AME) <- c("type", "factor", "level", "estimate", "se")

  # Combine STD estimate
  table_STD <- AME.fit.STD.se.sep(formula = formula,
                                  data = data,
                                  pair = pair,
                                  marginal_dist = marginal_dist,
                                  marginal_dist_u_base = marginal_dist_u_base)

  uniq_fac <- unique(table_AME$factor)
  table_AME_full <- c()
  for(i in 1:length(uniq_fac)){
    table_AME_main <- table_AME[table_AME$factor == uniq_fac[i], ]
    uniq_level <- unique(table_AME_main$level)
    for(j in 1:length(uniq_level)){
      table_AME_main_m <- table_AME_main[table_AME_main$level == uniq_level[j], ]
      STD_base <- table_STD[table_STD$fac == uniq_fac[i] & table_STD$level == uniq_level[j], ]
      table_AME_main_m <- rbind(STD_base, table_AME_main_m)

      if(difference == TRUE){
        table_AME_add <- data.frame(matrix(NA, ncol = 0, nrow = length(marginal_type)))
        dif_est  <- table_AME_main_m[2:(1+length(marginal_type)), "estimate"] - STD_base$estimate
        dif_se   <- sqrt((table_AME_main_m[2:(1+length(marginal_type)), "se"])^2 + STD_base$se^2)
        table_AME_add$type <- paste(marginal_type,"-sample AMCE",sep="")
        table_AME_add$factor <- rep(uniq_fac[i], length(marginal_type))
        table_AME_add$level <- rep(uniq_level[j], length(marginal_type))
        table_AME_add$estimate <- dif_est
        table_AME_add$se <- dif_se
        table_AME_main_m <- rbind(table_AME_main_m, table_AME_add)
      }
      table_AME_full <- rbind(table_AME_full, table_AME_main_m)
    }
  }

  ##  bootstrap coefficients
  boot_coef <- rmvnorm(n = boot, mean = coefInt, sigma = vcovInt)



  ## For Each Factor
  AME <- list()
  for(g in 1:length(unique(table_AME_full$factor))){
    AME[[g]] <- table_AME_full[table_AME_full$factor == unique(table_AME_full$factor)[g], ]
    AME[[g]]$low.95ci <- AME[[g]]$estimate - 1.96*AME[[g]]$se
    AME[[g]]$high.95ci <- AME[[g]]$estimate + 1.96*AME[[g]]$se
  }
  names(AME) <- unique(table_AME_full$factor)
  type_all   <- unique(table_AME_full$type)
  type_difference   <- setdiff(unique(table_AME_full$type), marginal_type)

  input  <- list("formula" = formula, "data" = data,
                 "pair" = pair, "pair_id" = pair_id, "cross_int" = cross_int,
                 "cluster" = cluster_original, "marginal_dist" = marginal_dist,
                 "marginal_type" = marginal_type, "difference" = difference)

  output <- list("AMCE" = AME, "baseline" = baseline, "coef" = coefInt,
                 "type_all" = type_all, "type_difference" = type_difference,
                 "boot_coef" = boot_coef,
                 "input" = input)
  return(output)
}

