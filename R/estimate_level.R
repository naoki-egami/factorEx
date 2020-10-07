#' Estimating the population AMCE using a model-based approach
#' @param object The output from \code{model_pAMCE} with \code{reg = FALSE}.
#' @export

estimate_level <- function(object){

  marginal_dist <- object$input$marginal_dist
  target_dist <- object$input$target_dist
  target_type <- object$input$target_type
  formula <- object$input$formula
  data <- object$input$data
  marginal_dist_u_base <- cbind(marginal_dist[[1]][, 1], object$input$marginal_dist_u_base)
  colnames(marginal_dist_u_base)[1] <- "factor"

  marginal_dist_u_list <- object$input$marginal_dist_u_list

  coefInt <- object$coef
  vcovInt <- object$vcov

  target_name <- c("sample", names(target_dist))
  Sample_Mar <- createDist(formula = formula, target_data = data,
                           exp_data = data, type  = "marginal")
  target_dist  <- c(list(Sample_Mar), target_dist)
  names(target_dist)[1] <- "sample"  # (marginal distributions used for randomization)
  target_type <- c("marginal", target_type)

  # ########################################
  # Check each target profile distribution #
  # ########################################
  for(i in 1:length(target_dist)){
    target_dist[[i]] <- checkDist(target_dist[[i]], type = target_type[i], formula = formula, data = data)
  }

  # Prepare Joint Distribution
  joint_dist <- list()
  for(i in 1:length(target_dist)){
    if(target_type[i]  == "marginal"){joint_dist[[i]] <- Marginal2Joint(target_dist[[i]])}
    if(target_type[i]  == "target_data"){joint_dist[[i]] <- createDist(formula = formula,
                                                                       target_data = target_dist[[i]],
                                                                       exp_data = data, type  = "joint")}
    joint_dist[[i]] <- checkDist(joint_dist[[i]], type = "joint", formula = formula, data = data)
  }

  joint_dist_u_list <- list()
  for(i in 1:length(joint_dist)){
    temp <- do.call("rbind", joint_dist[[i]])
    joint_dist_u_list[[i]] <- data.frame(matrix(NA, ncol=0, nrow=nrow(temp)))
    joint_dist_u_list[[i]]$level_1 <- paste(temp[,1], temp[,3],sep="")
    joint_dist_u_list[[i]]$level_2 <- paste(temp[,2], temp[,4],sep="")
    joint_dist_u_list[[i]]$prop  <- temp[,5]
  }

  # Estimate AMEs ----------
  table_AME <- c()
  for(m in 1:nrow(marginal_dist_u_base)){

    main_fac <- which(marginal_dist_u_base$factor == marginal_dist_u_base$factor[m])
    coef_focus0 <- c()
    for(q in 1:length(main_fac)){
      coef_focus0 <- c(coefInt[grep(marginal_dist_u_base$level[main_fac[q]], names(coefInt), fixed = T)],
                       coef_focus0)
    }
    coef_other <- coefInt[is.element(names(coefInt), names(coef_focus0)) == FALSE]

    coef_focus1 <- coefInt[grep(marginal_dist_u_base$level[m], names(coefInt), fixed = T)]
    coef_focus <- c(coef_focus1, coef_other)
    match_ord <- match(names(coef_focus), names(coefInt))
    coef_focus <- coef_focus[order(match_ord)]

    vcov_focus <- vcovInt[names(coef_focus), names(coef_focus)]

    estNames <- gsub(paste(marginal_dist_u_base$level[m], ":", sep = ""), "", names(coef_focus), fixed = T)
    estNames <- gsub(paste(":", marginal_dist_u_base$level[m], sep = ""), "", estNames, fixed = T)

    # if(cross_int == TRUE){
    #   estNames <- sub(paste(marginal_dist[[1]]$factor[m],"_rp", sep = ""),
    #                   marginal_dist[[1]]$factor[m], estNames)
    # }
    # allow for three-ways
    #if(three_way == TRUE){
    estL <- strsplit(estNames, ":")
    estL <- do.call("rbind", lapply(estL, function(x) if(length(x)==1){c(x[1], NA)}else{x}))
    #}

    table_AME_m <- c()
    # For each marginal distribution,
    for(z in 1:length(marginal_dist_u_list)){
      marginal_dist_u <- marginal_dist_u_list[[z]]
      joint_dist_u <- joint_dist_u_list[[z]]
      # Find weights
      #if(three_way == FALSE){
      #  coef_prop <- c(1, as.numeric(as.character(marginal_dist_u[match(estNames, marginal_dist_u[, "level"]), "prop"]))[-1])
      #}else if(three_way == TRUE){
      # allow for three-ways
      ind_two  <- which(apply(estL, 1, function(x) any(is.na(x))))[-1]  # Removes (Intercept)
      prop_two <- marginal_dist_u[match(estL[ind_two,1], marginal_dist_u[, "level"]), "prop"]
      ind_three  <- which(apply(estL, 1, function(x) all(is.na(x) == FALSE)))
      prop_three <- joint_dist_u[row.match(as.data.frame(matrix(estL[ind_three, ], ncol = 2)),
                                           as.data.frame(joint_dist_u[, 1:2])), "prop"]
      coef_prop  <- rep(NA, nrow(estL))
      coef_prop[1] <- 1  # Intercept
      coef_prop[ind_two]   <-  prop_two
      coef_prop[ind_three] <- prop_three
      coef_prop[estNames == marginal_dist_u_base$level[m]] <- 1 # Change for the main level
      #}

      # Compute AMEs
      coef_AME <- sum(coef_focus * coef_prop)
      se_AME <- sqrt(coef_prop%*%vcov_focus%*%coef_prop)
      AME <- data.frame(matrix(NA, ncol = 0, nrow=1))
      AME$type <- target_name[z]
      AME$factor   <- marginal_dist[[z]][m,1]; AME$level <- marginal_dist[[z]][m,2]
      AME$estimate <- coef_AME;
      AME$se <- se_AME
      table_AME_m <- rbind(table_AME_m, AME)
    }
    table_AME <- rbind(table_AME, table_AME_m)

  }
  colnames(table_AME) <- c("type", "factor", "level", "estimate", "se")
  return(table_AME)
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
