#' Computing weights for estimating the population AMCE using a design-based approach
#' @param formula Formula
#' @param factor_name A factor for which the function computes weights
#' @param data Data
#' @param pair Whether we use a paired-choice conjoint design
#' @param pair_id Unique identifiers for pairs in the paired-choice conjoint design  (optional)
#' @param cross_int Include interactions across profiles. Default is FALSE
#' @param target_dist Target profile distributions to be used. This argument should be `list`
#' @param target_type Types of target profile distributions. `marginal`, `partial_joint` or `target_data`. See Examples for details
#' @param partial_joint_name Names of factors representing partial joint distributions. See Examples for details.
#' @export

weights_pAMCE <- function(formula, factor_name, data, pair, pair_id, cross_int,
                          target_dist,  target_type,
                          partial_joint_name){

  # Housekeeping
  if(all(target_type %in% c("marginal", "partial_joint", "target_data")) == FALSE){
    stop(" 'target_type' should be 'marginal', 'partial_joint', or 'target_data' ")
  }
  if(missing(pair_id) == TRUE) pair_id <- NULL
  if(pair==TRUE & is.null(pair_id)==TRUE){
    stop("When 'pair=TRUE', specify 'pair_id'.")
  }
  if(pair==TRUE & all(table(pair_id)==2)==FALSE){
    stop("When 'pair=TRUE', each of 'pair_id' should have two observations.")
  }
  if(is.null(pair_id) == FALSE){ pair <- TRUE }
  if(pair == FALSE){ cross_int <- FALSE }



  # Case 1: both are marginal
  if(target_type  == "marginal"){

    if(!inherits(target_dist, "list")){
      stop(" if 'target_type = marginal', class(target_dist) should be 'list' ")
    }

    randomize_dist <- createDist(formula, target_data = data, exp_data  = data, type = "marginal")

    out <- weights_pAMCE_mar(formula = formula, factor_name = factor_name,
                             data = data, pair = pair, pair_id =  pair_id, cross_int = cross_int,
                             randomize_dist = randomize_dist,
                             target_dist =  target_dist)
  }
  # Case 2: both are data
  if(target_type  == "target_data"){

    if(!inherits(target_dist, "data.frame")){
      stop(" if 'target_type = target_data', class(target_dist) should be 'data.frame' ")
    }

    factor_use <- all.vars(formula)[-1]
    randomize_dist <- model.frame(formula, data = data)[, -1]

    if(setequal(factor_use, colnames(target_dist)) == FALSE){
      stop(" colnames(target_dist) should be the same as names of factors in 'formula' ")
    }

    out <- weights_pAMCE_data(formula = formula, factor_name = factor_name,
                              data = data, pair = pair, pair_id =  pair_id, cross_int = cross_int,
                              randomize_dist = randomize_dist,
                              target_dist =  target_dist)
  }
  # Case 3: both are joint partial
  if(target_type  == "partial_joint"){

    if(missing(partial_joint_name) == TRUE){
      stop(" if 'target_type = partial_joint', please specify 'partial_joint_name' ")
    }

    if(!inherits(target_dist, "list")){
      stop(" if 'target_type = partial_joint', class(target_dist) should be 'list' ")
    }

    randomize_dist <- list()
    for(i in 1:length(partial_joint_name)){
      randomize_dist[[i]]  <- table(data[,  partial_joint_name[[i]]])/nrow(data)
    }

    out <- weights_pAMCE_partial(formula = formula, factor_name = factor_name,
                                 data = data, pair = pair, pair_id =  pair_id, cross_int = cross_int,
                                 randomize_dist = randomize_dist,
                                 target_dist =  target_dist, partial_joint_name = partial_joint_name)
  }

  return(out)


}

weights_pAMCE_partial <- function(formula, factor_name, data, pair, pair_id, cross_int =FALSE,
                                  randomize_dist, target_dist,
                                  partial_joint_name){

  data$internal_id <- seq(1:nrow(data))
  if(pair == TRUE){data$pair_id <- pair_id}

  # HouseKeeping (target_data Case)
  factor_use <- all.vars(formula)[-1]
  level_use  <- lapply(data[, factor_use], levels)

  if(length(unique(unlist(partial_joint_name))) !=  length(unlist(partial_joint_name))){
    stop(" factors should not be overlapped across partial joint distributions. unique(unlist(partial_joint_name)) should be equal to unlist(partial_joint_name). ")
  }

  if(setequal(factor_use, unlist(partial_joint_name))  == FALSE){
    stop(" colnames(randomize_dist) should be the same as names of factors in 'formula' ")
  }

  # check factor names
  partial_joint_name_use <-  c()
  for(i in 1:length(partial_joint_name)){
    if(length(partial_joint_name[[i]])  >  1){
      partial_joint_name_use[i] <- paste(partial_joint_name[[i]],  collapse =  ":")
    }else{
      partial_joint_name_use[i] <- partial_joint_name[[i]]
    }
  }
  names(randomize_dist) <- names(target_dist)  <- partial_joint_name_use
  randomize_dist <- randomize_dist[match(partial_joint_name_use, names(randomize_dist))]
  target_dist <- target_dist[match(partial_joint_name_use, names(target_dist))]

  # check level names
  for(z in  1:length(randomize_dist)){
    for(j in 1:length(dim(randomize_dist[[z]]))){
      if(length(dim(randomize_dist[[z]])) == 1){
        if(all(names(randomize_dist[[z]]) == names(target_dist[[z]]))  == FALSE){
          stop(" level names of 'target_dist' should be the same order as the order of level names in 'data' ")
        }
      }else{
        if(all(dimnames(randomize_dist[[z]])[[j]] == dimnames(target_dist[[z]])[[j]])  == FALSE){
          stop(" level names of 'target_dist' should be the same order as the order of level names in 'data' ")
        }
      }
    }
  }


  # Case 3: When both are partial joints
  randomize_prob  <- target_prob <- target_prob_s <- target_dist_s <- list()
  partial_joint_name_s <- partial_joint_name
  for(i in 1:length(partial_joint_name)){

    randomize_prob[[i]] <- c(randomize_dist[[i]])
    target_prob[[i]] <- c(target_dist[[i]])

    if((factor_name %in%  partial_joint_name[[i]]) ==TRUE){
      if(length(partial_joint_name[[i]]) == 1){
        target_prob_s[[i]] <- NA
        partial_joint_name_s[[i]] <- NA
        target_dist_s[[i]] <- NA
      }else{
        margin <- setdiff(seq(1:length(partial_joint_name[[i]])), which(partial_joint_name[[i]] == factor_name))
        target_prob_s[[i]]  <- c(apply(target_dist[[i]], margin, sum))
        partial_joint_name_s[[i]] <- setdiff(partial_joint_name[[i]],  factor_name)
        target_dist_s[[i]] <- as.table(apply(target_dist[[i]], margin, sum))
      }
    }else{
      target_prob_s[[i]] <- c(target_dist[[i]])
      target_dist_s[[i]] <-  as.table(target_dist[[i]])
    }
  }
  target_dist_s[is.na(target_dist_s)] <- NULL
  target_prob_s[is.na(target_prob_s)] <- NULL
  partial_joint_name_s[is.na(partial_joint_name_s)] <- NULL

  # keep track of names (long)
  name_match  <- list()
  for(i in 1:length(target_dist)){
    if(length(dim(randomize_dist[[i]])) == 1){
      name_match[[i]] <- paste(partial_joint_name[[i]], names(randomize_dist[[i]]), sep  = "")
    }else{
      name_l <- list()
      for(j in 1:length(dimnames(target_dist[[i]]))){
        name_l[[j]] <- as.character(paste(partial_joint_name[[i]][j], dimnames(target_dist[[i]])[[j]], sep  = ""))
      }
      exp_b <- expand.grid(name_l); for(m in 1:ncol(exp_b)){exp_b[,m] <- as.character(exp_b[,m])}
      exp_name <- c()
      for(z in 1:nrow(exp_b)){
        exp_name[z] <- paste(exp_b[z, ], collapse = ":")
      }
      name_match[[i]] <- exp_name
    }
  }
  # keep track of names (short)
  name_match_s  <- list()
  for(i in 1:length(target_dist_s)){
    if(length(dim(target_dist_s[[i]])) == 1){
      name_match_s[[i]] <- paste(partial_joint_name_s[[i]], names(target_dist_s[[i]]), sep  = "")
    }else{
      name_l_s <- list()
      for(j in 1:length(dimnames(target_dist_s[[i]]))){
        name_l_s[[j]] <- as.character(paste(partial_joint_name_s[[i]][j],
                                            dimnames(target_dist_s[[i]])[[j]], sep  = ""))
      }
      exp_b <- expand.grid(name_l_s); for(m in 1:ncol(exp_b)){exp_b[,m] <- as.character(exp_b[,m])}
      exp_name <- c()
      for(z in 1:nrow(exp_b)){
        exp_name[z] <- paste(exp_b[z, ], collapse = ":")
      }
      name_match_s[[i]] <- exp_name
    }
  }
  partial_joint_name_use_s <-  c()
  for(i in 1:length(partial_joint_name_s)){
    if(length(partial_joint_name_s[[i]])  >  1){
      partial_joint_name_use_s[i] <- paste(partial_joint_name_s[[i]],  collapse =  ":")
    }else{
      partial_joint_name_use_s[i] <- partial_joint_name_s[[i]]
    }
  }
  names(target_dist_s)  <- partial_joint_name_use_s

  # how to specify formula
  formula_x <- as.formula(paste("~", paste(partial_joint_name_use, collapse = "+"), sep = ""))
  formula_x_s <- as.formula(paste("~", paste(partial_joint_name_use_s, collapse = "+"), sep = ""))

  # the case of cross_int ==  FALSE
  if(cross_int == FALSE){
    new_id <- data$internal_id

    # Target Prob
    data_x_s <-  model.matrixBayes(formula_x_s, data = data, keep.order = TRUE, drop.baseline=FALSE)
    data_x_s <- data_x_s[, match(unlist(name_match_s), colnames(data_x_s))]
    data_x_s[data_x_s ==  0] <- NA

    data_p_s <- t(t(data_x_s)*unlist(target_prob_s))
    target_prob_f <- apply(data_p_s, 1, function(x) prod(x,  na.rm = TRUE))

    # Randomizing Prob
    data_x <-  model.matrixBayes(formula_x, data = data, keep.order = TRUE, drop.baseline=FALSE)
    data_x <- data_x[, match(unlist(name_match), colnames(data_x))]
    data_x[data_x ==  0] <- NA

    data_p_r <- t(t(data_x)*unlist(randomize_prob))
    randomize_prob_f <- apply(data_p_r, 1, function(x) prod(x,  na.rm = TRUE))

  }
  if(cross_int == TRUE){
    data0 <- data[order(data$pair_id),]
    side <- rep(c(1,0), times=nrow(data0)/2)
    data1 <- data0[side==1,]
    data2 <- data0[side==0,]
    new_id <- c(data1$internal_id, data2$internal_id)

    data1_x <-  model.matrixBayes(formula_x, data = data1)
    data2_x <-  model.matrixBayes(formula_x, data = data2)
    data1_x_s <-  model.matrixBayes(formula_x_s, data = data1)
    data2_x_s <-  model.matrixBayes(formula_x_s, data = data2)

    # combine
    data_x_s <- rbind(data1_x_s, data2_x_s)
    data_x <- rbind(data1_x, data2_x)
    data_x_rp <- rbind(data2_x, data1_x)
    data_x_s[data_x_s ==  0] <- NA
    data_x[data_x ==  0] <- NA
    data_x_rp[data_x_rp ==  0] <- NA

    # Prob in Target Distribution
    data_p_s <- t(t(data_x_s)*unlist(target_prob_s))
    target_prob_f1 <- apply(data_p_s, 1, function(x) prod(x,  na.rm = TRUE))

    data_p_s2 <- t(t(data_x_rp)*unlist(target_prob))
    target_prob_f2 <- apply(data_p_s2, 1, function(x) prod(x,  na.rm = TRUE))
    target_prob_f <- target_prob_f1*target_prob_f2

    # Prob in Randomizing Distribution
    data_p_r1 <- t(t(data_x)*unlist(randomize_prob))
    randomize_prob_f1 <- apply(data_p_r1, 1, function(x) prod(x,  na.rm = TRUE))

    data_p_r2 <- t(t(data_x_rp)*unlist(randomize_prob))
    randomize_prob_f2 <- apply(data_p_r2, 1, function(x) prod(x,  na.rm = TRUE))

    randomize_prob_f <- randomize_prob_f1*randomize_prob_f2
  }

  design_w <- target_prob_f/randomize_prob_f
  design_w[design_w > 1000] <- 0
  design_w[design_w > 10] <- 10
  design_weight <- design_w[match(data$internal_id, new_id)]

  newdata <- data
  newdata$design_weight <- design_weight


  output <- list("design_weight" = design_weight, "newdata" = newdata)
  return(output)
}





weights_pAMCE_data <- function(formula, factor_name, data, pair, pair_id, cross_int,
                               randomize_dist, target_dist){

  data$internal_id <- seq(1:nrow(data))
  if(pair == TRUE){data$pair_id <- pair_id}

  # HouseKeeping (target_data Case)
  factor_use <- all.vars(formula)[-1]
  level_use  <- lapply(data[, factor_use], levels)

  if(setequal(factor_use, colnames(randomize_dist))  == FALSE){
    stop(" colnames(randomize_dist) should be the same as names of factors in 'formula' ")
  }
  if(setequal(factor_use, colnames(target_dist))  == FALSE){
    stop(" colnames(target_dist) should be the same as names of factors in 'formula' ")
  }

  randomize_dist <- randomize_dist[, match(factor_use, colnames(randomize_dist))]
  target_dist <- target_dist[, match(factor_use, colnames(target_dist))]
  randomize_level_use <- lapply(randomize_dist, levels)
  target_level_use <- lapply(target_dist, levels)
  for(z in  1:length(randomize_dist)){
    if(setequal(level_use[[z]], randomize_level_use[[z]])  == FALSE){
      stop(" level names of 'randomize_dist' should be the same as level names in 'data' ")
    }
    if(setequal(level_use[[z]], randomize_level_use[[z]])  == FALSE){
      stop(" level names of 'target_dist' should be the same as level names in 'data' ")
    }
  }


  # Case 2: When both are data
  formula_x <- as.formula(paste(as.character(formula)[1], as.character(formula)[3], sep =  ""))
  formula_x_s <- as.formula(paste(as.character(formula)[1],
                                  paste(setdiff(all.vars(formula)[-1], factor_name), collapse = "+"),
                                  sep =  ""))

  # the case of cross_int ==  FALSE
  if(cross_int == FALSE){
    new_id <- data$internal_id
    data_x <-  model.matrixBayes(formula_x, data = data)
    randomize_dist_x <-  model.matrixBayes(formula_x, data = randomize_dist)
    data_x_s <-  model.matrixBayes(formula_x_s, data = data)
    target_dist_x_s <-  model.matrixBayes(formula_x_s, data = target_dist)
    all_equal <- length(all.vars(formula_x))
    all_equal_s <- length(all.vars(formula_x_s))

    # Prob in Target Distribution (use data_x_s)
    count_x0 <- data_x_s %*% t(target_dist_x_s)
    count_x  <- apply(count_x0, 1, function(i) sum(i == all_equal_s))
    target_prob <-  count_x/nrow(target_dist_x_s)

    # Prob in Randomizing Distribution (use data_x)
    count_x0_r <- data_x %*% t(randomize_dist_x)
    count_x_r  <- apply(count_x0_r, 1, function(i) sum(i == all_equal))
    randomize_prob <-  count_x_r/nrow(randomize_dist_x)
  }
  if(cross_int == TRUE){
    data0 <- data[order(data$pair_id),]
    side <- rep(c(1,0), times=nrow(data0)/2)
    data1 <- data0[side==1,]
    data2 <- data0[side==0,]
    new_id <- c(data1$internal_id, data2$internal_id)

    data1_x <-  model.matrixBayes(formula_x, data = data1)
    data2_x <-  model.matrixBayes(formula_x, data = data2)
    data1_x_s <-  model.matrixBayes(formula_x_s, data = data1)
    data2_x_s <-  model.matrixBayes(formula_x_s, data = data2)

    # combine
    data_x_s <- rbind(data1_x_s, data2_x_s)
    data_x <- rbind(data1_x, data2_x)
    data_x_rp <- rbind(data2_x, data1_x)

    # setup distributions
    randomize_dist_x <-  model.matrixBayes(formula_x, data = randomize_dist)
    target_dist_x <-  model.matrixBayes(formula_x, data = target_dist)
    target_dist_x_s <-  model.matrixBayes(formula_x_s, data = target_dist)
    all_equal <- length(all.vars(formula_x))
    all_equal_s <- length(all.vars(formula_x_s))

    # Prob in Target Distribution
    count_x0_t1 <- data_x_s %*% t(target_dist_x_s)
    count_x_t1  <- apply(count_x0_t1, 1, function(i) sum(i == all_equal_s))
    target_prob_t1 <-  count_x_t1/nrow(target_dist_x_s)

    count_x0_t2 <- data_x_rp %*% t(target_dist_x)
    count_x_t2  <- apply(count_x0_t2, 1, function(i) sum(i == all_equal))
    target_prob_t2 <-  count_x_t2/nrow(target_dist_x)
    target_prob <- target_prob_t1*target_prob_t2

    # Prob in Randomizing Distribution
    count_x0_r1 <- data_x %*% t(randomize_dist_x)
    count_x_r1  <- apply(count_x0_r1, 1, function(i) sum(i == all_equal))
    randomize_prob_r1 <-  count_x_r1/nrow(randomize_dist_x)

    count_x0_r2 <- data_x_rp %*% t(randomize_dist_x)
    count_x_r2  <- apply(count_x0_r2, 1, function(i) sum(i == all_equal))
    randomize_prob_r2 <-  count_x_r2/nrow(randomize_dist_x)

    randomize_prob <- randomize_prob_r1*randomize_prob_r2
  }

  design_w <- target_prob/randomize_prob

  design_w[design_w > 1000] <- 0
  design_w[design_w > 10] <- 10
  design_weight <- design_w[match(data$internal_id, new_id)]

  newdata <- data
  newdata$design_weight <- design_weight

  output <- list("design_weight" = design_weight, "newdata" = newdata)
  return(output)
}


weights_pAMCE_mar <- function(formula, factor_name, data, pair, pair_id, cross_int,
                              randomize_dist, target_dist){

  data$internal_id <- seq(1:nrow(data))
  if(pair == TRUE){data$pair_id <- pair_id}

  # HouseKeeping (Marginal Case)
  factor_use <- all.vars(formula)[-1]
  level_use  <- lapply(data[, factor_use], levels)

  if(setequal(factor_use, names(randomize_dist))  == FALSE){
    stop(" names(randomize_dist) should be the same as names of factors in 'formula' ")
  }
  if(setequal(factor_use, names(target_dist))  == FALSE){
    stop(" names(target_dist) should be the same as names of factors in 'formula' ")
  }

  randomize_dist <- randomize_dist[match(factor_use, names(randomize_dist))]
  target_dist <- target_dist[match(factor_use, names(target_dist))]
  for(z in  1:length(randomize_dist)){
    if(setequal(level_use[[z]], names(randomize_dist[[z]]))  == FALSE){
      stop(" level names of 'randomize_dist' should be the same as level names in 'data' ")
    }
    if(setequal(level_use[[z]], names(target_dist[[z]]))  == FALSE){
      stop(" level names of 'target_dist' should be the same as level names in 'data' ")
    }
    randomize_dist[[z]] <- randomize_dist[[z]][match(level_use[[z]], names(randomize_dist[[z]]))]
    target_dist[[z]] <- target_dist[[z]][match(level_use[[z]], names(target_dist[[z]]))]
  }


  # Case 1: When both are marginal
  weight_dist <- list()
  for(i in 1:length(factor_use)){
    if(factor_use[i] == factor_name){
      weight_dist[[i]] <- 1/randomize_dist[[factor_name]]
    }else{
      weight_dist[[i]] <- target_dist[[factor_use[i]]]/randomize_dist[[factor_use[i]]]
    }
    names(weight_dist[[i]]) <-  paste(factor_use[i], names(weight_dist[[i]]), sep = "")
  }
  if(cross_int == TRUE){
    weight_dist_cross <- list()
    for(i in 1:length(factor_use)){
      weight_dist_cross[[i]] <- target_dist[[factor_use[i]]]/randomize_dist[[factor_use[i]]]
      names(weight_dist_cross[[i]]) <-  paste(factor_use[i], names(weight_dist_cross[[i]]), sep = "")
    }
    weight_vec <- c(unlist(weight_dist), unlist(weight_dist_cross))
  }else{
    weight_vec <- unlist(weight_dist)
  }

  # combine with data
  if(cross_int == FALSE){
    data_x_exp <- model.matrixBayes(formula, data = data); data_x_exp[data_x_exp ==  0] <- NA
    new_id <- data$internal_id
  }else{
    data0 <- data[order(data$pair_id),]
    side <- rep(c(1,0), times=nrow(data0)/2)
    data1 <- data0[side==1,]
    data2 <- data0[side==0,]
    data1_x_exp <- model.matrixBayes(formula, data = data1)
    data2_x_exp <- model.matrixBayes(formula, data = data2)
    data_x_exp  <- rbind(cbind(data1_x_exp, data2_x_exp), cbind(data2_x_exp, data1_x_exp))
    data_x_exp[data_x_exp ==  0] <- NA
    new_id <- c(data1$internal_id, data2$internal_id)
  }
  data_x_exp <-  data_x_exp[, match(names(weight_vec), colnames(data_x_exp))]
  data_w <- t(t(data_x_exp)*weight_vec)
  design_w <- apply(data_w, 1, function(x) prod(x,  na.rm = TRUE))
  design_w[design_w > 1000] <- 0
  design_w[design_w > 10] <- 10
  design_weight <- design_w[match(data$internal_id, new_id)]

  newdata <- data
  newdata$design_weight <- design_weight

  output <- list("design_weight" = design_weight, "newdata" = newdata)
  return(output)
}
