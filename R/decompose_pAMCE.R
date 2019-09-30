#' Decompose the difference between the pAMCEs
#' @param out An object of class "pAMCE", a result of a call to 'model_pAMCE'.
#' @param effect_name Effect for which the function decomposes the difference in the pAMCEs. The first element should be a factor name and the second element should be a level name.
#' @param target_diff Two target profile distributions for which the function compares the pAMCEs. If missing, the function comapres the first target profile distribution and the in-sample profile distribution.
#' @return \code{decompose_pAMCE} returns `data.frame` showing the decomposition of the difference between the pAMCEs.
#' @description See examples in `model_pAMCE`.
#' @export

decompose_pAMCE <- function(out, effect_name, target_diff){

  factor_name <- effect_name[1]
  level_name  <- effect_name[2]

  marginal_dist <- out$input$marginal_dist
  marginal_dist_u_list  <- out$input$marginal_dist_u_list
  marginal_dist_u_base  <- out$input$marginal_dist_u_base
  marginal_type  <- out$input$marginal_type
  cross_int <- out$input$cross_int

  if(missing(target_diff)){
    if(length(out$input$marginal_type) == 1){
      stop("Cannot compare two target profile distributions because `model_pAMCE` only includes one target profile distribution.")
    }else{
      target_diff <- out$input$marginal_type[c(2, 1)]
    }
  }

  if(all(is.element(target_diff, marginal_type)) == FALSE){
    stop("`target_diff` can only take names in `target_dist` specified in 'model_pAMCE'.")
  }

  if(any(out$coef_order >  2)){
    stop("Note: Decomposition is available only for two-way interaction models.")
  }

  if(all(is.element(factor_name, names(out$AMCE))) == FALSE){
    stop(" 'effect_name[1]' should be one of the factors estimated ")
  }

  if(length(factor_name) !=1 ){
    stop(" The function can take only one 'factor_name'")
  }

  boot_coef <- out$boot_coef
  use_name <- paste(factor_name, level_name, sep = "")


  # Estimate Marginal Contributions to AMEs ----------
  coef_focus <- boot_coef[, grep(use_name, colnames(boot_coef), fixed = T)]
  # vcov_focus <- vcovInt[grep(marginal_dist_u_base$level[m], rownames(vcovInt), fixed = T),
  #                       grep(marginal_dist_u_base$level[m], colnames(vcovInt), fixed = T)]
  if((length(coef_focus) > 0) == FALSE){
    stop(" 'effect_name[2]' is a name of the level and cannot take the baseline level or undefined levels of the specified factor")
  }
  estNames <- gsub(paste(use_name, ":", sep = ""), "", colnames(coef_focus), fixed = T)
  estNames <- gsub(paste(":", colnames(coef_focus)[1], sep = ""), "", estNames, fixed = T)

  if(cross_int == TRUE){
    estNames <- sub(paste(factor_name,"_rp", sep = ""), factor_name, estNames)
  }

  table_AME_diff <- c()
  # For each marginal distribution,
  ind1 <- which(marginal_type == target_diff[1])
  ind2 <- which(marginal_type == target_diff[2])
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
    AME_diff$type <- paste(target_diff[1], target_diff[2], sep = " - ")
    AME_diff$factor   <- marginal_factor[m]
    AME_diff$estimate <- coef_AME_diff;
    AME_diff$se <- se_AME_diff
    AME_diff$low.95ci  <- low_AME_diff
    AME_diff$high.95ci <- high_AME_diff
    table_AME_diff <- rbind(table_AME_diff, AME_diff)
  }
  class(table_AME_diff)  <- c(class(table_AME_diff), "decompose")
  return(table_AME_diff)
}

#' Plot decomposition of the difference between pAMCEs
#' @param x An object of class "pAMCE", a result of a call to 'model_pAMCE'.
#' @param effect_name Effect for which the function decomposes the difference in the pAMCEs. The first element should be a factor name and the second element should be a level name.
#' @param target_diff Two target profile distributions for which the function compares the pAMCEs. If missing, the function comapres the first target profile distribution and the in-sample profile distribution.
#' @param mar Space on the left side of the plot. Default is 12.
#' @description  See examples in `model_pAMCE`.
#' @export

plot_decompose <- function(x, effect_name, target_diff, mar = 12){

  # reset parameters on exit
  opar <- par(no.readonly =TRUE)
  on.exit(par(opar))

  factor_name <- effect_name[1]
  level_name  <- effect_name[2]

  cross_int <- x$input$cross_int


  if(missing(target_diff)){
    if(length(x$input$marginal_type) == 1){
      stop("Cannot compare two target profile distributions because `model_pAMCE` only includes one target profile distribution.")
    }else{
      target_diff <- x$input$marginal_type[c(2, 1)]
    }
  }

  dec_tab <- decompose_pAMCE(out = x, effect_name = effect_name, target_diff = target_diff)

  point <- rev(dec_tab[,3])
  low  <- rev(dec_tab[,5])
  high   <- rev(dec_tab[,6])
  fac_name_p <- rev(dec_tab[,2])

  xmin <- min(low); xmax <- max(high)

  par(mar = c(4, mar, 6, 4))
  plot(point, seq(1:nrow(dec_tab)), pch = 19, ylim = c(0.5, nrow(dec_tab)+0.5), yaxt = "n",
       xlim = c(xmin, xmax), ylab =  "", xlab = "Difference in Popuation AMCEs",
       main = paste("Decompose Difference in Population AMCEs:\n", unique(dec_tab$type), sep =""))
  segments(low, seq(1:nrow(dec_tab)), high, seq(1:nrow(dec_tab)), lwd = 2)
  abline(v = 0, lty = 2)
  Axis(side = 2, at = seq(1:nrow(dec_tab)), labels = fac_name_p, las = 1)

  return(dec_tab)
}
