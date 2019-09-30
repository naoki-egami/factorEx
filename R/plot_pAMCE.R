#' Plotting the estimated population AMCEs
#' @param x An object of class "pAMCE", a result of a call to 'model_pAMCE' or 'design_pAMCE'
#' @param factor_name Factors for which the function visualizes the pAMCEs
#' @param target_dist_name Names of the target profile distributions to be used
#' @param legend_pos Position of the legend. Default is 'topright'
#' @param main Title of the plot
#' @param xlim Range for the x-axis
#' @param mar Space on the left side of the plot. Default is 12
#' @param diagnose Whether we plot diagnostic checks recommended in de la Cuesta, Egami, and Imai (2019). Default is FALSE
#' @param ... Other graphical parameters
#' @export

plot.pAMCE <- function(x, factor_name, target_dist_name,
                       legend_pos = "topright",
                       main = "Estimated population AMCEs",
                       xlim, mar = 12, diagnose = FALSE, ...){

  # reset parameters on exit
  opar <- par(no.readonly =TRUE)
  on.exit(par(opar))

  if(missing(factor_name) == TRUE){
    factor_name  <- names(x$AMCE)[1:2]
  }else{
    if(all(is.element(factor_name, names(x$AMCE))) == FALSE){
      stop(" 'factor_name' can only take a subset of factors estimated ")
    }
  }

  if(x$approach  == "design_based"){
    diagnose <- FALSE
  }

  if(missing(target_dist_name) == FALSE){
    if(all(is.element(target_dist_name, unique(x$AMCE[[1]]$type))) == FALSE){
      stop(" 'target_dist_name' can only take a subset of target profile distributions used. ")
    }
    pl_l <- length(target_dist_name)
  }else{
    if(x$approach  == "model_based"){
      target_dist_name <- setdiff(unique(x$AMCE[[1]]$type), "sample")
      pl_l <- length(unique(x$AMCE[[1]]$type)) - 1
    }
    if(x$approach  == "design_based"){
      target_dist_name <- "target"
      pl_l <- 1
    }
  }

  col <- palette()[1:pl_l]
  pch <- rep(19, pl_l)

  if(missing(mar)){
    mar <- 12
  }

  plot_name <- target_dist_name
  # type_difference <- setdiff(x$type_difference, "sample AMCE")
  # plot_name[plot_name == "Sample"] <- "sample AMCE"

  ## Correct all esimates ----------
  p_coef_full <- p_high_full <- p_low_full <- c()
  p_name_full <- p_name_f_full <- p_col_full <- p_pch_full <- c()
  for(g in 1:length(factor_name)){
    p_AME  <- x$AMCE[[factor_name[g]]]
    p_AME <- p_AME[order(factor(p_AME$level, levels = unique(p_AME$level))),] # updated on 12/27 (Naoki)

    # p_AME   <- p_AME[(p_AME$type %in% type_difference) == FALSE, ]

    if(missing(target_dist_name) == FALSE){
      p_AME  <- p_AME[(p_AME$type %in% target_dist_name) == TRUE, ]
      p_AME  <- p_AME[order(factor(p_AME$level, levels = unique(p_AME$level)),
                            factor(p_AME$type, levels = target_dist_name)),]
    }

    p_coef <- c(NA, p_AME$estimate)
    p_se   <- c(NA, p_AME$se)
    p_high <- p_coef + 1.96*p_se
    p_low  <- p_coef - 1.96*p_se
    p_name_t <- paste(factor_name[g], " (", x$baseline[factor_name[g]], "):   ", sep="")
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
  if(missing(target_dist_name) == FALSE) p_type <- target_dist_name
  if(missing(xlim)){
    p_x <- c(min(p_low_full, na.rm=TRUE), max(p_high_full, na.rm = TRUE))
  }else{
    p_x <- xlim
  }

  cex <- 1

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
  abline(v = 0, lty = 2)
  if(length(plot_name) > 1){
    if(is.character(legend_pos[1]) == TRUE)  legend(legend_pos, plot_name, col= col, pch = pch)
    if(is.character(legend_pos[1]) == FALSE) legend(x=legend_pos[1], y=legend_pos[2],
                                                    plot_name, col= col, pch = pch)
  }

  if(diagnose ==  TRUE){
    par(ask = TRUE)
    plot_diagnose(x = x, factor_name = factor_name, legend_pos = legend_pos,
                  target_dist_name = target_dist_name, xlim =  xlim, mar = mar)
  }
}

plot_pAMCE_base <- function(x, factor_name, target_dist_name,target_dist_name_use,
                            legend_pos = "topright",
                            main = "Estimated population AMCEs",
                            xlim, mar = 12){

  if(all(is.element(factor_name, names(x$AMCE))) == FALSE){
    stop(" 'factor_name' can only take a subset of factors estimated ")
  }

  if(missing(target_dist_name) == FALSE){
    pl_l <- length(target_dist_name)
  }else{
    target_dist_name <- unique(x$AMCE[[1]]$type)
    pl_l <- length(unique(x$AMCE[[1]]$type))
  }

  col <- palette()[1:pl_l]
  pch <- rep(19, pl_l)

  if(missing(mar)){
    mar <- 12
  }

  if(missing(target_dist_name_use)){
    target_dist_name_use <- target_dist_name
  }

  plot_name <- target_dist_name_use
  type_difference <- setdiff(x$type_difference, "sample AMCE")

  ## Correct all esimates ----------
  p_coef_full <- p_high_full <- p_low_full <- c()
  p_name_full <- p_name_f_full <- p_col_full <- p_pch_full <- c()
  for(g in 1:length(factor_name)){
    p_AME  <- x$AMCE[[factor_name[g]]]
    p_AME <- p_AME[order(factor(p_AME$level, levels = unique(p_AME$level))),] # updated on 12/27 (Naoki)

    p_AME   <- p_AME[(p_AME$type %in% type_difference) == FALSE, ]

    if(missing(target_dist_name) == FALSE){
      p_AME  <- p_AME[(p_AME$type %in% target_dist_name) == TRUE, ]
      p_AME  <- p_AME[order(factor(p_AME$level, levels = unique(p_AME$level)),
                            factor(p_AME$type, levels = target_dist_name)),]
    }

    p_coef <- c(NA, p_AME$estimate)
    p_se   <- c(NA, p_AME$se)
    p_high <- p_coef + 1.96*p_se
    p_low  <- p_coef - 1.96*p_se
    p_name_t <- paste(factor_name[g], " (", x$baseline[factor_name[g]], "):   ", sep="")
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
  if(missing(target_dist_name) == FALSE) p_type <- target_dist_name
  if(missing(xlim)){
    p_x <- c(min(p_low_full, na.rm=TRUE), max(p_high_full, na.rm = TRUE))
  }else{
    p_x <- xlim
  }

  cex <- 1

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
  abline(v = 0, lty = 2)
  if(length(plot_name) > 1){
    if(is.character(legend_pos[1]) == TRUE) legend(legend_pos, plot_name, col= col, pch = pch)
    if(is.character(legend_pos[1]) == FALSE) legend(x=legend_pos[1], y=legend_pos[2],
                                                    plot_name, col= col, pch = pch)
  }
}
