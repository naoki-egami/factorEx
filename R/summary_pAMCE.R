#' Summarizing the estimated population AMCEs
#' @param object An object of class "pAMCE", a result of a call to 'model_pAMCE' or 'design_pAMCE'
#' @param factor_name Factors for which the function visualizes the pAMCEs
#' @param target_dist_name Names of the target profile distributions to be used
#' @param sample Whether to print the sample AMCEs, which is estimated based on the profile distribution used for randomization
#' @param ... Other parameters
#' @export

summary.pAMCE <- function(object, factor_name, target_dist_name, sample = FALSE, ...){

  # #################
  # Housekeeping
  # #################
  if(missing(factor_name)){
    factor_name <- all.vars(object$input$formula)[-1]
  }else{
    if(all(is.element(factor_name, names(object$AMCE))) == FALSE){
      stop(" 'factor_name' can only take a subset of factors estimated ")
    }
  }

  if(missing(target_dist_name)){
    if(object$approach  == "model_based")  target_dist_name <- names(object$input$target_dist)
    if(object$approach  == "design_based") target_dist_name <- "target"
  }else{
    if(all(is.element(target_dist_name, unique(object$AMCE[[1]]$type))) == FALSE){
      stop(" 'target_dist_name' can only take a subset of target profile distributions used. ")
    }
  }
  if(sample == TRUE){
    target_dist_name <- c(target_dist_name, "sample AMCE")
  }

  AMCE_table <- do.call("rbind", object$AMCE)

  AMCE_table_u <- AMCE_table[AMCE_table$type %in% target_dist_name, ]
  AMCE_table_u <- AMCE_table_u[AMCE_table_u$factor %in% factor_name, ]
  AMCE_table_u$pv <- round(2*(1- pnorm(abs(AMCE_table_u$estimate/AMCE_table_u$se))), 3)
  AMCE_print <- AMCE_table_u[, c(1, 2, 3, 4, 5, 8)]
  colnames(AMCE_print) <- c("target_dist", "factor", "level", "Estimate", "Std. Error", "p value")

  # add significance
  sig <- rep("", length(AMCE_print[,6]))
  sig[AMCE_print[ , 6]  < 0.001] <- "***"
  sig[AMCE_print[ , 6]  >= 0.001 & AMCE_print[ , 6]  < 0.01] <- "**"
  sig[AMCE_print[ , 6]  >= 0.01 & AMCE_print[ , 6]  < 0.05] <- "*"
  sig[AMCE_print[ , 6]  >= 0.05 & AMCE_print[ , 6]  < 0.1] <- "."
  AMCE_print <- as.data.frame(AMCE_print)
  AMCE_print$Sig <- sig
  colnames(AMCE_print)[7] <- ""

  cat("\n")
  cat("----------------\n")
  cat("Population AMCEs:\n")
  cat("----------------\n")
  print(AMCE_print, row.names = FALSE)
  cat("---\n")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")

  rownames(AMCE_print) <- NULL
  invisible(AMCE_print)
}
