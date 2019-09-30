#' Diagnose modeling assumptions for the model-based approach
#' @param x An object of class "pAMCE", a result of a call to 'model_pAMCE'.
#' @param factor_name A factor for which the function diagnoses modeling assumptions.
#' @export

diagnose_pAMCE <- function(x, factor_name){

  if(all(is.element(factor_name, names(x$AMCE))) == FALSE){
    stop(" 'factor_name' can only take a subset of factors estimated ")
  }

  # Diagnostic 1: Specification Test
  st <- list()
  for(i in 1:length(factor_name)){
    st_0 <- x$AMCE[[factor_name[i]]]
    st[[i]] <- st_0[st_0$type  %in% c("sample AMCE", "sample"), ]
  }
  names(st) <- factor_name

  # Diagnostic 2: Check Bootstrap Distributions
  if(x$input$reg  == FALSE){
    bd <- "no regularization"
  }else if(x$input$reg  == TRUE){
    table_AMCE <-  do.call("rbind", x$AMCE)
    table_AMCE_type <- setdiff(unique(table_AMCE$type), c("sample AMCE"))
    bd <- list()
    for(i in 1:length(factor_name)){
      ind_i <- which((table_AMCE$factor  == factor_name[i]) & (table_AMCE$type %in%  table_AMCE_type))
      name_i <- paste("Effect of ", table_AMCE[ind_i, "level"], " ~ ", table_AMCE[ind_i, "type"], sep = "")
      bd[[i]] <- x$boot_AMCE[ind_i, ]
      attr(bd[[i]], "level") <- table_AMCE[ind_i, "level"]
      attr(bd[[i]], "dist")  <- table_AMCE[ind_i, "type"]
    }
    names(bd) <- factor_name
  }

  out <- list("st" = st, "bd" = bd)

}

#' Plotting diagnostic checks
#' @param x An object of class "pAMCE", a result of a call to 'model_pAMCE.'
#' @param factor_name A factor for which the function diagnoses modeling assumptions.
#' @param legend_pos Position of the legend. Default is 'topright'.
#' @param target_dist_name Names of the target profile distributions to be used.
#' @param xlim Range for the x-axis.
#' @param mar Space on the left side of the plot. Default is 12.
#' @export

plot_diagnose <- function(x, factor_name, legend_pos = "topright", target_dist_name, xlim, mar){

  # reset parameters on exit
  opar <- par(no.readonly =TRUE)
  on.exit(par(opar))

  if(all(is.element(factor_name, names(x$AMCE))) == FALSE){
    stop(" 'factor_name' can only take a subset of factors estimated ")
  }

  d_result  <-  diagnose_pAMCE(x = x, factor_name =  factor_name)

  # Specification Test
  plot_pAMCE_base(x = x, factor_name = factor_name, legend_pos = legend_pos,
                  target_dist_name = c("sample AMCE", "sample"),
                  target_dist_name_use = c("Design-based", "Model-based"),
                  main = "Diagnostic 1: Specification Test", xlim = xlim, mar = mar)

  # Bootstrap Distribution
  bd <- d_result$bd
  if(x$input$reg ==  TRUE){
    total <- sum(unlist(lapply(1:length(bd), function(w) sum(attr(bd[[w]], "dist") %in% target_dist_name))))
    count <- 0
    boot <- ncol(bd[[1]])
    par(mar = c(4, 4, 6, 2))
    for(i in 1:length(factor_name)){
      keep <- which(attr(bd[[i]], "dist") %in% target_dist_name)
      at_dist <- attr(bd[[i]], "dist")[keep]
      at_level <- attr(bd[[i]], "level")[keep]
      bd[[i]]  <- bd[[i]][keep, , drop=FALSE]
      attr(bd[[i]], "dist")  <- at_dist
      attr(bd[[i]], "level")  <- at_level
      for(j in 1:nrow(bd[[i]])){
        count <-  count + 1
        par(ask = TRUE)
        plot(density(bd[[i]][j, ]), main = paste("Diagnostic 2: Check Bootstrap Distributions (", count ,"/",total,") \n",
                                                 "[factor, level] = [", names(bd)[i], ",", attr(bd[[i]], "level")[j], "]",
                                                 ",  ", attr(bd[[i]], "dist")[j], sep = ""),
             xlab = "", lwd = 2)

        if(count == 1 & boot < 500){
          message("Note: suggest 'boot' greater than 500 for final results.")
        }
      }
    }
  }
  par(ask = FALSE)
}
