#' Estimating PAMCE withithout regularization
#' @param formula formula
#' @param data data
#' @param pair whether we use the paired-conjoint design
#' @param pair_id id for paired-conjoint design. required when 'pair = TRUE'
#' @export


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
