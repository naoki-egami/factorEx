#' #' #' Estimating the population AMCE using the design-based estimator
#' #' #' @param formula Formula
#' #' #' @export
#' #'
#' #'
#' #'
#'
#'
#' weightDesign <- function(factor_name, target_dist, randomize_dist, sample, sim_size = 1000){
#'
#'   # we average over factor levels
#'   # The standard case:
#'   ## For each potential outcomes, you use 1/D_ell x N sample
#'   ## We sume them over D_ell, so we get to use N sample
#'
#'   # Simulation-based Calculation
#'   # Only Marginal Distributions
#'
#'   # check distributions carefully
#'
#'   # they should have the same order
#'   randomize_dist <- randomize_dist[match(names(target_dist), names(randomize_dist))]
#'
#'   # Compute weights
#'   all_fac <- names(target_dist)
#'   weight_dist <- list()
#'   for(i in 1:length(all_fac)){
#'     if(all_fac[i] == factor_name){
#'       weight_dist[[i]] <- 1/randomize_dist[[factor_name]]
#'     }else{
#'       weight_dist[[i]] <- target_dist[[all_fac[i]]]/randomize_dist[[all_fac[i]]]
#'     }
#'   }
#'   names(weight_dist) <- all_fac
#'   weight_vec <- unlist(weight_dist)
#'
#'   # Simulate Data
#'   seq_level <- lapply(target_dist, function(x) seq(1:length(x)))
#'   ss <- c()
#'   cat("Simulation:")
#'   for(s in 1:sim_size){
#'     sample_w <- rep(1, times = sample)
#'     for(i in 1:length(all_fac)){
#'       sample_w0 <- weight_dist[[i]][sample(seq_level[[i]], size = sample, replace =  TRUE, prob = randomize_dist[[i]])]
#'       sample_w  <-  sample_w*sample_w0
#'     }
#'     sample_w[sample_w > 10] <- 10
#'     ss[s] <-  (sum(sample_w)^2)/sum(sample_w^2)
#'     if((s %% ((sim_size)/5)) == 0){
#'       per <- s/((sim_size)/5)*20
#'       cat(paste(per, "%..", sep  = ""))
#'     }
#'   }
#'   effective_ss <- mean(ss)
#'
#'   out <-  list("effective_ss" = effective_ss, "ss" = ss)
#'   return(out)
#' }
#'
#'
#'
#'
#' # weightDesign <- function(factor_name, target_dist, randomize_dist, sample, pair, sim_size = 1000){
#' #
#' #   # check distributions carefully
#' #
#' #   # they should have the same order
#' #   randomize_dist <- randomize_dist[match(names(target_dist), names(randomize_dist))]
#' #
#' #   # Compute weights
#' #   all_fac <- names(target_dist)
#' #   weight_dist <- list()
#' #   for(i in 1:length(all_fac)){
#' #     if(all_fac[i] == factor_name){
#' #       weight_dist[[i]] <- randomize_dist[[factor_name]]
#' #     }else{
#' #       weight_dist[[i]] <- target_dist[[all_fac[i]]]/randomize_dist[[all_fac[i]]]
#' #     }
#' #   }
#' #   names(weight_dist) <- all_fac
#' #   weight_vec <- unlist(weight_dist)
#' #
#' #   # Simulate Data
#' #   all_levels <- lapply(target_dist, function(x) names(x))
#' #   ss <- c()
#' #   for(s in 1:sim_size){
#' #     data_list <- as.data.frame(matrix(NA, ncol = 0, nrow = sample))
#' #     for(i in 1:length(all_fac)){
#' #       data_list[, all_fac[i]] <- factor(sample(all_levels[[i]], size = sample,
#' #                                                replace =  TRUE,
#' #                                                prob = randomize_dist[[i]]), levels = all_levels[[i]])
#' #     }
#' #     data_x <- model.matrixBayes(~., data = data_list); data_x[data_x ==  0] <- NA
#' #     data_w <- t(t(data_x)*weight_vec)
#' #     sample_w <- apply(data_w, 1, function(x) prod(x,  na.rm = TRUE))
#' #
#' #     ss[s] <-  (sum(sample_w)^2)/sum(sample_w^2)
#' #   }
#' #   effective_ss <- mean(ss)
#' #
#' #   out <-  list("effective_ss" = effective_ss, "ss" = ss)
#' #   return(out)
#' # }
