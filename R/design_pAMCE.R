weights_pAMCE <- function(factor_name, data, pair_id,
                          randomize_dist, randomize_type = "marginal",
                          target_dist, target_type = "marginal"){

  # HouseKeeping
  randomize_dist <- randomize_dist[match(names(target_dist), names(randomize_dist))]

  # Case 1: When both are marginal
  all_fac <- names(target_dist)
  weight_dist <- list()
  for(i in 1:length(all_fac)){
    if(all_fac[i] == factor_name){
      weight_dist[[i]] <- 1/randomize_dist[[factor_name]]
    }else{
      weight_dist[[i]] <- target_dist[[all_fac[i]]]/randomize_dist[[all_fac[i]]]
    }
  }
  names(weight_dist) <- all_fac
  weight_vec <- unlist(weight_dist)

}
