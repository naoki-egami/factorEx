rm(list=ls())
library(factorEx)
data("OnoBurden")

# we focus on target profile distributions based on Democratic legislators.
# See de la Cuesta, Egami, and Imai (2019+) for details.
target_dist_marginal <- OnoBurden$target_dist_marginal

OnoBurden_data <- OnoBurden$OnoBurden_data # randomization based on uniform
# due to large sample size, focus on "congressional candidates" for this example
OnoBurden_data_cong <- OnoBurden_data[OnoBurden_data$office == "Congress", ]

# without pair_id: the function assume non-paird factorial experiment.
out_model <-
  model_pAMCE(formula = Y ~  gender + age + family + race + experience + party + pos_security,
              data = OnoBurden_data_cong,
              reg =  FALSE,
              cluster_id = OnoBurden_data_cong$id,
              target_dist  = target_dist_marginal, target_type = "marginal")
# ##################################################
# the same as the original function up to here
# ##################################################

# ##################################################
# New function is just this!
# ##################################################
mean_outcomes <- estimate_level(out_model)


# ##################################################
# Checking that the function is working properly
# ##################################################
# Check whether AMCE is the same as the difference in mean_outcomes
check <- mean_outcomes[mean_outcomes$type == "target_1", ]
uniq_fac <- unique(check$factor)

est_fun <- est_level <- list()
for(z in 1:length(uniq_fac)){
    AMCE_use <- out_model$AMCE[names(out_model$AMCE) == uniq_fac[z]][[1]]
    est_fun[[z]] <- AMCE_use[AMCE_use$type == "target_1", "estimate"]
    use_ind <- which(check$factor == uniq_fac[z])
    est_level[[z]] <- check[use_ind[-1], "estimate"] - check[use_ind[1], "estimate"]

    cat(z, all(round(est_fun[[z]],6) == round(est_level[[z]],6)))
}
