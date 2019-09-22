rm(list = ls())
load("../data/ono_data_analysis.RData")

load("../data/OnoBurden_data_pr.rdata")
OnoBurden_data_pr <- Pop.all
## OnoBurden_data_pr <- Pop.all[is.element(Pop.all$id, seq(1:500)),]

# setup the main data
OnoBurden_data <- dfOno
use_id <- sample(dfOno$id, size = 500, replace = FALSE)
# OnoBurden_data <- OnoBurden_data[OnoBurden_data$office == "Congress", ]
# OnoBurden_data <- dfOno[dfOno$id %in% use_id,]

formula_u <-  Y ~ gender + age + family + race + experience + party + pos_security
# formula_u <- as.formula(Y ~ gender + age + race + trait + party + family + fav_rating + experience + policy_expertise +
#                           pos_security + pos_immigrants + pos_abortion + pos_deficit)

## setup marginal
base <- marUse[[3]]

name_u <-  names(base)
target_dist_marginal  <- list()
for(i in  1:length(base)){
  target_dist_marginal[[i]] <- base[[i]]
}
names(target_dist_marginal) <- name_u
target_dist_marginal <- target_dist_marginal[all.vars(formula_u)[-1]]


target_data <- dfJoint[, is.element(colnames(dfJoint), all.vars(formula_u)) == TRUE]
factor_use <- all.vars(formula_u)[-1]
target_data_rep <- target_data[target_data$party == "Rep", ]
target_data_dem <- target_data[target_data$party == "Dem", ]
# target_data_dem$trait <- dfOno$trait[1:nrow(target_data_dem)]

## setup target_data
target_dist_data <- target_data_dem[, factor_use]

## setup partial Joints
## (gender, age, family)
g_a_f <- table(target_dist_data[, c("gender", "age", "family")])/nrow(target_dist_data)
## c("experience", "pos_security")
e_p <- table(target_dist_data[, c("experience", "pos_security")])/nrow(target_dist_data)

target_dist_partial <-
  list(g_a_f,
       target_dist_marginal[[c("race")]],
       target_dist_marginal[["party"]],
       e_p)
names(target_dist_partial) <- c("gender:age:family", c("race", "party"), "experience:pos_security")

# ## setup partial Joints
# ## (gender, age, family)
# g_a_f <- table(target_dist_data[, c("gender", "age", "family")])/nrow(target_dist_data)
# ## ("fav_rating", "experience")
# f_e <- table(target_dist_data[, c("fav_rating", "experience")])/nrow(target_dist_data)
#
# target_dist_partial <-
#   list(g_a_f, target_dist_marginal[[c("race")]],
#        target_dist_marginal[["trait"]],
#        target_dist_marginal[["party"]],
#        f_e,
#        target_dist_marginal[["policy_expertise"]],
#        target_dist_marginal[["pos_security"]],
#        target_dist_marginal[["pos_immigrants"]],
#        target_dist_marginal[["pos_abortion"]],
#        target_dist_marginal[["pos_deficit"]])
# names(target_dist_partial) <- c("gender:age:family", c("race", "trait", "party"),
#                                 "fav_rating:experience",
#                                 c("policy_expertise", "pos_security",
#                                   "pos_immigrants", "pos_abortion", "pos_deficit"))

OnoBurden <- list("OnoBurden_data"  = OnoBurden_data, "OnoBurden_data_pr"  = OnoBurden_data_pr,
                  "target_dist_marginal" = target_dist_marginal,
                  "target_dist_partial"  = target_dist_partial,
                  "target_dist_data"  = target_dist_data)
# save(OnoBurden, file = "../data/OnoBurden.rdata")
