

randomize_dist <-  exp_marginal
target_dist <- exp_marginal
formula = as.formula(formula_u)
data = dfOnoRep
data$pair_id <- dfOnoRep$pair_id
cross_int <- TRUE

setequal(names(marginal_rep), all.vars(as.formula(formula_u)) [-1])

name_u <-  names(marginal_dem)
marginal_dem_l  <- list()
for(i in  1:length(marginal_dem)){
  marginal_dem_l[[i]] <- marginal_dem[[i]]
}
names(marginal_dem_l) <- name_u

library(estimatr)
design_out <- design_pAMCE(formula = as.formula(formula_u),
                           data = dfOnoRep, factor_name = "gender",
                           pair_id = dfOnoRep$pair_id,
                           cluster_id = dfOnoRep$id,
                           target_dist  = marginal_dem_l,
                           target_type = "marginal")
design_out$AMCE$gender

model_out <- model_pAMCE(formula = as.formula(formula_u),
                         data = dfOnoRep, reg = FALSE,
                         pair = TRUE, pair_id = dfOnoRep$pair_id,
                         cluster = dfOnoRep$id,
                         target_dist = target_dist,
                         target_type = c("marginal", "marginal" , "target_data"),
                         boot = 100)
model_out$AMCE$gender



# $gender
# type factor level    estimate         se     low.95ci  high.95ci
# 1      target gender  Male -0.01149564 0.05783735 -0.125982632 0.10299136
# 2 sample AMCE gender  Male  0.01841705 0.01041972 -0.002005604 0.03883971


target_marginal <- createDist(formula_u, target_data = target_data[, factor_use],
                              exp_data  = dfOnoRep, type = "marginal")
target_joint <- createDist(formula_u, target_data = target_data[, factor_use],
                           exp_data  = dfOnoRep, type = "joint")

target_marginal <- createDist(formula_u, target_data = dfOnoRep,
                              exp_data  = dfOnoRep, type = "marginal")
target_joint <- createDist(formula_u, target_data = dfOnoRep,
                           exp_data  = dfOnoRep, type = "joint")


name_use <-  c()
target_dist <- list()
for(i in 1:length(partial_joint_name)){
  if(length(partial_joint_name[[i]])  >  1){
    name_use[i] <- paste(partial_joint_name[[i]],  collapse =  ":")
    dim1 <- length(unique(target_joint[[name_use[i]]][,3]))
    dim2 <- length(unique(target_joint[[name_use[i]]][,4]))
    dim1_name <- unique(target_joint[[name_use[i]]][,3])
    dim2_name <- unique(target_joint[[name_use[i]]][,4])
    target_dist[[i]]  <- array(target_joint[[name_use[i]]][, 5],
                               dim = c(dim1, dim2),
                               dimnames  = list(dim1_name,  dim2_name))
  }else{
    name_use[i] <- partial_joint_name[[i]]
    target_dist[[i]]  <- target_marginal[[name_use[i]]]
  }
}

partial_joint_name
target_dist

randomize_dist <- list()
for(i in 1:length(partial_joint_name)){
  randomize_dist[[i]]  <- table(data_u[,  partial_joint_name[[i]]])/nrow(data_u)
}

design_p <- weights_pAMCE_partial(formula = as.formula(formula_u),
                                  factor_name  = "gender",
                                  data = dfOnoRep, pair_id = dfOnoRep$pair_id,
                                  randomize_dist = randomize_dist, randomize_type = "marginal",
                                  target_dist  = target_dist, target_type = "marginal",
                                  partial_joint_name = partial_joint_name)


library(estimatr)
a <- lm_robust(Y  ~ gender, data  =  design_p$newdata, weights = design_weight)
a_2 <- lm_robust(Y  ~ gender, data  =  data)

summary(a)
summary(a_2)
wdm(y = data$Y, treat = as.numeric(data$gender == "Male"), weights = data$design_weight)
summary(ameOut, factor_name = "gender")

wdm <- function(y, treat, weights){
  est <- sum(y[treat == 1]*weights[treat ==  1])/sum(weights[treat == 1]) -
    sum(y[treat == 0]*weights[treat ==  0])/sum(weights[treat == 0])
  return(est)
}

summary(a)
summary(a_2)

summary(ameOut,  factor_name = "gender")

formula_u <- as.formula(Y ~ gender + age + race + trait + party + family + fav_rating + experience + policy_expertise + pos_security + pos_immigrants +
                          pos_abortion + pos_deficit)
target_data <- dfJoint[, is.element(colnames(dfJoint), all.vars(formula_u)) == TRUE]
levels(target_data$age)  <-  levels(dfOnoRep$age)
levels(target_data$race) <-  levels(dfOnoRep$race)
levels(target_data$experience)   <- levels(dfOnoRep$experience)
levels(target_data$pos_security) <- levels(dfOnoRep$pos_security)
levels(target_data$pos_deficit)  <- levels(dfOnoRep$pos_deficit)
target_data$trait <- data$trait[1:nrow(target_data)]
target_data_rep <- target_data[target_data$party == "Rep", ]
target_data_dem <- target_data[target_data$party == "Dem", ]

data_u <- data[, factor_use]
target_dist <- target_data[, factor_use]

formula_u2 <- as.formula(~ gender + age + race + trait + party + family + fav_rating + experience + policy_expertise + pos_security + pos_immigrants +
                           pos_abortion + pos_deficit)

data_u_x <-  model.matrixBayes(formula_u2, data = data_u)
target_data_u_x <-  model.matrixBayes(formula_u2, data = target_data_u)
all.vars(formula_u2)
all_equal <- 13

count_x0 <- data_u_x %*% t(target_data_u_x)
count_x  <- apply(count_x0, 1, function(i) sum(i == 13))
design_prob <-  count_x/nrow(target_data_u_x)
design_w <- 1/design_prob
design_w[design_w > 1000] <- 0
design_w[design_w > 10] <- 10


library(fastLink)
a <- fastLink(data_u, target_data_u, varnames = factor_use)

randomize_dist <- data[, factor_use]
target_dist <- target_data[, factor_use]

formula_u <- Y ~ age + family + race + experience + trait + party +
  policy_expertise + pos_security + pos_immigrants + pos_abortion +
  pos_deficit + fav_rating + gender

# Create Experimental Joint Distribution
target_marginal <- createDist(formula_u, target_data = target_data[, factor_use],
                              exp_data  = dfOnoRep, type = "marginal")
target_joint <- createDist(formula_u, target_data = target_data[, factor_use],
                           exp_data  = dfOnoRep, type = "joint")

partial_joint_name  <- list(c("age", "family"), "race",  c("experience", "trait"),
                            c("party", "pos_immigrants"), c("policy_expertise", "pos_security"),
                            "pos_abortion", "pos_deficit", "fav_rating", "gender")
setequal(factor_use, unlist(partial_joint_name))

name_use <-  c()
target_dist <- list()
for(i in 1:length(partial_joint_name)){
  if(length(partial_joint_name[[i]])  >  1){
    name_use[i] <- paste(partial_joint_name[[i]],  collapse =  ":")
    dim1 <- length(unique(target_joint[[name_use[i]]][,3]))
    dim2 <- length(unique(target_joint[[name_use[i]]][,4]))
    dim1_name <- unique(target_joint[[name_use[i]]][,3])
    dim2_name <- unique(target_joint[[name_use[i]]][,4])
    target_dist[[i]]  <- array(target_joint[[name_use[i]]][, 5],
                               dim = c(dim1, dim2),
                               dimnames  = list(dim1_name,  dim2_name))
  }else{
    name_use[i] <- partial_joint_name[[i]]
    target_dist[[i]]  <- target_marginal[[name_use[i]]]
  }
}

partial_joint_name
target_dist

randomize_dist <- list()
for(i in 1:length(partial_joint_name)){
  randomize_dist[[i]]  <- table(data_u[,  partial_joint_name[[i]]])/nrow(data_u)
}



sum(table(data_u$age, data_u$family)/nrow(data_u))

array(exp_joint[[1]][,5], dim =  c(2, 4, 2))

weights_pAMCE_partial(formula = as.formula(formula_u),
                      factor_name  = "gender",
                      data, pair_id,
                      randomize_dist, randomize_type = "marginal",
                      target_dist, target_type = "marginal",
                      partial_joint_name)
