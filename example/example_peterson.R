rm(list=ls())
#dev.off() # Clear plot window so the only thing that is saved at end of scripts are plots that are run in this script
basePar <- par() # Save this for later to use before plotting AMEs
plotNames <- c()

require(haven)
require(PopCon)
require(genlasso)
require(FindIt)
require(sjlabelled)
require(magrittr)
require(tidyverse)
source("example/example_code_internal.R")


## TO-DO: Make separate distributions by gender
############ Peterson (2017) ######################

load("example/estimation.data.conjoint.RData")
dfExp <- estimation.data.conjoint; rm(estimation.data.conjoint)
#dfExp <- dfExp[dfExp$information.condition <= 7, ] # Remove the ones where only one level was shown
## Do some misc variable generation

dfExp$party.candidate <- factor(ifelse(dfExp$party.id == "Democrat" & dfExp$copartisan == 1, "Democrat",
                                ifelse(dfExp$party.id == "Republican" & dfExp$copartisan == 1, "Republican",
                                       ifelse(dfExp$party.id == "Democrat" & dfExp$copartisan == 0, "Republican", "Democrat"))))
dfExp$copartisan <- factor(dfExp$copartisan)
dfExp$id <- as.character(dfExp$participant.id)
dfExp[dfExp$id[which(table(dfExp$id) != 6)], ] # Should be no rows
dfExp$profile <- dfExp$candidate.order
dfExp$task <- dfExp$choice.task.number
dfExp$pair_id <- paste(dfExp$id, dfExp$task, sep = "-")



facVars <-  c("gender.candidate", "profession.candidate", "age.candidate",
              "family.status.candidate", "race.candidate", "military.service.candidate",
              "education.candidate", "abortion.stance.candidate", "spending.stance.candidate",
              "party.candidate", "copartisan")
facVars <-  c("gender.candidate", "profession.candidate", "age.candidate",
              "family.status.candidate", "race.candidate", "military.service.candidate",
              "education.candidate", "abortion.stance.candidate", "spending.stance.candidate",
              "party.candidate", "copartisan")
facNames <- c("Gender", "Profession", "Age", "Family Status", "Race", "Military Services",
              "Education", "Abortion", "Spending", "Party", "Copartisan")
lapply(dfExp[, facVars], levels)

dfNew <- dfExp
facOrders <- list(c(2, 3, 1), c(6, 5, 3, 2, 4, 7, 1),
                  c(2:10, 1), c(3, 2, 4, 8, 7, 5, 6, 1),
                  c(5, 2, 3, 4, 1), c(2, 3, 1),
                  c(6, 2, 5, 4, 3, 1), c(2, 4, 3, 1),
                  c(4, 5, 6, 3, 2, 1), c(1, 2), c(1, 2))

for(i in 1:length(facVars)){
  dfNew[, facVars[i]] <- droplevels(factor(dfNew[, facVars[i]]))
  dfNew[, facVars[i]] <- factor(dfNew[, facVars[i]],
                                levels = levels(dfNew[, facVars[i]])[facOrders[[i]]])
}


lapply(dfNew[, facVars], levels)
lapply(dfNew[, facVars], FUN = function(x) prop.table(table(x)))

## Declare experimental variables

expVars <- facVars
depVar <- "chosen.candidate"
idVars <- c("id", "profile", "task", "pair_id")

## Format variables correctly
dfNew[, expVars][sapply(dfNew[, expVars], is.numeric)] <- lapply(dfNew[, expVars][sapply(dfNew[, expVars],
                                                                                         is.numeric)], as.factor)

lapply(dfNew[, facVars], levels)

## Create formula and experimental distribution


require(openxlsx)
rmVars <- c("copartisan")
expEq <- as.formula(make.equation(depVar, expVars[!expVars %in% rmVars]))
expDist <- createDist(expEq, data = dfNew)
#write.csv(expDist, createPath("peterson_expdist.csv", "Re-analysis/data"))

repHigh <- read.xlsx("example/peterson_popdist_all.xlsx", sheet = 1, cols = 1:3)
repMed <- read.xlsx("example/peterson_popdist_all.xlsx", sheet = 2, cols = 1:3)
repLow <- read.xlsx("example/peterson_popdist_all.xlsx", sheet = 3, cols = 1:3)
demHigh <- read.xlsx("example/peterson_popdist_all.xlsx", sheet = 4, cols = 1:3)
demMed <- read.xlsx("example/peterson_popdist_all.xlsx", sheet = 5, cols = 1:3)
demLow <- read.xlsx("example/peterson_popdist_all.xlsx", sheet = 6, cols = 1:3)
genHigh <- read.xlsx("example/peterson_popdist_all.xlsx", sheet = 7, cols = 1:3)
genMed <- read.xlsx("example/peterson_popdist_all.xlsx", sheet = 8, cols = 1:3)
genLow <- read.xlsx("example/peterson_popdist_all.xlsx", sheet = 9, cols = 1:3)

## Create pure uniform distribution

facInd <- unique(repHigh$factor)
unifPure <- repHigh
unifPure$prop <- unlist(tapply(repHigh$prop, repHigh$factor,
                               function(x) rep(1/length(x), length(x)))[facInd])

## Create separate uniform distributions for high, low and medium visibilities

## High visibility
facInd <- unique(repHigh$factor)
unifHigh <- repHigh
unifHigh$prop[unifHigh$levels != "NOT.SEEN"] <- unlist(tapply(repHigh$prop[unifHigh$levels != "NOT.SEEN"],
                                                              repHigh$factor[unifHigh$levels != "NOT.SEEN"],
                                                              function(x) rep(1/length(x), length(x)))[facInd])
for(i in 1:length(facInd)){
  if(any(unifHigh$levels[unifHigh$factor == facInd[i]] %in% c("NOT.SEEN"))){
    notProp <- unifHigh$prop[unifHigh$factor == facInd[i] & unifHigh$levels == "NOT.SEEN"]
    levProp <- unifHigh$prop[unifHigh$factor == facInd[i] & unifHigh$levels != "NOT.SEEN"]
    unifHigh$prop[unifHigh$factor == facInd[i] & unifHigh$levels != "NOT.SEEN"] <- levProp*(1 - notProp)
  }
}

## Medium visibility
facInd <- unique(repMed$factor)
unifMed <- repMed
unifMed$prop[unifMed$levels != "NOT.SEEN"] <- unlist(tapply(repMed$prop[unifMed$levels != "NOT.SEEN"],
                                                              repMed$factor[unifMed$levels != "NOT.SEEN"],
                                                              function(x) rep(1/length(x), length(x)))[facInd])
for(i in 1:length(facInd)){
  if(any(unifMed$levels[unifMed$factor == facInd[i]] %in% c("NOT.SEEN"))){
    notProp <- unifMed$prop[unifMed$factor == facInd[i] & unifMed$levels == "NOT.SEEN"]
    levProp <- unifMed$prop[unifMed$factor == facInd[i] & unifMed$levels != "NOT.SEEN"]
    unifMed$prop[unifMed$factor == facInd[i] & unifMed$levels != "NOT.SEEN"] <- levProp*(1 - notProp)
  }
}

## Low visibility
facInd <- unique(repLow$factor)
unifLow <- repLow
unifLow$prop[unifLow$levels != "NOT.SEEN"] <- unlist(tapply(repLow$prop[unifLow$levels != "NOT.SEEN"],
                                                            repLow$factor[unifLow$levels != "NOT.SEEN"],
                                                            function(x) rep(1/length(x), length(x)))[facInd])
for(i in 1:length(facInd)){
  if(any(unifLow$levels[unifLow$factor == facInd[i]] %in% c("NOT.SEEN"))){
    notProp <- unifLow$prop[unifLow$factor == facInd[i] & unifLow$levels == "NOT.SEEN"]
    levProp <- unifLow$prop[unifLow$factor == facInd[i] & unifLow$levels != "NOT.SEEN"]
    unifLow$prop[unifLow$factor == facInd[i] & unifLow$levels != "NOT.SEEN"] <- levProp*(1 - notProp)
  }
}

## Create distributions that keep ratio of experimental distribution but alter visibility

## Low visibility

visProp <- 0.1
expLow <- expDist
facInd <- unique(expLow$factor)
expLow$prop[expLow$levels == "NOT.SEEN"] <- visProp


for(i in 1:length(facInd)){
  if(any(expLow$levels[expLow$factor == facInd[i]] %in% c("NOT.SEEN"))){
    notProp <- expLow$prop[expLow$factor == facInd[i] & expLow$levels == "NOT.SEEN"]
    levProp <- expLow$prop[expLow$factor == facInd[i] & expLow$levels != "NOT.SEEN"]
    levShare <- levProp/sum(levProp)
    expLow$prop[expLow$factor == facInd[i] & expLow$levels != "NOT.SEEN"] <- levShare * (1 - notProp)
  }
}


## High visibility

visProp <- 0.9
expHigh <- expDist
facInd <- unique(expHigh$factor)
expHigh$prop[expHigh$levels == "NOT.SEEN"] <- visProp


for(i in 1:length(facInd)){
  if(any(expHigh$levels[expHigh$factor == facInd[i]] %in% c("NOT.SEEN"))){
    notProp <- expHigh$prop[expHigh$factor == facInd[i] & expHigh$levels == "NOT.SEEN"]
    levProp <- expHigh$prop[expHigh$factor == facInd[i] & expHigh$levels != "NOT.SEEN"]
    levShare <- levProp/sum(levProp)
    expHigh$prop[expHigh$factor == facInd[i] & expHigh$levels != "NOT.SEEN"] <- levShare * (1 - notProp)
  }
}



## Build marginal distribution list for function

marDist <- prepDist(expDist, list(repHigh, repMed, repLow,
                                  demHigh, demMed, demLow,
                                  genHigh, genMed, genLow,
                                  expHigh, expLow,
                                  unifHigh, unifMed, unifLow, unifPure), name.ref = "Exp",
                    name.clean = c("Rep High", "Rep Med", "Rep Low",
                                   "Dem High", "Dem Med", "Dem Low",
                                   "Gen High", "Gen Med", "Gen Low",
                                   "Exp High", "Exp Low",
                                   "Unif High", "Unif Med", "Unif Low", "Unif Pure"))#[c(2, 1)] # Change the order so exp is first

##################################### ESTIMATE MODEL AND PLOT AMEs ###########################################

numBoot <- 100

############################# AMEs

###### DEMOCRATS

dfUse <- as.data.frame(na.omit(dfNew[dfNew$party.id == "Democrat",
                                     c(expVars, depVar, idVars)]))
dfUse$id <- as.numeric(factor(dfUse$id))

### Fix factor marginals at uniform, alter visibility
formula_u <- as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"]))
data <- dfUse
data0 <- data[order(data$pair_id),]
side <- rep(c(1,0), times=nrow(data0)/2)
data1 <- data0[side==1,]
data2 <- data0[side==0,]
X1 <- model.matrix(formula_u, data=data1)[ ,-1]
X2 <- model.matrix(formula_u, data=data2)[ ,-1]
X <- cbind(1, X1 - X2)
cluster <- rep(seq(1:(1540/2)), each = 2)

marUse <- marDist[c("Unif High", "Unif Med", "Unif Low")]
marUse <- marDist[c("Exp", "Unif Med")]
formual_u <- chosen.candidate ~

ameVisDem1 <- AME_estimate_full(formula = as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"])),
                                data = dfUse, pair = FALSE,
                                # pair = T, pair_id = dfUse$pair_id,
                                cluster = dfUse$id,
                                marginal_dist = marUse,
                                marginal_type = names(marUse), boot = numBoot,
                                cv.type = "cv.min", type = "No-Reg", difference = F)
ameVisDem1$AME

lapply(model.frame(as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"])),
       data = dfUse), levels)

ameVisDem2 <- AME_estimate_full(formula = as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"])),
                                data = dfUse,
                                pair = TRUE, pair_id = dfUse$pair_id,
                                cluster = dfUse$id,
                                marginal_dist = marUse,
                                marginal_type = names(marUse), boot = numBoot,
                                cv.type = "cv.min", type = "No-Reg", difference = F)
ameVisDem2$AME

data <- dfUse
data0 <- data[order(data$pair_id),]
side <- rep(c(1,0), times=nrow(data0)/2)
data1 <- data0[side==1,]
data2 <- data0[side==0,]
ameVisDem3 <- AME_estimate_full(formula = as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"])),
                                 data = data1,
                                 pair = FALSE,
                                 cluster = data1$id,
                                 marginal_dist = marUse,
                                 marginal_type = names(marUse), boot = numBoot,
                                 cv.type = "cv.min", type = "No-Reg", difference = F)


AME1 <- do.call("rbind", ameVisDem1$AME)
AME2 <- do.call("rbind", ameVisDem2$AME)
AME3 <- do.call("rbind", ameVisDem3$AME)

AME1[1:10, 3:5]
AME2[1:10, 3:5]
AME3[1:10, 3:5]

plot(AME1[,4], AME2[,4])
abline(0, 1)


ameVisDem2 <- AME_estimate_full(formula = as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"])),
                                data = dfUse,
                                pair = T, pair_id = dfUse$pair_id,
                                cluster = dfUse$id,
                                marginal_dist = marUse,
                                marginal_type = names(marUse), boot = numBoot,
                                cv.type = "cv.min", type = "No-Reg", difference = F)

ameVisDem1$AME$party.candidate
ameVisDem10$AME$party.candidate
ameVisDem2$AME$party.candidate

ameVisDem1$AME$education.candidate
ameVisDem10$AME$education.candidate
ameVisDem2$AME$education.candidate

plot(ameVisDem1$AME$education.candidate[,4], ameVisDem2$AME$education.candidate[,4])



formula_u <- as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"]))
factor_l <- length(all.vars(formula_u)[-1])
combMat <- combn(factor_l,2); intNames <- c()
for(k in 1:ncol(combMat)){
  intNames[k] <- paste(all.vars(formula_u)[-1][combMat[1,k]], "*", all.vars(formula_u)[-1][combMat[2,k]], sep = "")
}
formula_full <- as.formula(paste(all.vars(formula_u)[1], "~", paste(intNames, collapse = "+"), sep=""))

x1 <- model.matrix(formula_full, data = data1)
x2 <- model.matrix(formula_full, data = data2)
head(x1)

colnames(x1)[apply(x1 - x2, 2, function(x) all(x == 0))]

apply(x1[, apply(x1 - x2, 2, function(x) all(x == 0))], 2, mean)

ameVisDem1$AME$party.candidate
ameVisDem2$AME$party.candidate

ameVisDem1$AME$education.candidate
ameVisDem2$AME$education.candidate

plot(ameVisDem1$AME$education.candidate[,4], ameVisDem2$AME$education.candidate[,4])

ameVisDem <- AME_estimate_full(formula = as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"])),
                            data = dfUse,
                            pair = T, pair_id = dfUse$pair_id,
                            # cluster = dfUse$id,
                            marginal_dist = marUse,
                            marginal_type = names(marUse), boot = numBoot,
                            cv.type = "cv.min",
                            difference = F)

### Fix level ratios at experimental, alter visibility

marUse <- marDist[c("Exp Low", "Exp High")]
ameExpDem <- AME_estimate_full(formula = as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"])),
                               data = dfUse,
                               #pair = T, pair_id = dfUse$pair_id, # Note that the differencing is off here
                               #cluster = dfUse$id,
                               marginal_dist = marUse,
                               marginal_type = names(marUse), boot = numBoot,
                               cv.type = "cv.min",
                               difference = F)


## Fix visibility, change other marginals

marUse <- marDist[c("Rep Med", "Dem Med")]
ameMarDem <- AME_estimate_full(formula = as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"])),
                            data = dfUse,
                            #pair = T, pair_id = dfUse$pair_id, # Note that the differencing is off here
                            #cluster = dfUse$id,
                            marginal_dist = marUse,
                            marginal_type = names(marUse), boot = numBoot,
                            cv.type = "cv.min",
                            difference = F)


###### REPUBLICANS

dfUse <- as.data.frame(na.omit(dfNew[dfNew$party.id == "Republican",
                                     c(expVars, depVar, idVars)]))
dfUse$id <- as.numeric(factor(dfUse$id))

### Fix factor marginals at uniform, alter visibility

marUse <- marDist[c("Unif High", "Unif Med", "Unif Low")]
ameVisRep <- AME_estimate_full(formula = as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"])),
                               data = dfUse,
                               #pair = T, pair_id = dfUse$pair_id, # Note that the differencing is off here
                               #cluster = dfUse$id,
                               marginal_dist = marUse,
                               marginal_type = names(marUse), boot = numBoot,
                               cv.type = "cv.min",
                               difference = F)

### Fix level ratios at experimental, alter visibility

marUse <- marDist[c("Exp Low", "Exp High")]
ameExpRep <- AME_estimate_full(formula = as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"])),
                               data = dfUse,
                               #pair = T, pair_id = dfUse$pair_id, # Note that the differencing is off here
                               #cluster = dfUse$id,
                               marginal_dist = marUse,
                               marginal_type = names(marUse), boot = numBoot,
                               cv.type = "cv.min",
                               difference = F)

## Fix visibility, change other marginals

marUse <- marDist[c("Rep Med", "Dem Med")]
ameMarRep <- AME_estimate_full(formula = as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"])),
                               data = dfUse,
                               #pair = T, pair_id = dfUse$pair_id, # Note that the differencing is off here
                               #cluster = dfUse$id,
                               marginal_dist = marUse,
                               marginal_type = names(marUse), boot = numBoot,
                               cv.type = "cv.min",
                               difference = F)


###################### PLOT



######## VISIBILITY

facInc <- c("party.candidate")
#facInc <- names(ameUse$AME)
compLevs <- c("Black", "NOT.SEEN",
              "NOT.SEEN", "Republican") # This is for the CAMEs
par(mfrow = c(1, 2))

## Democrats

ameUse <- ameVisDem
sampUse <- "Democratic Subjects"
visUse <- ""

plot.AME(ameUse, factor_name = facInc, plot_difference = "no",
         col = c("black", "red", "blue", "darkgreen"),
         main = c("Democratic Subjects"),
         mar = 10)


## Republicans

ameUse <- ameVisRep
sampUse <- "Republican Subjects"
visUse <- ""

plot.AME(ameUse, factor_name = facInc, plot_difference = "no",
         col = c("black", "red", "blue", "darkgreen"),
         main = c("Republican Subjects"),
         mar = 10)


###### CANDIDATE MARGINALS

facInc <- c("spending.stance.candidate")

par(mfrow = c(2, 1))


## Democrats

ameUse <- ameMarDem
sampUse <- "Democratic Subjects"
visUse <- ""
plot.AME(ameUse, factor_name = facInc, plot_difference = "no",
         col = c("black", "red", "blue"),
         main = c("Democratic Subjects"),
         mar = 15, xlim = c(-0.2, 0.4))

## Republicans

ameUse <- ameMarRep
sampUse <- "Republican Subjects"
visUse <- ""
plot.AME(ameUse, factor_name = facInc, plot_difference = "no",
         col = c("black", "red", "blue"),
         main = c("Republican Subjects"),
         mar = 15, xlim = c(-0.2, 0.4))


############################## SELECTED DECOMPOSITIONS

###### REPUBLICANS

dfUse <- as.data.frame(na.omit(dfNew[dfNew$party.id == "Republican",
                                     c(expVars, depVar, idVars)]))
dfUse$id <- as.numeric(factor(dfUse$id))
marUse <- marDist[c("Exp", "Unif High", "Unif Low")]
ameDecompRep <- AME_estimate_full(formula = as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"])),
                               data = dfUse,
                               #pair = T, pair_id = dfUse$pair_id, # Note that the differencing is off here
                               #cluster = dfUse$id,
                               marginal_dist = marUse,
                               marginal_type = names(marUse), boot = numBoot,
                               cv.type = "cv.min",
                               difference = F)

###### DEMOCRATS

dfUse <- as.data.frame(na.omit(dfNew[dfNew$party.id == "Democrat",
                                     c(expVars, depVar, idVars)]))
dfUse$id <- as.numeric(factor(dfUse$id))
marUse <- marDist[c("Exp", "Unif High", "Unif Low")]
ameDecompDem <- AME_estimate_full(formula = as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"])),
                                  data = dfUse,
                                  #pair = T, pair_id = dfUse$pair_id, # Note that the differencing is off here
                                  #cluster = dfUse$id,
                                  marginal_dist = marUse,
                                  marginal_type = names(marUse), boot = numBoot,
                                  cv.type = "cv.min",
                                  difference = F)


########### COMPARISONS

## Change for high visibility

par(mfrow = c(1,2), oma = c(0.5, 0.5, 2, 0.5))
plot_decompose_alt(ameDecompRep, factor_name = "party.candidate",
                   level_name = "Republican", type = "decompose",
                   mar = 15, marginal_diff = c("Exp", "Unif High"),
                   main = "Republican Subjects")
plot_decompose_alt(ameDecompDem, factor_name = "party.candidate",
                   level_name = "Republican", type = "decompose",
                   mar = 15, marginal_diff = c("Exp", "Unif High"),
                   main = "Democratic Subjects")
mtext("Decomposition of AME Change: High Visibility", line = -1, cex = 1.8, outer = T)

## Change for low visibility

par(mfrow = c(1,2), oma = c(0.5, 0.5, 2, 0.5))
plot_decompose_alt(ameDecompRep, factor_name = "party.candidate",
                   level_name = "Republican", type = "decompose",
                   mar = 15, marginal_diff = c("Exp", "Unif Low"),
                   main = "Republican Subjects")
plot_decompose_alt(ameDecompDem, factor_name = "party.candidate",
                   level_name = "Republican", type = "decompose",
                   mar = 15, marginal_diff = c("Exp", "Unif Low"),
                   main = "Democratic Subjects")
mtext("Decomposition of AME Change: Low Visibility", line = -1, cex = 1.8, outer = T)


###################### SELECTED CAMES

plotFac <- c("abortion.stance.candidate", "spending.stance.candidate")

############## VISIBILITY

whichFac <- "party.candidate"
whichLev <- c("Republican")
ameUse <- ameDecompRep
cameUse <- cAME_from_boot(ameUse, factor_name = whichFac,
                          level_name = whichLev, difference = F)
#par(oma = c(0.5, 0.5, 5, 0.5))
plot_cAME_alt(cameUse,
              factor_name = plotFac[!plotFac %in% whichFac],
              main = c("Marginal Distributions", "Conditional AMEs of Candidate \n Being Republican"),
              col = c("black", "red", "blue"),
              marginal_effect = "Exp",
              legend_pos = "topright",
              plot_all = T,
              plot_difference = "no",
              cex = 1, mar = 22, topmar = 3)
mtext("Republican Subjects", outer = T, cex = 1.8, line = 1)


ameUse <- ameDecompDem
cameUse <- cAME_from_boot(ameUse, factor_name = whichFac,
                          level_name = whichLev, difference = F)
#par(oma = c(0.5, 0.5, 5, 0.5))
plot_cAME_alt(cameUse,
              factor_name = plotFac[!plotFac %in% whichFac],
              main = c("Marginal Distributions", "Conditional AMEs of Candidate \n Being Republican"),
              col = c("black", "red", "blue"),
              marginal_effect = "Exp",
              legend_pos = "topright",
              plot_all = T,
              plot_difference = "no",
              cex = 1, mar = 22, topmar = 3)
mtext("Democratic Subjects", outer = T, cex = 1.8, line = 1)


############################################## ADDITIONAL ANALYSIS ###############################################


## For each factor, set NOT.SEEN to 0, re-assign the weight according to experimental distribution to remaining levels
## then set remaining factors at a fixed visibility, compare this to the standard estimate

######################## DEMOCRATS

dfDem <- as.data.frame(na.omit(dfNew[dfNew$party.id == "Democrat",
                                     c(expVars, depVar, idVars)]))
dfDem$id <- as.numeric(factor(dfDem$id))
dfRep <- as.data.frame(na.omit(dfNew[dfNew$party.id == "Republican",
                                     c(expVars, depVar, idVars)]))
dfRep$id <- as.numeric(factor(dfRep$id))
ameListDem <- ameListRep <- list()


for(i in 1:length(facInd[-10])){ # Don't include party here because it's always seen (and this is the AME you're interested in anyway)
  distHigh <- expHigh
  distLow <- expLow
  distMed <- expDist
  notProp <- 0


  ## Change high visibility distribution for factor i
  distHigh$prop[distHigh$factor == facInd[i] & distHigh$levels == "NOT.SEEN"] <- notProp
  levProp <- distHigh$prop[distHigh$factor == facInd[i] & distHigh$levels != "NOT.SEEN"]
  levShare <- levProp/sum(levProp)
  distHigh$prop[distHigh$factor == facInd[i] & distHigh$levels != "NOT.SEEN"] <- levShare * (1 - notProp)

  ## Change medium visibility distribution for factor i
  distMed$prop[distMed$factor == facInd[i] & distMed$levels == "NOT.SEEN"] <- notProp
  levProp <- distMed$prop[distMed$factor == facInd[i] & distMed$levels != "NOT.SEEN"]
  levShare <- levProp/sum(levProp)
  distMed$prop[distMed$factor == facInd[i] & distMed$levels != "NOT.SEEN"] <- levShare * (1 - notProp)

  ## Change low visibility distribution for factor i
  distLow$prop[distLow$factor == facInd[i] & distLow$levels == "NOT.SEEN"] <- notProp
  levProp <- distLow$prop[distLow$factor == facInd[i] & distLow$levels != "NOT.SEEN"]
  levShare <- levProp/sum(levProp)
  distLow$prop[distLow$factor == facInd[i] & distLow$levels != "NOT.SEEN"] <- levShare * (1 - notProp)

  ## Estimate AMEs with new distributions

  marUse <- list(expDist, distHigh, distMed, distLow)
  names(marUse) <- c("Exp", "Exp High", "Exp Med", "Exp Low")

  ## Democrats
  ameListDem[[i]] <- AME_estimate_full(formula = as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"])),
                                    data = dfDem,
                                    #pair = T, pair_id = dfUse$pair_id, # Note that the differencing is off here
                                    cluster = dfDem$id,
                                    marginal_dist = marUse,
                                    marginal_type = names(marUse), boot = numBoot,
                                    cv.type = "cv.min",
                                    difference = F)$AME # Note just taking AME object here

  ## Republicans
  ameListRep[[i]] <- AME_estimate_full(formula = as.formula(make.equation(depVar, expVars[!expVars %in% "copartisan"])),
                                       data = dfRep,
                                       #pair = T, pair_id = dfUse$pair_id, # Note that the differencing is off here
                                       #cluster = dfUse$id,
                                       marginal_dist = marUse,
                                       marginal_type = names(marUse), boot = numBoot,
                                       cv.type = "cv.min",
                                       difference = F)$AME # Note just taking AME object here

}

names(ameListDem) <- names(ameListRep) <- facInd[-10]


### PLOT

############## REPUBLICAN SUBJECTS

par(mfrow = c(1,3))
## Low visibility

## Create AME Object
facUse <- facInd[-10]
ameTmp <- ameUse

for(i in 1:length(facUse)){
  ameTmp$AME[[i]] <- ameListRep[[facUse[i]]]$party.candidate[ameListRep[[facUse[i]]]$party.candidate$type == "Exp Low", ]
  ameTmp$AME[[i]]$level <- ""
}
names(ameTmp$AME) <- facUse

plot.AME(ameTmp, factor_name = facUse, plot_difference = "no",
         main = c("Effect of Republican Candidate \n Republican Subjects, Low Visibility"),
         mar  = 15)

## Medium visibility

## Create AME Object
facUse <- facInd[-10]
ameTmp <- ameUse

for(i in 1:length(facUse)){
  ameTmp$AME[[i]] <- ameListRep[[facUse[i]]]$party.candidate[ameListRep[[facUse[i]]]$party.candidate$type == "Exp Med", ]
  ameTmp$AME[[i]]$level <- ""
}
names(ameTmp$AME) <- facUse

plot.AME(ameTmp, factor_name = facUse, plot_difference = "no",
         main = c("Effect of Republican Candidate \n Republican Subjects, Medium Visibility"),
         mar  = 15)



## High visibility

## Create AME Object
facUse <- facInd[-10]
ameTmp <- ameUse

for(i in 1:length(facUse)){
  ameTmp$AME[[i]] <- ameListRep[[facUse[i]]]$party.candidate[ameListRep[[facUse[i]]]$party.candidate$type == "Exp High", ]
  ameTmp$AME[[i]]$level <- ""
}
names(ameTmp$AME) <- facUse

plot.AME(ameTmp, factor_name = facUse, plot_difference = "no",
         main = c("Effect of Republican Candidate \n Republican Subjects, High Visibility"),
         mar  = 15)


##################### DEMOCRATS

par(mfrow = c(1, 3))
## Low visibility

## Create AME Object
facUse <- facInd[-10]
ameTmp <- ameUse

for(i in 1:length(facUse)){
  ameTmp$AME[[i]] <- ameListDem[[facUse[i]]]$party.candidate[ameListDem[[facUse[i]]]$party.candidate$type == "Exp Low", ]
  ameTmp$AME[[i]]$level <- ""
}
names(ameTmp$AME) <- facUse

plot.AME(ameTmp, factor_name = facUse, plot_difference = "no",
         main = c("Effect of Republican Candidate \n Democratic Subjects, Low Visibility"),
         mar = 15)

## Medium visibility

## Create AME Object
facUse <- facInd[-10]
ameTmp <- ameUse

for(i in 1:length(facUse)){
  ameTmp$AME[[i]] <- ameListDem[[facUse[i]]]$party.candidate[ameListDem[[facUse[i]]]$party.candidate$type == "Exp Med", ]
  ameTmp$AME[[i]]$level <- ""
}
names(ameTmp$AME) <- facUse

plot.AME(ameTmp, factor_name = facUse, plot_difference = "no",
         main = c("Effect of Republican Candidate \n Democratic Subjects, Medium Visibility"),
         mar = 15)


## High visibility

## Create AME Object
facUse <- facInd[-10]
ameTmp <- ameUse

for(i in 1:length(facUse)){
  ameTmp$AME[[i]] <- ameListDem[[facUse[i]]]$party.candidate[ameListDem[[facUse[i]]]$party.candidate$type == "Exp High", ]
  ameTmp$AME[[i]]$level <- ""
}
names(ameTmp$AME) <- facUse

plot.AME(ameTmp, factor_name = facUse, plot_difference = "no",
         main = c("Effect of Republican Candidate \n Democratic Subjects, High Visibility"),
         mar = 15)



