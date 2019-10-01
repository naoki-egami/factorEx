#' Estimating the population AMCE using a design-based approach
#' @param formula Formula
#' @param data Data
#' @param factor_name Factors for which the function estimates the pAMCEs. If not specified, the function estimates for all factors.
#' @param pair Whether we use a paired-choice conjoint design.
#' @param pair_id Unique identifiers for pairs in the paired-choice conjoint design  (optional).
#' @param cluster_id Unique identifiers for computing cluster standard errors (optional).
#' @param cross_int Include interactions across profiles. Default is FALSE.
#' @param target_dist Target profile distributions to be used. See Examples in the GitHub page for details.
#' @param target_type Types of target profile distributions. `marginal`, 'partial_joint', or `target_data`.
#' @param partial_joint_name (when `target_type = "partial_joint"`) Names of factors representing partial joint distributions. See Examples in the GitHub page for details.
#' @import arm
#' @importFrom estimatr lm_robust
#' @importFrom sandwich sandwich estfun
#' @return \code{design_pAMCE} returns an object of \code{pAMCE} class.
#'  \itemize{
#'    \item \code{AMCE}: Estimates of the pAMCE for all factors.
#'    \item \code{design_weight}: Weight for each observation used in the weighted difference-in-means.
#'    \item \code{approach}: "design_based"
#'    \item \code{input}: Input into the function.
#'    \item \code{...}: Values for internal use.
#'  }
#' @examples
#' # Small example
#' data("OnoBurden")
#' OnoBurden_data_pr_s <- OnoBurden$OnoBurden_data_pr[1:500, ]
#' # randomization based on marginal population design
#' target_dist_marginal_s <- OnoBurden$target_dist_marginal[c("gender", "age")]
#'
#' # design-based estimation
#' out_design_mar_s <-
#'   design_pAMCE(formula = Y ~ gender + age,
#'                factor_name = "gender",
#'                data = OnoBurden_data_pr_s,
#'                pair_id = OnoBurden_data_pr_s$pair_id,
#'                cluster_id = OnoBurden_data_pr_s$id,
#'                target_dist  = target_dist_marginal_s, target_type = "marginal")
#' summary(out_design_mar_s)
#' \donttest{
#'   #  Example
#'   data("OnoBurden")
#'   OnoBurden_data_pr <- OnoBurden$OnoBurden_data_pr
#'   # randomization based on marginal population design
#'   target_dist_marginal <- OnoBurden$target_dist_marginal
#'
#'   # design-based estimation
#'   out_design_mar <-
#'     design_pAMCE(formula = Y ~ gender + age + family + race + experience + party + pos_security,
#'                  factor_name = c("gender", "age", "experience"),
#'                  data = OnoBurden_data_pr,
#'                  pair_id = OnoBurden_data_pr$pair_id,
#'                  cluster_id = OnoBurden_data_pr$id,
#'                  target_dist  = target_dist_marginal, target_type = "marginal")
#'  summary(out_design_mar)
#' }
#' @description \code{design_pAMCE} implements the design-based approach to estimate the pAMCE. See de la Cuesta, Egami, and Imai (2019+) for details. More examples are available at the GitHub page of \code{factorEx}.
#' @references de la Cuesta, Egami, and Imai. (2019+). Improving the External Validity of Conjoint Analysis: The Essential Role of Profile Distribution. (Working Paper). Available at \url{https://scholar.princeton.edu/sites/default/files/negami/files/conjoint_profile.pdf}.
#' @references Egami and Imai. (2019). Causal Interaction in Factorial Experiments: Application to Conjoint Analysis. Journal of the American Statistical Association, Vol.114, No.526 (June), pp. 529–540. Available at \url{https://scholar.princeton.edu/sites/default/files/negami/files/causalint.pdf}.
#' @export

design_pAMCE <- function(formula, factor_name,
                         data, pair = FALSE, pair_id = NULL,
                         cross_int = FALSE,
                         cluster_id = NULL,
                         target_dist, target_type,
                         partial_joint_name){

  ##################
  ## HouseKeeping ##
  ##################
  if(missing(pair_id) == TRUE) pair_id <- NULL
  if(missing(cluster_id) == TRUE) cluster_id <- seq(1:nrow(data))
  if(missing(target_type) == TRUE){ target_type <- "marginal" }

  if((target_type %in% c("marginal", "partial_joint", "target_data")) == FALSE){
    warning(" 'target_type' should be 'marginal', 'partial_joint', or 'target_data' ")
  }
  if((target_type %in% c("partial_joint")) == FALSE){
    partial_joint_name <-  NULL
  }

  if(any(is.na(data))==TRUE){
    stop("Remove NA before using the function.")
  }
  if(pair==TRUE & is.null(pair_id)==TRUE){
    stop("When 'pair=TRUE', specify 'pair_id'.")
  }
  if(pair==TRUE & all(table(pair_id)==2)==FALSE){
    stop("When 'pair=TRUE', each of 'pair_id' should have two observations.")
  }
  if(is.null(pair_id) == FALSE) pair <- TRUE

  if(pair == FALSE){
    cross_int <- FALSE
  }
  if(missing(cluster_id) == TRUE) cluster_id <- seq(1:nrow(data))

  ## Store Baselines
  baseline <- lapply(model.frame(formula, data = data)[,-1], FUN = function(x) levels(x)[1])

  ## Compute weights
  if(missing(factor_name)){ factor_name <-  all.vars(formula)[-1] }
  data$cluster_id <- cluster_id

  AMCE_list  <- list()
  design_weight_list <- list()

  message("Estimaing the pAMCEs for ", appendLF = FALSE)
  for(z in 1:length(factor_name)){
    factor_name_z <- factor_name[z]
    out_w <- weights_pAMCE(formula =  formula, factor_name = factor_name_z,
                           data = data, pair = pair, pair_id = pair_id, cross_int = cross_int,
                           target_dist =  target_dist,  target_type  =  target_type,
                           partial_joint_name = partial_joint_name)

    formula_w <- as.formula(paste(as.character(formula)[2], "~", factor_name_z, sep =  ""))
    message(paste(factor_name_z,  "...", sep = ""), appendLF = FALSE)
    design_weight <- out_w$design_weight
    suppressWarnings(lm_w <- lm_robust(formula_w,
                                       data  =  out_w$newdata,
                                       weights = design_weight, clusters = cluster_id))

    lm_no_w <- lm(formula_w, data  =  data)
    se_no_w <- sqrt(diag(cluster_se_glm(lm_no_w, cluster = as.factor(data$cluster_id)))[-1])
    coef_no_w <- cbind(lm_no_w$coefficients[-1], se_no_w, lm_no_w$coefficients[-1] - 1.96*se_no_w,
                       lm_no_w$coefficients[-1] + 1.96*se_no_w)
    colnames(coef_no_w)  <- c("estimate", "se", "low.95ci", "high.95ci")

    design_weight_list[[z]]  <- out_w$design_weight
    AMCE0 <- as.data.frame(summary(lm_w)$coef[, c("Estimate", "Std. Error", "CI Lower", "CI Upper")])
    colnames(AMCE0)  <- c("estimate", "se", "low.95ci", "high.95ci")
    AMCE0 <- AMCE0[-1, ]

    AMCE <- rbind(AMCE0, coef_no_w)

    AMCE$level <- rep(gsub(factor_name_z, "", rownames(AMCE0)), 2)
    AMCE$factor <-  rep(factor_name_z, times = nrow(AMCE))
    AMCE$type <-  c(rep("target", times = nrow(AMCE)/2), rep("sample AMCE", times = nrow(AMCE)/2))
    AMCE <- AMCE[, c("type", "factor", "level", "estimate", "se", "low.95ci", "high.95ci")]
    rownames(AMCE) <- NULL

    AMCE <- AMCE[order(AMCE$level),  ]
    AMCE_list[[z]] <- AMCE
  }
  names(AMCE_list) <- names(design_weight_list) <-  factor_name

  # out_w <- weights_pAMCE(formula =  formula, factor_name = factor_name,
  #                        data = data, pair = pair, pair_id = pair_id, cross_int = cross_int,
  #                        target_dist =  target_dist,  target_type  =  target_type,
  #                        partial_joint_name = partial_joint_name)
  #
  # formula_w <- as.formula(paste(as.character(formula)[2], "~", factor_name, sep =  ""))
  #
  # cat("Estimaing the pAMCEs...")
  # suppressWarnings(lm_w <- lm_robust(formula_w, data  =  out_w$newdata, weights = design_weight, clusters = cluster_id))
  #
  # lm_no_w <- lm(formula_w, data  =  data)
  # se_no_w <- sqrt(diag(cluster_se_glm(lm_no_w, cluster = as.factor(data$cluster_id)))[-1])
  # coef_no_w <- cbind(lm_no_w$coefficients[-1], se_no_w, lm_no_w$coefficients[-1] - 1.96*se_no_w,
  #                    lm_no_w$coefficients[-1] + 1.96*se_no_w)
  # colnames(coef_no_w)  <- c("estimate", "se", "low.95ci", "high.95ci")
  #
  # design_weight_list  <- out_w$design_weight
  # AMCE0 <- as.data.frame(summary(lm_w)$coef[, c("Estimate", "Std. Error", "CI Lower", "CI Upper")])
  # colnames(AMCE0)  <- c("estimate", "se", "low.95ci", "high.95ci")
  # AMCE0 <- AMCE0[-1, ]
  #
  # AMCE <- rbind(AMCE0, coef_no_w)
  #
  # AMCE$level <- rep(gsub(factor_name, "", rownames(AMCE0)), 2)
  # AMCE$factor <-  rep(factor_name, times = nrow(AMCE))
  # AMCE$type <-  c(rep("target", times = nrow(AMCE)/2), rep("sample AMCE", times = nrow(AMCE)/2))
  # AMCE <- AMCE[, c("type", "factor", "level", "estimate", "se", "low.95ci", "high.95ci")]
  # rownames(AMCE) <- NULL
  #
  # AMCE <- AMCE[order(AMCE$level),  ]
  #
  # AMCE_list  <- AMCE

  # input
  input <- list("formula" =  formula, "data" =  data,
                "target_dist" =  target_dist, "target_dist_name" = "target")
  out <- list("AMCE" = AMCE_list, "design_weight" = design_weight_list, "baseline" = baseline,
              "approach" = "design_based")
  out$input <- input

  class(out)  <- c(class(out), "pAMCE")
  return(out)

}
