a <- model.frame(formula_full, data=data1)
dim(attr(terms(a), "factors"))
terms(a)
attr(terms(a), "variables")

pair <- TRUE

pair_var <- "fav_rating"

pair_var <- c("fav_rating", "pos_deficit")

terms(a)

terms(Xa)
Xa <- model.matrix(a)
attr(terms(a), "term.labels")
attr(Xa, "assign")


Xf <- model.frame(formula_full, data=data1)
Xf_o <- as.numeric(attr(terms(Xf), "order") == 1) + 1
Xf <- attr(terms(Xf), "factors")
Xf_ind <- which(apply(t(t(matrix(Xf[pair_var, ], ncol = ncol(Xf))*Xf_o)), 2, sum)  == 2)
X01 <- model.matrix(formula_full, data=data1)
X02 <- model.matrix(formula_full, data=data2)
pair_var_ind <- is.element(attr(X01, "assign"), Xf_ind)
X0 <- X01 - X02
X0[, pair_var_ind == TRUE] <- X01[, pair_var_ind == TRUE]
X <- cbind(1, X0[,-1])
ind_b <- attr(X01, "assign")[-1]
