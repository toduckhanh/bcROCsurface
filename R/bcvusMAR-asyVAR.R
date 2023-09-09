####========================================================================####
## This file consists of some functions that are related to the computation   ##
## asymptotic variance of VUS                                                 ##
##                																														##
####========================================================================####
##
#' @title Asymptotic variance estimation for VUS
#'
#' @description \code{asy_var_vus} computes the asymptotic variance of full data (FULL) and bias-corrected estimators (i.e. full imputation, mean score imputation, inverse probability weighting, semiparametric efficient and K nearest neighbor) of VUS.
#'
#' @param obj_vus  a result of a call to \code{\link{vus_mar}}.
#' @param diag_test  a numeric vector containing the diagnostic test values. \code{NA} values of \code{diag_test} are not accepted.
#' @param dise_vec  a n * 3  binary matrix with three columns, corresponding to the three classes of the disease status. In row i, 1 in column j indicates that the i-th subject belongs to class j, with j = 1, 2, 3. A row of \code{NA} values indicates a non-verified subject.
#' @param veri_stat  a binary vector containing the verification status (1 verified, 0 not verified).
#' @param rho_est  a result of a call to \code{\link{rho_mlogit}} of \code{\link{rho_knn}} to fit the disease model.
#' @param pi_est  a result of a call to \code{\link{psglm}} to fit the verification model.
#' @param boot a logical value. Default = \code{FALSE}. If set to \code{TRUE}, a bootstrap resampling is employed to estimate the asymptotic variance of the bias-corrected VUS estimators.
#' @param n_boot  the number of bootstrap replicates, which is used for FULL or KNN estimators, or option \code{boot = TRUE}. The defaut is 250.
#' @param parallel  a logical value. If \code{TRUE}, a parallel computing is employed in the bootstrap resampling process.
#' @param ncpus  number of processes to be used in parallel computing. Default is half of available cores.
#'
#'
#' @details
#' For the FULL estimator, a bootstrap resampling process or Jackknife approach is used to estimate the asymptotic variance, whereas, a bootstrap resampling process is employed to obtain the asymptotic variance of K nearest neighbor estimator.
#'
#' For the full imputation, mean score imputation, inverse probability weighting and semiparametric efficient estimators of VUS, the asymptotic variances are computed by using the explicit form. Furthermore, a bootstrap procedure is also available, useful in case of small sample sizes.
#'
#' @return \code{asy_var_vus} returns a estimated value of the asymptotic variance.
#'
#' @references
#' To Duc, K., Chiogna, M. and Adimari, G. (2020)
#' Nonparametric estimation of ROC surfaces in presence of verification bias.
#' \emph{REVSTAT-Statistical Journal}. \bold{18}, 5, 697–720.
#'
#' To Duc, K., Chiogna, M. and Adimari, G. (2016)
#' Bias-corrected methods for estimating the receiver operating characteristic surface of continuous diagnostic tests.
#' \emph{Electronic Journal of Statistics}, \bold{10}, 3063-3113.
#'
#' Guangming, P., Xiping, W. and Wang, Z. (2013)
#' Non-parameteric statistical inference for $P(X < Y < Z)$.
#' \emph{Sankhya A}, \bold{75}, 1, 118-138.
#'
#'
#' @examples
#' data(EOC)
#'
#' # Preparing the missing disease status
#' dise_na <- pre_data(EOC$D, EOC$CA125)
#' dise_vec_na <- dise_na$dise_vec
#' dise_fact_na <- dise_na$dise
#'
#' rho_out <- rho_mlogit(dise_fact_na ~ CA125 + CA153 + Age, data = EOC,
#'                       test = TRUE)
#' vus_fi <- vus_mar("fi", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'                   veri_stat = EOC$V, rho_est = rho_out, ci = FALSE)
#' var_fi <- asy_var_vus(vus_fi, diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'                       veri_stat = EOC$V, rho_est = rho_out)
#'
#'
#'
#' @importFrom Rcpp evalCpp
#' @import parallel
#' @import boot
#' @export
asy_var_vus <- function(obj_vus, diag_test, dise_vec, veri_stat = NULL,
                        rho_est = NULL, pi_est = NULL, boot = FALSE,
                        n_boot = 250, parallel = FALSE,
                        ncpus = ifelse(parallel, detectCores() / 2, NULL)) {
  if (!inherits(obj_vus, "vus_mar"))
    stop("The argument \"obj_vus\" is not a result of vus_mar()")
  ## checking the argument diag_test
  if (missing(diag_test)) stop("argument \"diag_test\" is missing \n")
  if (!inherits(diag_test, "numeric") || any(is.na(diag_test)))
    stop("\"diag_test\" must be a numeric vector and not include NA values")
  ## checking dise_vec
  if (missing(dise_vec)) stop("argument \"dise_vec\" is missing \n")
  if (!inherits(dise_vec, "matrix") || ncol(dise_vec) != 3 ||
        !all(is.element(na.omit(dise_vec), c(0, 1))))
    stop("variable \"dise_vec\" must be a binary matrix with 3 columns")
  if (length(diag_test) != nrow(dise_vec))
    stop(gettextf("arguments imply differing number of observation: %d",
                  length(diag_test)), gettextf(", %d", nrow(dise_vec)),
         domain = NA)
  vus_obj <- obj_vus$vus_fit
  meth <- tolower(attr(vus_obj, "name"))
  if (!is.element(meth, c("full", "fi")) && is.null(veri_stat))
    stop("argument \"veri_stat\" is missing \n")
  if (is.element(meth, c("fi", "msi", "knn", "spe")) && is.null(rho_est))
    stop("argument \"rho_est\" is missing \n")
  if (is.element(meth, c("ipw", "spe")) && is.null(pi_est))
    stop("argument \"pi_est\" is missing \n")
  ##
  if (!is.null(rho_est)) {
    if (!is.element(class(rho_est), c("prob_dise", "prob_dise_knn")))
      stop("\"rho_est\" not a result of rhoMlogit or rho_knn \n")
    if (is.element(class(rho_est), c("prob_dise_knn"))) {
      if (is.null(rho_est$K))
        stop("Please, choose one value of K\n")
    }
  }
  if (!is.null(pi_est)) {
    if (!is.element(class(pi_est), c("prob_veri")))
      stop("\"pi_est\" not a result of psglm \n")
  }
  n <- length(diag_test)
  dise_vec_temp <- dise_vec
  if (any(is.na(dise_vec))) dise_vec_temp[is.na(dise_vec)] <- 99
  if (meth == "fi") {
    if (!boot) {
      score <- mlog_est_func(dise_vec_temp, veri_stat, rho_est$X,
                             rho_est$values)
      hess <- solve(rho_est$Hess)
      der_rho1 <- t(rho_deriv(rho_est$X, rho_est$values, ref_level = "1"))
      der_rho2 <- t(rho_deriv(rho_est$X, rho_est$values, ref_level = "2"))
      der_rho3 <- - (der_rho1 + der_rho2)
      q_fi <- asy_var_vus_c(diag_test, rho_est$values, vus_obj, score, hess,
                            der_rho1, der_rho2, der_rho3)
      ans_var <- n^4 * sum(q_fi^2) / prod(colSums(rho_est$values)^2)
    } else {
      bst_fi <- function(dt, inds, formula) {
        dat <- dt[inds, ]
        out <- rho_mlogit(formula, data = dat)
        vus_c(dat[, 1], out$values)
      }
      name_var <- as.character(attributes(terms(rho_est$formula))$variables)[-1]
      data <- data.frame(diag_test, rho_est$D, rho_est$X[, name_var[-1]])
      names(data) <- c("Temp", name_var)
      if (parallel) {
        res_bst <- boot(data, bst_fi, R = n_boot, formula = rho_est$formula,
                        parallel = "snow", ncpus = ncpus)
      } else {
        res_bst <- boot(data, bst_fi, R = n_boot, formula = rho_est$formula)
      }
      ans_var <- as.numeric(var(res_bst$t, na.rm = TRUE))
    }
  } else if (meth == "msi") {
    if (!boot) {
      d_msi <- as.matrix(veri_stat * dise_vec_temp +
                           (1 - veri_stat) * rho_est$values)
      score <- mlog_est_func(dise_vec_temp, veri_stat, rho_est$X,
                             rho_est$values)
      hess <- solve(rho_est$Hess)
      der_d_msi1 <- t((1 - veri_stat) * rho_deriv(rho_est$X, rho_est$values,
                                          ref_level = "1"))
      der_d_msi2 <- t((1 - veri_stat) * rho_deriv(rho_est$X, rho_est$values,
                                          ref_level = "2"))
      der_d_msi3 <- - (der_d_msi1 + der_d_msi2)
      q_msi <- asy_var_vus_c(diag_test, d_msi, vus_obj, score, hess, der_d_msi1,
                           der_d_msi2, der_d_msi3)
      ans_var <- n^4 * sum(q_msi^2) / prod(colSums(d_msi)^2)
    } else {
      bst_msi <- function(dt, inds, formula) {
        dat <- dt[inds, ]
        out <- rho_mlogit(formula, data = dat)
        d_msi <- as.matrix(dat[, c("D1", "D2", "D3")] * dat[, "V"] +
                            (1 - dat[, "V"]) * out$values)
        vus_c(dat[, 1], d_msi)
      }
      name_var <- as.character(attributes(terms(rho_est$formula))$variables)[-1]
      data <- data.frame(diag_test, rho_est$D, rho_est$X[, name_var[-1]],
                         veri_stat, dise_vec_temp)
      names(data) <- c("Temp", name_var, "V", "D1", "D2", "D3")
      if (parallel) {
        res_bst <- boot(data, bst_msi, R = n_boot, formula = rho_est$formula,
                        parallel = "snow", ncpus = ncpus)
      } else {
        res_bst <- boot(data, bst_msi, R = n_boot, formula = rho_est$formula)
      }
      ans_var <- as.numeric(var(res_bst$t, na.rm = TRUE))
    }
  } else if (meth == "ipw") {
    if (!boot) {
      d_ipw <- as.matrix(veri_stat * dise_vec_temp / pi_est$values)
      score <- pi_est_func(veri_stat, pi_est$X, pi_est$values, pi_est$coeff,
                          pi_est$model)
      hess <- solve(pi_est$Hess)
      der_pi_inv <- pi_inv_deriv(pi_est$X, pi_est$values, pi_est$coeff,
                                 pi_est$model)
      der_d_ipw1 <- t(veri_stat * dise_vec_temp[, 1] * der_pi_inv)
      der_d_ipw2 <- t(veri_stat * dise_vec_temp[, 2] * der_pi_inv)
      der_d_ipw3 <- t(veri_stat * dise_vec_temp[, 3] * der_pi_inv)
      q_ipw <- asy_var_vus_c(diag_test, d_ipw, vus_obj, score, hess, der_d_ipw1,
                           der_d_ipw2, der_d_ipw3)
      ans_var <- sum(q_ipw^2) /
        prod((colSums(d_ipw) / sum(veri_stat / pi_est$values))^2) / n^2
    } else {
      bst_ipw <- function(dt, inds, formula) {
        dat <- dt[inds, ]
        out <- psglm(formula, data = dat, test = FALSE, trace = FALSE)
        d_ipw <- as.matrix(dat[, c("D1", "D2", "D3")] * dat[, 2] / out$values)
        vus_c(dat[, 1], d_ipw)
      }
      name_var <- as.character(attributes(terms(pi_est$formula))$variables)[-1]
      data <- data.frame(diag_test, veri_stat, pi_est$X[, name_var[-1]],
                         dise_vec_temp)
      names(data) <- c("Temp", name_var, "D1", "D2", "D3")
      if (parallel) {
        res_bst <- boot(data, bst_ipw, R = n_boot, formula = pi_est$formula,
                        parallel = "snow", ncpus = ncpus)
      } else {
        res_bst <- boot(data, bst_ipw, R = n_boot, formula = pi_est$formula)
      }
      ans_var <- as.numeric(var(res_bst$t, na.rm = TRUE))
    }
  } else if (meth == "spe") {
    if (!boot) {
      temp1 <- veri_stat / pi_est$values
      d_spe <- as.matrix(temp1 * dise_vec_temp - (temp1 - 1) * rho_est$values)
      score1 <- mlog_est_func(dise_vec_temp, veri_stat, rho_est$X,
                              rho_est$values)
      score2 <- pi_est_func(veri_stat, pi_est$X, pi_est$values, pi_est$coeff,
                            pi_est$model)
      score <- cbind(score1, score2)
      hess1 <- cbind(solve(rho_est$Hess),
                     matrix(0, nrow = nrow(rho_est$Hess),
                            ncol = ncol(pi_est$Hess)))
      hess2 <- cbind(matrix(0, nrow = nrow(pi_est$Hess),
                            ncol = ncol(rho_est$Hess)),
                     solve(pi_est$Hess))
      hess <- rbind(hess1, hess2)
      der_rho1 <- rho_deriv(rho_est$X, rho_est$values, ref_level = "1")
      der_rho2 <- rho_deriv(rho_est$X, rho_est$values, ref_level = "2")
      der_pi_inv <- pi_inv_deriv(pi_est$X, pi_est$values, pi_est$coeff,
                                 pi_est$model)
      der_d_spe1 <- t(
        cbind(-(temp1 - 1) * der_rho1,
              veri_stat * (dise_vec_temp[, 1] - rho_est$values[, 1]) *
                der_pi_inv))
      der_d_spe2 <- t(
        cbind(-(temp1 - 1) * der_rho2,
              veri_stat * (dise_vec_temp[, 2] - rho_est$values[, 2]) *
                der_pi_inv))
      der_d_spe3 <- -(der_d_spe1 + der_d_spe2)
      q_spe <- asy_var_vus_c(diag_test, d_spe, vus_obj, score, hess, der_d_spe1,
                           der_d_spe2, der_d_spe3)
      ans_var <- n^4 * sum(q_spe^2) / prod(colSums(d_spe)^2)
    } else {
      bst_spe <- function(dt, inds, formula_rho, formula_pi) {
        dat <- dt[inds, ]
        out_rho <- rho_mlogit(formula_rho, data = dat)
        out_pi <- psglm(formula_pi, data = dat, test = FALSE, trace = FALSE)
        d_spe <- as.matrix(dat[, c("D1", "D2", "D3")] * dat[, "V"] /
                            out_pi$values - (dat[, "V"] / out_pi$values - 1) *
                            out_rho$values)
        vus_c(dat[, 1], d_spe)
      }
      name_var_rho <- as.character(
        attributes(terms(rho_est$formula))$variables)[-1]
      name_var_pi <- as.character(
        attributes(terms(pi_est$formula))$variables)[-1]
      data_rho <- data.frame(rho_est$X[, name_var_rho[-1]])
      data_pi <- data.frame(pi_est$X[, name_var_pi[-1]])
      data <- data.frame(diag_test, veri_stat, rho_est$D,
                         union(data_rho, data_pi), dise_vec_temp)
      names(data) <- c("Temp", name_var_pi[1], name_var_rho[1],
                       union(name_var_rho[-1], name_var_pi[-1]), "D1", "D2",
                       "D3")
      if (parallel) {
        res_bst <- boot(data, bst_spe, R = n_boot,
                        formula_rho = rho_est$formula,
                        formula_pi = pi_est$formula, parallel = "snow",
                        ncpus = ncpus)
      } else {
        res_bst <- boot(data, bst_spe, R = n_boot,
                        formula_rho = rho_est$formula,
                        formula_pi = pi_est$formula)
      }
      ans_var <- as.numeric(var(res_bst$t, na.rm = TRUE))
    }
  } else if (meth == "knn") {
    bst_knn <- function(dt, inds, k, type) {
      dat <- dt[inds, ]
      x_mat <- as.matrix(dat[, -c(1:5)])
      dise_bt <- as.matrix(dat[, c(2:4)])
      dise_bt[dise_bt == 99] <- NA
      rho_knn <- rho_knn(x_mat, dise_bt, dat[, 5], k = k, type = type)
      d_knn <- as.matrix(dat[, c(2:4)] * dat[, 5] +
                           (1 - dat[, 5]) * rho_knn$values)
      vus_c(dat[, 1], d_knn)
    }
    data <- data.frame(diag_test, dise_vec_temp, veri_stat, rho_est$X)
    if (parallel) {
      res_bst <- boot(data, bst_knn, R = n_boot, k = rho_est$K,
                      type = rho_est$type, parallel = "snow", ncpus = ncpus)
    } else {
      res_bst <- boot(data, bst_knn, R = n_boot, k = rho_est$K,
                      type = rho_est$type)
    }
    ans_var <- as.numeric(var(res_bst$t, na.rm = TRUE))
  } else {
    if (!boot) {
      diag_test1 <- diag_test[dise_vec[, 1] == 1]
      diag_test2 <- diag_test[dise_vec[, 2] == 1]
      diag_test3 <- diag_test[dise_vec[, 3] == 1]
      n1 <- length(diag_test1)
      n2 <- length(diag_test2)
      n3 <- length(diag_test3)
      vus_est_mult <- vus_obj * n1 * n2 * n3
      var_term <- vus_c_full_core(diag_test1, diag_test2, diag_test3)
      var_term_i <- sapply(1:n1, function(x) sum(var_term$ind1[-x]))
      var_term_j <- sapply(1:n2, function(x) sum(var_term$ind2[-x]))
      var_term_k <- sapply(1:n3, function(x) sum(var_term$ind3[-x]))
      ans_var <- var((vus_est_mult - var_term_i) / (n2 * n3)) / n1 +
        var((vus_est_mult - var_term_j) / (n1 * n3)) / n2 +
        var((vus_est_mult - var_term_k) / (n1 * n2)) / n3
    } else {
      bst_full <- function(dt, inds) {
        dat <- dt[inds, ]
        vus_c(dat[, 1], as.matrix(dat[, c(2:4)]))
      }
      data <- data.frame(diag_test, dise_vec_temp)
      if (parallel) {
        res_bst <- boot(data, bst_full, R = n_boot, parallel = "snow",
                        ncpus = ncpus)
      } else {
        res_bst <- boot(data, bst_full, R = n_boot)
      }
      ans_var <- as.numeric(var(res_bst$t, na.rm = TRUE))
    }
  }
  return(ans_var)
}
