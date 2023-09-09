####========================================================================####
## The R code for the asymptotic covariance matrix of true class fractions.   ##
##                                                                            ##
##                                                                            ##
####========================================================================####
#'
#' @title Asymptotic variance-covariance estimation for True Class Fractions (TCFs) at the cut point \eqn{(c_1, c_2)}
#'
#' @description
#' \code{asy_cov_tcf} computes the asymptotic variance-covariance matrix of full data (FULL) and bias-corrected estimators (i.e. full imputation, mean score imputation, inverse probability weighting, semiparametric efficient and K nearest neighbor) of TCFs.
#'
#'
#' @param obj_tcf  a result of a call to \code{\link{rocs.tcf}}.
#' @param diag_test  a numeric vector containing the diagnostic test values. \code{NA} values of \code{diag_test} are not accepted.
#' @param dise_vec a n * 3  binary matrix with three columns, corresponding to the three classes of the disease status. In row i, 1 in column j indicates that the i-th subject belongs to class j, with j = 1, 2, 3. A row of \code{NA} values indicates a non-verified subject.
#' @param veri_stat  a binary vector containing the verification status (1 verified, 0 not verified).
#' @param rho_est  a result of a call to \code{\link{rho_mlogit}} of \code{\link{rho_knn}} to fit the disease model.
#' @param pi_est  a result of a call to \code{\link{psglm}} to fit the verification model.
#' @param boot  a logical value. Default = \code{FALSE}. If set to \code{TRUE}, a bootstrap resampling is employed to estimate the asymptotic variance-covariance matrix of bias-corrected TCFs.
#' @param n_boot  the number of bootstrap replicates, used when \code{boot = TRUE} or for FULL estimator. Usually this will be a single positive integer. Default 250.
#' @param parallel  a logical value. If \code{TRUE}, a parallel computing is employed in the bootstrap resampling process.
#' @param ncpus  number of processes to be used in parallel computing. Default is half of available cores.
#'
#' @details For bias-corrected estimators of TCFs, the asymptotic variance-covariance matrix at a fixed cut point is estimated by using the Delta method. The function \code{asy_cov_tcf} implements the explicit forms presented in To Duc et al. (2016, 2020). In addition, the bootstrap procedure is also available.
#'
#' For FULL estimator, the asymptotic variance-covariance matrix is computed via bootstrap only.
#'
#' @return This function returns an estimated asymptotic variance-covariance matrix for FULL estimator and bias-corrected estimators of TCFs at a fixed cut point.
#'
#' @references
#'
#' To Duc, K., Chiogna, M. and Adimari, G. (2016)
#' Bias-corrected methods for estimating the receiver operating characteristic surface of continuous diagnostic tests.
#' \emph{Electronic Journal of Statistics}, \bold{10}, 3063-3113.
#'
#' To Duc, K., Chiogna, M. and Adimari, G. (2020)
#' Nonparametric estimation of ROC surfaces in presence of verification bias.
#' \emph{REVSTAT-Statistical Journal}. \bold{18}, 5, 697–720.
#'
#'
#' @examples
#' data(EOC)
#'
#' # FULL data estimator
#' dise_full <- pre_data(EOC$D.full, EOC$CA125)
#' dise_vec_full <- dise_full$dise_vec
#'
#' full_tcf <- rocs.tcf("full", diag_test = EOC$CA125, dise_vec = dise_vec_full,
#'                      cps = c(2, 4))
#' full_var <- asy_cov_tcf(full_tcf, diag_test = EOC$CA125,
#'                         dise_vec = dise_vec_full)
#'
#' # Preparing the missing disease status
#' dise_na <- pre_data(EOC$D, EOC$CA125)
#' dise_vec_na <- dise_na$dise_vec
#' dise_fact_na <- dise_na$dise
#'
#' rho_out <- rho_mlogit(dise_fact_na ~ CA125 + CA153 + Age, data = EOC,
#'                       test = TRUE)
#'
#' ## FI estimator
#' fi_tcf <- rocs.tcf("fi", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'                    veri_stat = EOC$V, rho_est = rho_out, cps = c(2, 4))
#' fi_var <- asy_cov_tcf(fi_tcf, diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'                       veri_stat = EOC$V, rho_est = rho_out)
#'
#' ## MSI estimator
#' msi_tcf <- rocs.tcf("msi", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'                     veri_stat = EOC$V, rho_est = rho_out, cps = c(2, 4))
#' msi_var <- asy_cov_tcf(msi_tcf, diag_test = EOC$CA125,
#'                        dise_vec = dise_vec_na, veri_stat = EOC$V,
#'                        rho_est = rho_out)
#'
#' ## IPW estimator
#' pi_out <- psglm(V ~ CA125 + CA153 + Age, data = EOC, test = TRUE)
#'
#' ipw_tcf <- rocs.tcf("ipw", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'                     veri_stat = EOC$V, pi_est = pi_out, cps = c(2, 4))
#' ipw_var <- asy_cov_tcf(ipw_tcf, diag_test = EOC$CA125,
#'                        dise_vec = dise_vec_na, veri_stat = EOC$V,
#'                        pi_est = pi_out)
#'
#' ## SPE estimator
#' spe_tcf <- rocs.tcf("spe", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'                     veri_stat = EOC$V, rho_est = rho_out, pi_est = pi_out,
#'                     cps = c(2, 4))
#' spe_var <- asy_cov_tcf(spe_tcf, diag_test = EOC$CA125,
#'                        dise_vec = dise_vec_na, veri_stat = EOC$V,
#'                        rho_est = rho_out, pi_est = pi_out)
#'
#' ## KNN estimators
#' x_mat <- cbind(EOC$CA125, EOC$CA153, EOC$Age)
#' rho_1nn <- rho_knn(x_mat = x_mat, dise_vec = dise_vec_na, veri_stat = EOC$V,
#'                    k = 1, type = "mahala")
#' knn_tcf <- rocs.tcf("knn", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'                     veri_stat = EOC$V, rho_est = rho_1nn, cps = c(2, 4))
#' knn_var <- asy_cov_tcf(knn_tcf, diag_test = EOC$CA125,
#'                        dise_vec = dise_vec_na, veri_stat = EOC$V,
#'                        rho_est = rho_1nn)
#'
#'
#' @import parallel
#' @export
asy_cov_tcf <- function(obj_tcf, diag_test, dise_vec, veri_stat = NULL,
                        rho_est = NULL, pi_est = NULL, boot = FALSE,
                        n_boot = 250, parallel = FALSE,
                        ncpus = ifelse(parallel, detectCores() / 2, NULL)) {
  if (is.list(obj_tcf) && length(obj_tcf) > 1) {
    check_tcf <- function(x) {
      class(x) == "tcfs"
    }
    if (all(sapply(obj_tcf, check_tcf)))
      stop("The \"obj_tcf\" is a list containing results of TCFs at more than 1 cut point! \n Please, choose one result of \"obj_tcf\" \n")
    else stop("The argument \"obj_tcf\" is not a result of roc.tcf()")
  }
  if (!inherits(obj_tcf, "tcfs"))
    stop("The argument \"obj_tcf\" is not a result of a call to roc.tcf()")
  ## checking the argument diag_test
  if (missing(diag_test)) stop("argument \"diag_test\" is missing \n")
  if (!inherits(diag_test, "numeric") || any(is.na(diag_test)))
    stop("\"diag_test\" must be a numeric vector and not include NA values")
  ## checking dise_vec
  if (missing(dise_vec)) stop("argument \"dise_vec\" is missing \n")
  if (!inherits(dise_vec, "matrix") || ncol(dise_vec) != 3 ||
      !all(is.element(na.omit(dise_vec), c(0, 1))))
    stop("\"dise_vec\" must be a binary matrix with 3 columns")
  if (length(diag_test) != nrow(dise_vec))
    stop(gettextf("arguments imply differing number of observation: %d",
                  length(diag_test)), gettextf(", %d", nrow(dise_vec)),
         domain = NA)
  ##
  method <- tolower(attr(obj_tcf, "name"))
  if (!is.element(method, c("full", "fi")) && is.null(veri_stat))
    stop("argument \"veri_stat\" is missing \n")
  if (is.element(method, c("fi", "msi", "knn", "spe")) && is.null(rho_est))
    stop("argument \"rho_est\" is missing \n")
  if (is.element(method, c("ipw", "spe")) && is.null(pi_est))
    stop("argument \"pi_est\" is missing \n")
  ##
  if (!is.null(rho_est)) {
    if (!is.element(class(rho_est), c("prob_dise", "prob_dise_knn")))
      stop("\"rho_est\" not a result of rho_mlogit or rho_knn \n")
    if (is.element(class(rho_est),  c("prob_dise_knn"))) {
      if (is.null(rho_est$K))
        stop("Please, choose one value of K \n")
    }
  }
  if (!is.null(pi_est)) {
    if (!is.element(class(pi_est), c("prob_veri")))
      stop("\"pi_est\" not a result of psglm \n")
  }
  ##
  dise_vec_temp <- dise_vec
  if (any(is.na(dise_vec))) dise_vec_temp[is.na(dise_vec)] <- 99
  tcf_orgi <- obj_tcf
  tcf_thet <- attr(tcf_orgi, "theta")
  tcf_bet <- attr(tcf_orgi, "beta")
  cp <- attr(tcf_orgi, "cp")
  if (method == "full") {
    bst_full <- function(dt, inds, cp) {
      dat <- dt[inds, ]
      rocs.tcf(method = "full", diag_test = dat[, 1],
               dise_vec = as.matrix(dat[, c(2:4)]), cps = cp)
    }
    data <- data.frame(diag_test, dise_vec)
    if (parallel) {
      res_bst <- boot(data, bst_full, R = n_boot, cp = cp, parallel = "snow",
                      ncpus = ncpus)
    } else {
      res_bst <- boot(data, bst_full, R = n_boot, cp = cp)
    }
    ans <- var(res_bst$t)
  } else if (method == "knn") {
    if (!boot) {
      x_mat <- rho_est$X
      rho_knn_est <- rho_knn(x_mat, dise_vec, veri_stat, k = 2,
                             type = rho_est$type)$values
      pi_knn_est <- psknn(x_mat, veri_stat, rho_est$type)
      ans <- asy_cov_knn(diag_test, tcf_thet, tcf_bet, cp, pi_knn_est,
                         rho_knn_est, rho_est$K)
    } else {
      bst_knn <- function(dt, inds, cp, k, type) {
        dat <- dt[inds, ]
        x_mat <- as.matrix(dat[, -c(1:5)])
        rho_knn_est <- rho_knn(x_mat, as.matrix(dat[, c(2:4)]), dat[, 5], k = k,
                               type = type)
        rocs.tcf(method = "knn", diag_test = dat[, 1],
                 dise_vec = as.matrix(dat[, c(2:4)]),
                 veri_stat = dat[, 5], rho_est = rho_knn_est, cps = cp)
      }
      data <- data.frame(diag_test, dise_vec, veri_stat, rho_est$X)
      if (parallel) {
        res_bst <- boot(data, bst_knn, R = n_boot, cp = cp, k = rho_est$K,
                        type = rho_est$type, parallel = "snow", ncpus = ncpus)
      } else {
        res_bst <- boot(data, bst_knn, R = n_boot, cp = cp, k = rho_est$K,
                        type = rho_est$type)
      }
      ans <- var(res_bst$t)
    }
  } else if (method == "fi") {
    if (!boot) {
      term1 <- est_func_ie_deriv(veri_stat, diag_test, rho_est$X, cp,
                                 rho_est$coeff, rho_est$values,
                                 rho_est$Hess, tcf_thet, tcf_bet, m = 0)
      term2 <- est_func_ie(dise_vec_temp, veri_stat, diag_test, rho_est$X, cp,
                           rho_est$coeff, rho_est$values, tcf_thet,
                           tcf_bet, m = 0)
      term1_inv <- solve(term1)
      sig <- term1_inv %*% (t(term2) %*% term2) %*% t(term1_inv)
      hh <- h_deriv(tcf_thet, tcf_bet, length(rho_est$coeff))
      ans <- hh %*% sig %*% t(hh)
    } else {
      bst_fi <- function(dt, inds, formula, cp) {
        dat <- dt[inds, ]
        out <- rho_mlogit(formula, data = dat)
        rocs.tcf(method = "fi", diag_test = dat[, 1],
                 dise_vec = as.matrix(dat[, c("D1", "D2", "D3")]),
                 rho_est = out, cps = cp)
      }
      name_var <- as.character(attributes(terms(rho_est$formula))$variables)[-1]
      data <- data.frame(diag_test, rho_est$D, rho_est$X[, name_var[-1]],
                         dise_vec)
      names(data) <- c("Temp", name_var, "D1", "D2", "D3")
      if (parallel) {
        res_bst <- boot(data, bst_fi, R = n_boot, formula = rho_est$formula,
                        cp = cp, parallel = "snow", ncpus = ncpus)
      } else {
        res_bst <- boot(data, bst_fi, R = n_boot, formula = rho_est$formula,
                        cp = cp)
      }
      ans <- var(res_bst$t)
    }
  } else if (method == "msi") {
    if (!boot) {
      term1 <- est_func_ie_deriv(veri_stat, diag_test, rho_est$X, cp,
                                 rho_est$coeff, rho_est$values, rho_est$Hess,
                                 tcf_thet, tcf_bet, m = 1)
      term2 <- est_func_ie(dise_vec_temp, veri_stat, diag_test, rho_est$X, cp,
                           rho_est$coeff, rho_est$values, tcf_thet,
                           tcf_bet, m = 1)
      term1_inv <- solve(term1)
      sig <- term1_inv %*% (t(term2) %*% term2) %*% t(term1_inv)
      hh <- h_deriv(tcf_thet, tcf_bet, length(rho_est$coeff))
      ans <- hh %*% sig %*% t(hh)
    } else {
      bst_msi <- function(dt, inds, formula, cp) {
        dat <- dt[inds, ]
        out <- rho_mlogit(formula, data = dat)
        rocs.tcf(method = "msi", diag_test = dat[, 1],
                 dise_vec = as.matrix(dat[, c("D1", "D2", "D3")]),
                 veri_stat = dat[, "V"], rho_est = out, cps = cp)
      }
      name_var <- as.character(attributes(terms(rho_est$formula))$variables)[-1]
      data <- data.frame(diag_test, rho_est$D, rho_est$X[, name_var[-1]],
                         veri_stat, dise_vec)
      names(data) <- c("Temp", name_var, "V", "D1", "D2", "D3")
      if (parallel) {
        res_bst <- boot(data, bst_msi, R = n_boot, formula = rho_est$formula,
                        cp = cp, parallel = "snow", ncpus = ncpus)
      } else {
        res_bst <- boot(data, bst_msi, R = n_boot, formula = rho_est$formula,
                        cp = cp)
      }
      ans <- var(res_bst$t)
    }
  } else if (method == "ipw") {
    if (!boot) {
      term1 <- est_func_ipw_deriv(dise_vec_temp, veri_stat, diag_test, pi_est$X,
                                  cp, pi_est$coeff, pi_est$values, pi_est$Hess,
                                  tcf_thet, tcf_bet, pi_est$model)
      term2 <- est_func_ipw(dise_vec_temp, veri_stat, diag_test, pi_est$X, cp,
                            pi_est$coeff, pi_est$values, tcf_thet, tcf_bet,
                            pi_est$model)
      term1_inv <- solve(term1)
      sig <- term1_inv %*% (t(term2) %*% term2) %*% t(term1_inv)
      hh <- h_deriv(tcf_thet, tcf_bet, length(pi_est$coeff))
      ans <- hh %*% sig %*% t(hh)
    } else {
      bst_ipw <- function(dt, inds, formula, cp) {
        dat <- dt[inds, ]
        out <- psglm(formula, data = dat, test = FALSE, trace = FALSE)
        rocs.tcf(method = "ipw", diag_test = dat[, 1],
                 dise_vec = as.matrix(dat[, c("D1", "D2", "D3")]),
                 veri_stat = dat[, 2], pi_est = out, cps = cp)
      }
      name_var <- as.character(attributes(terms(pi_est$formula))$variables)[-1]
      data <- data.frame(diag_test, veri_stat, pi_est$X[, name_var[-1]],
                         dise_vec)
      names(data) <- c("Temp", name_var, "D1", "D2", "D3")
      if (parallel) {
        res_bst <- boot(data, bst_ipw, R = n_boot, formula = pi_est$formula,
                        cp = cp, parallel = "snow", ncpus = ncpus)
      } else {
        res_bst <- boot(data, bst_ipw, R = n_boot, formula = pi_est$formula,
                        cp = cp)
      }
      ans <- var(res_bst$t)
    }
  } else {
    if (!boot) {
      term1 <- est_func_spe_deriv(dise_vec_temp, veri_stat, diag_test,
                                  rho_est$X, pi_est$X, cp, rho_est$coeff,
                                  rho_est$values, rho_est$Hess, pi_est$coeff,
                                  pi_est$values, pi_est$Hess, tcf_thet,
                                  tcf_bet, pi_est$model)
      term2 <- est_func_spe(dise_vec_temp, veri_stat, diag_test, rho_est$X,
                            pi_est$X, cp, rho_est$coeff, rho_est$values,
                            pi_est$coeff, pi_est$values, tcf_thet,
                            tcf_bet, pi_est$model)
      term1_inv <- solve(term1)
      sig <- term1_inv %*% (t(term2) %*% term2) %*% t(term1_inv)
      hh <- h_deriv(tcf_thet, tcf_bet, length(rho_est$coeff) +
                      length(pi_est$coeff))
      ans <- hh %*% sig %*% t(hh)
    } else {
      bst_spe <- function(dt, inds, formula_rho, formula_pi, cp) {
        dat <- dt[inds, ]
        out_rho <- rho_mlogit(formula_rho, data = dat)
        out_pi <- psglm(formula_pi, data = dat, test = FALSE, trace = FALSE)
        rocs.tcf(method = "spe", diag_test = dat[, 1],
                 dise_vec = as.matrix(dat[, c("D1", "D2", "D3")]),
                 veri_stat = dat[, "V"], rho_est = out_rho, pi_est = out_pi,
                 cps = cp)
      }
      name_var_rho <- as.character(
        attributes(terms(rho_est$formula))$variables)[-1]
      name_var_pi <- as.character(
        attributes(terms(pi_est$formula))$variables)[-1]
      data_rho <- data.frame(rho_est$X[, name_var_rho[-1]])
      data_pi <- data.frame(pi_est$X[, name_var_pi[-1]])
      data <- data.frame(diag_test, veri_stat, rho_est$D,
                         union(data_rho, data_pi), dise_vec)
      names(data) <- c("Temp", name_var_pi[1], name_var_rho[1],
                       union(name_var_rho[-1], name_var_pi[-1]),
                       "D1", "D2", "D3")
      if (parallel) {
        res_bst <- boot(data, bst_spe, R = n_boot,
                        formula_rho = rho_est$formula,
                        formula_pi = pi_est$formula, cp = cp, parallel = "snow",
                        ncpus = ncpus)
      } else {
        res_bst <- boot(data, bst_spe, R = n_boot,
                        formula_rho = rho_est$formula,
                        formula_pi = pi_est$formula, cp = cp)
      }
      ans <- var(res_bst$t)
    }
  }
  rownames(ans) <- colnames(ans) <- paste("TCF", c(1:3), sep = "")
  return(ans)
}
