####========================================================================####
## This file consists of some functions that are related to the computation   ##
## asymptotic variance of VUS                                                 ##
##                																														##
####========================================================================####
##
#' @title Asymptotic variance estimation for VUS
#'
#' @description \code{asyVarVUS} computes the asymptotic variance of full data (FULL) and bias-corrected estimators (i.e. full imputation, mean score imputation, inverse probability weighting, semiparametric efficient and K nearest neighbor) of VUS.
#'
#' @param obj_vus  a result of a call to \code{\link{vus}}.
#' @param Test  a numeric vector containing the diagnostic test values. \code{NA} values of \code{Test} are not accepted.
#' @param Dvec  a n * 3  binary matrix with three columns, corresponding to the three classes of the disease status. In row i, 1 in column j indicates that the i-th subject belongs to class j, with j = 1, 2, 3. A row of \code{NA} values indicates a non-verified subject.
#' @param V  a binary vector containing the verification status (1 verified, 0 not verified).
#' @param rhoEst  a result of a call to \code{\link{rhoMLogit}} of \code{\link{rhoKNN}} to fit the disease model.
#' @param piEst  a result of a call to \code{\link{psglm}} to fit the verification model.
#' @param BOOT a logical value. Default = \code{FALSE}. If set to \code{TRUE}, a bootstrap resampling is employed to estimate the asymptotic variance of the bias-corrected VUS estimators.
#' @param nR  the number of bootstrap replicates, which is used for FULL or KNN estimators, or option \code{BOOT = TRUE}. The defaut is 250.
#' @param parallel  a logical value. If \code{TRUE}, a parallel computing is employed in the bootstrap resampling process.
#' @param ncpus  number of processes to be used in parallel computing. Default is half of available cores.
#'
#'
#' @details
#' For the FULL estimator, a bootstrap resampling process or Jackknife approach is used to estimate the asymptotic variance, whereas, a bootstrap resampling process is employed to obtain the asymptotic variance of K nearest neighbor estimator.
#'
#' For the full imputation, mean score imputation, inverse probability weighting and semiparametric efficient estimators of VUS, the asymptotic variances are computed by using the explicit form. Furthermore, a bootstrap procedure is also available, useful in case of small sample sizes.
#'
#' @return \code{asyVarVUS} returns a estimated value of the asymptotic variance.
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
#' Dna <- preDATA(EOC$D, EOC$CA125)
#' Dfact_na <- Dna$D
#' Dvec_na <- Dna$Dvec
#'
#' rho_out <- rhoMLogit(Dfact_na ~ CA125 + CA153 + Age, data = EOC, test = TRUE)
#' vus_fi <- vus("fi", Test = EOC$CA125, Dvec = Dvec_na, V = EOC$V,
#'               rhoEst = rho_out, ci = FALSE)
#' var_fi <- asyVarVUS(vus_fi, Test = EOC$CA125, Dvec = Dvec_na, V = EOC$V,
#'                     rhoEst = rho_out)
#'
#' \dontrun{
#' var_bst_spe <- asyVarVUS(vus_spe, Test = EOC$CA125, Dvec = Dvec_na, V = EOC$V,
#'                          rhoEst = rho_out, piEst = pi_out, BOOT = TRUE,
#'                          parallel = TRUE)
#' }
#'
#'
#' @importFrom Rcpp evalCpp
#' @import parallel
#' @import boot
#' @export
asyVarVUS <- function(obj_vus, Test, Dvec, V = NULL, rhoEst = NULL,
                      piEst = NULL, BOOT = FALSE, nR = 250, parallel = FALSE,
                      ncpus = ifelse(parallel, detectCores() / 2, NULL)) {
  if (!inherits(obj_vus, "vus"))
    stop("The argument \"obj_vus\" is not a result of vus()")
  ## checking the argument Test
  if (missing(Test)) stop("argument \"Test\" is missing \n")
  if (!inherits(Test, "numeric") || any(is.na(Test)))
    stop("variable \"Test\" must be a numeric vector and not include NA values")
  ## checking Dvec
  if (missing(Dvec)) stop("argument \"Dvec\" is missing \n")
  if (!inherits(Dvec, "matrix") || ncol(Dvec) != 3 ||
        !all(is.element(na.omit(Dvec), c(0, 1))))
    stop("variable \"Dvec\" must be a binary matrix with 3 columns")
  if (length(Test) != nrow(Dvec))
    stop(gettextf("arguments imply differing number of observation: %d",
                  length(Test)), gettextf(", %d", nrow(Dvec)), domain = NA)
  VUS_obj <- obj_vus$vus_fit
  meth <- tolower(attr(VUS_obj, "name"))
  if (!is.element(meth, c("full", "fi")) && is.null(V))
    stop("argument \"V\" is missing \n")
  if (is.element(meth, c("fi", "msi", "knn", "spe")) && is.null(rhoEst))
    stop("argument \"rhoEst\" is missing \n")
  if (is.element(meth, c("ipw", "spe")) && is.null(piEst))
    stop("argument \"piEst\" is missing \n")
  ##
  if (!is.null(rhoEst)) {
    if (!is.element(class(rhoEst), c("prob_dise", "prob_dise_knn")))
      stop("\"rhoEst\" not a result of rhoMlogit or rhoKNN \n")
    if (is.element(class(rhoEst), c("prob_dise_knn"))) {
      if (is.null(rhoEst$K))
        stop("Please, choose one value of K\n")
    }
  }
  if (!is.null(piEst)) {
    if (!is.element(class(piEst), c("prob_veri")))
      stop("\"piEst\" not a result of psglm \n")
  }
  n <- length(Test)
  Dvec_temp <- Dvec
  if (any(is.na(Dvec))) Dvec_temp[is.na(Dvec)] <- 99
  if (meth == "fi") {
    if (!BOOT) {
      score <- mlogEstFunc(Dvec_temp, V, rhoEst$X, rhoEst$values)
      hess <- solve(rhoEst$Hess)
      der_rho1 <- t(rho_deriv(rhoEst$X, rhoEst$values, ref_level = "1"))
      der_rho2 <- t(rho_deriv(rhoEst$X, rhoEst$values, ref_level = "2"))
      der_rho3 <- - (der_rho1 + der_rho2)
      Q_Fi <- asyVarVUS_C(Test, rhoEst$values, VUS_obj, score, hess, der_rho1,
                          der_rho2, der_rho3)
      ans_var <- n ^ 4 * sum(Q_Fi ^ 2) / prod(colSums(rhoEst$values) ^ 2)
    } else {
      bst_fi <- function(dt, inds, formula) {
        dat <- dt[inds, ]
        out <- rhoMLogit(formula, data = dat)
        vusC(dat[, 1], out$values)
      }
      name_var <- as.character(attributes(terms(rhoEst$formula))$variables)[-1]
      data <- data.frame(Test, rhoEst$D, rhoEst$X[, name_var[-1]])
      names(data) <- c("Temp", name_var)
      if (parallel) {
        res_bst <- boot(data, bst_fi, R = nR, formula = rhoEst$formula,
                        parallel = "snow", ncpus = ncpus)
      } else {
        res_bst <- boot(data, bst_fi, R = nR, formula = rhoEst$formula)
      }
      ans_var <- as.numeric(var(res_bst$t, na.rm = TRUE))
    }
  } else if (meth == "msi") {
    if (!BOOT) {
      D_MSI <- as.matrix(V * Dvec_temp + (1 - V) * rhoEst$values)
      score <- mlogEstFunc(Dvec_temp, V, rhoEst$X, rhoEst$values)
      hess <- solve(rhoEst$Hess)
      der_D_MSI1 <- t((1 - V) * rho_deriv(rhoEst$X, rhoEst$values,
                                          ref_level = "1"))
      der_D_MSI2 <- t((1 - V) * rho_deriv(rhoEst$X, rhoEst$values,
                                          ref_level = "2"))
      der_D_MSI3 <- - (der_D_MSI1 + der_D_MSI2)
      Q_Msi <- asyVarVUS_C(Test, D_MSI, VUS_obj, score, hess, der_D_MSI1,
                           der_D_MSI2, der_D_MSI3)
      ans_var <- n ^ 4 * sum(Q_Msi ^ 2) / prod(colSums(D_MSI) ^ 2)
    } else {
      bst_msi <- function(dt, inds, formula) {
        dat <- dt[inds, ]
        out <- rhoMLogit(formula, data = dat)
        Dmsi <- as.matrix(dat[, c("D1", "D2", "D3")] * dat[, "V"] +
                            (1 - dat[, "V"]) * out$values)
        vusC(dat[, 1], Dmsi)
      }
      name_var <- as.character(attributes(terms(rhoEst$formula))$variables)[-1]
      data <- data.frame(Test, rhoEst$D, rhoEst$X[, name_var[-1]], V, Dvec_temp)
      names(data) <- c("Temp", name_var, "V", "D1", "D2", "D3")
      if (parallel) {
        res_bst <- boot(data, bst_msi, R = nR, formula = rhoEst$formula,
                        parallel = "snow", ncpus = ncpus)
      } else {
        res_bst <- boot(data, bst_msi, R = nR, formula = rhoEst$formula)
      }
      ans_var <- as.numeric(var(res_bst$t, na.rm = TRUE))
    }
  } else if (meth == "ipw") {
    if (!BOOT) {
      D_IPW <- as.matrix(V * Dvec_temp / piEst$values)
      score <- piEstFunc(V, piEst$X, piEst$values, piEst$coeff, piEst$model)
      hess <- solve(piEst$Hess)
      der_pi_inv <- pi_inv_deriv(piEst$X, piEst$values, piEst$coeff,
                                 piEst$model)
      der_D_IPW1 <- t(V * Dvec_temp[, 1] * der_pi_inv)
      der_D_IPW2 <- t(V * Dvec_temp[, 2] * der_pi_inv)
      der_D_IPW3 <- t(V * Dvec_temp[, 3] * der_pi_inv)
      Q_Ipw <- asyVarVUS_C(Test, D_IPW, VUS_obj, score, hess, der_D_IPW1,
                           der_D_IPW2, der_D_IPW3)
      ans_var <- sum(Q_Ipw ^ 2) / prod((colSums(D_IPW) /
                                          sum(V / piEst$values)) ^ 2) / n ^ 2
    } else {
      bst_ipw <- function(dt, inds, formula) {
        dat <- dt[inds, ]
        out <- psglm(formula, data = dat, test = FALSE, trace = FALSE)
        Dipw <- as.matrix(dat[, c("D1", "D2", "D3")] * dat[, 2] / out$values)
        vusC(dat[, 1], Dipw)
      }
      name_var <- as.character(attributes(terms(piEst$formula))$variables)[-1]
      data <- data.frame(Test, V, piEst$X[, name_var[-1]], Dvec_temp)
      names(data) <- c("Temp", name_var, "D1", "D2", "D3")
      if (parallel) {
        res_bst <- boot(data, bst_ipw, R = nR, formula = piEst$formula,
                        parallel = "snow", ncpus = ncpus)
      } else {
        res_bst <- boot(data, bst_ipw, R = nR, formula = piEst$formula)
      }
      ans_var <- as.numeric(var(res_bst$t, na.rm = TRUE))
    }
  } else if (meth == "spe") {
    if (!BOOT) {
      temp1 <- V / piEst$values
      D_SPE <- as.matrix(temp1 * Dvec_temp - (temp1 - 1) * rhoEst$values)
      score1 <- mlogEstFunc(Dvec_temp, V, rhoEst$X, rhoEst$values)
      score2 <- piEstFunc(V, piEst$X, piEst$values, piEst$coeff, piEst$model)
      score <- cbind(score1, score2)
      hess1 <- cbind(solve(rhoEst$Hess), matrix(0, nrow = nrow(rhoEst$Hess),
                                                ncol = ncol(piEst$Hess)))
      hess2 <- cbind(matrix(0, nrow = nrow(piEst$Hess),
                            ncol = ncol(rhoEst$Hess)),
                     solve(piEst$Hess))
      hess <- rbind(hess1, hess2)
      der_rho1 <- rho_deriv(rhoEst$X, rhoEst$values, ref_level = "1")
      der_rho2 <- rho_deriv(rhoEst$X, rhoEst$values, ref_level = "2")
      der_pi_inv <- pi_inv_deriv(piEst$X, piEst$values, piEst$coeff,
                                 piEst$model)
      der_D_SPE1 <- t(cbind(-(temp1 - 1) * der_rho1,
                            V * (Dvec_temp[, 1] - rhoEst$values[, 1]) *
                              der_pi_inv))
      der_D_SPE2 <- t(cbind(-(temp1 - 1) * der_rho2,
                            V * (Dvec_temp[, 2] - rhoEst$values[, 2]) *
                              der_pi_inv))
      der_D_SPE3 <- -(der_D_SPE1 + der_D_SPE2)
      Q_Spe <- asyVarVUS_C(Test, D_SPE, VUS_obj, score, hess, der_D_SPE1,
                           der_D_SPE2, der_D_SPE3)
      ans_var <- n ^ 4 * sum(Q_Spe ^ 2) / prod(colSums(D_SPE) ^ 2)
    } else {
      bst_spe <- function(dt, inds, formula_rho, formula_pi) {
        dat <- dt[inds, ]
        out_rho <- rhoMLogit(formula_rho, data = dat)
        out_pi <- psglm(formula_pi, data = dat, test = FALSE, trace = FALSE)
        Dspe <- as.matrix(dat[, c("D1", "D2", "D3")] * dat[, "V"] /
                            out_pi$values - (dat[, "V"] / out_pi$values - 1) *
                            out_rho$values)
        vusC(dat[, 1], Dspe)
      }
      name_var_rho <- as.character(
        attributes(terms(rhoEst$formula))$variables)[-1]
      name_var_pi <- as.character(
        attributes(terms(piEst$formula))$variables)[-1]
      data_rho <- data.frame(rhoEst$X[, name_var_rho[-1]])
      data_pi <- data.frame(piEst$X[, name_var_pi[-1]])
      data <- data.frame(Test, V, rhoEst$D, union(data_rho, data_pi), Dvec_temp)
      names(data) <- c("Temp", name_var_pi[1], name_var_rho[1],
                       union(name_var_rho[-1], name_var_pi[-1]), "D1", "D2",
                       "D3")
      if (parallel) {
        res_bst <- boot(data, bst_spe, R = nR, formula_rho = rhoEst$formula,
                        formula_pi = piEst$formula, parallel = "snow",
                        ncpus = ncpus)
      } else {
        res_bst <- boot(data, bst_spe, R = nR, formula_rho = rhoEst$formula,
                        formula_pi = piEst$formula)
      }
      ans_var <- as.numeric(var(res_bst$t, na.rm = TRUE))
    }
  } else if (meth == "knn") {
    bst_knn <- function(dt, inds, k, type) {
      dat <- dt[inds, ]
      XX <- as.matrix(dat[, -c(1:5)])
      dise_bt <- as.matrix(dat[, c(2:4)])
      dise_bt[dise_bt == 99] <- NA
      rho_knn <- rhoKNN(XX, dise_bt, dat[, 5], K = k, type = type)
      Dknn <- as.matrix(dat[, c(2:4)] * dat[, 5] +
                          (1 - dat[, 5]) * rho_knn$values)
      vusC(dat[, 1], Dknn)
    }
    data <- data.frame(Test, Dvec_temp, V, rhoEst$X)
    if (parallel) {
      res_bst <- boot(data, bst_knn, R = nR, k = rhoEst$K, type = rhoEst$type,
                      parallel = "snow", ncpus = ncpus)
    } else {
      res_bst <- boot(data, bst_knn, R = nR, k = rhoEst$K,
                      type = rhoEst$type)
    }
    ans_var <- as.numeric(var(res_bst$t, na.rm = TRUE))
  } else {
    if (!BOOT) {
      Test1 <- Test[Dvec[, 1] == 1]
      Test2 <- Test[Dvec[, 2] == 1]
      Test3 <- Test[Dvec[, 3] == 1]
      n1 <- length(Test1)
      n2 <- length(Test2)
      n3 <- length(Test3)
      vus_est_mult <- VUS_obj * n1 * n2 * n3
      var_term <- vusC_full_core(Test1, Test2, Test3)
      var_term_i <- sapply(1:n1, function(x) sum(var_term$ind1[-x]))
      var_term_j <- sapply(1:n2, function(x) sum(var_term$ind2[-x]))
      var_term_k <- sapply(1:n3, function(x) sum(var_term$ind3[-x]))
      ans_var <- var((vus_est_mult - var_term_i) / (n2 * n3)) / n1 +
        var((vus_est_mult - var_term_j) / (n1 * n3)) / n2 +
        var((vus_est_mult - var_term_k) / (n1 * n2)) / n3
    } else {
      bst_full <- function(dt, inds) {
        dat <- dt[inds, ]
        vusC(dat[, 1], as.matrix(dat[, c(2:4)]))
      }
      data <- data.frame(Test, Dvec_temp)
      if (parallel) {
        res_bst <- boot(data, bst_full, R = nR, parallel = "snow",
                        ncpus = ncpus)
      } else {
        res_bst <- boot(data, bst_full, R = nR)
      }
      ans_var <- as.numeric(var(res_bst$t, na.rm = TRUE))
    }
  }
  return(ans_var)
}
