####========================================================================####
## This file consists of some functions for computation VUS                   ##
## by using the VUS formula.                                                  ##
##                  																													##
####========================================================================####
##
#' @title Estimation methods for volume under ROC surface (VUS) under MAR
#'
#' @description \code{vus_mar} computes bias-corrected estimates of the volume under the ROC surface for evaluating the accuracy of a continuous diagnostic test.
#'
#' @param method  name of bias-corrected estimation method to be used for estimating the VUS in presence of verification bias. See \code{\link{rocs}} for more details.
#' @param diag_test  a numeric vector containing the diagnostic test values. \code{NA} values are not admitted.
#' @param dise_vec  a n * 3  binary matrix with the three columns, corresponding to three classes of the disease status. In row i, 1 in column j indicates that the i-th subject belongs to class j, with j = 1, 2, 3. A row of \code{NA} values indicates a non-verified subject.
#' @param veri_stat  a binary vector containing the verification status (1 verified, 0 not verified).
#' @param rho_est  a result of a call to \code{\link{rho_mlogit}} of \code{\link{rho_knn}} to fit the disease model.
#' @param pi_est  a result of a call to \code{\link{psglm}} to fit the verification model.
#' @param ci a logical value. If TRUE (default), computes an confidence interval of VUS and tests the null hypothesis H0: VUS = 1/6.
#' @param ci_level  an confidence level to be used for constructing the confidence interval; default 0.95.
#' @param boot a logical value. Default = \code{FALSE}. If set to \code{TRUE}, a bootstrap resampling is employed to estimate the asymptotic variance of bias-corrected VUS estimates. See \code{\link{asy_var_vus}}.
#' @param n_boot  the number of bootstrap replicates, which is used for FULL or KNN estimator, or option \code{boot = TRUE}. Usually this will be a single positive integer.
#' @param parallel  a logical value. If \code{TRUE}, a parallel computing is employed to the bootstrap resampling process.
#' @param ncpus  number of processes to be used in parallel computing. Default is a half of available cores.
#' @param trace a logical value. If \code{TRUE}, tracing information on the progress of the estimation is produced.
#'
#' @details
#' The function implements five bias-corrected estimation methods in To Duc et al. (2016, 2020) for estimating VUS of a three-class continuous diagnostic test in presence of verification bias. The estimators are full imputation (FI), mean score imputation (MSI), inverse probability weighted (IPW), semiparametric efficient (SPE) and K nearest-neighbor (KNN), see \code{\link{rocs}}. These estimators work under MAR assumption.
#'
#' The standard error of the estimates are obtained through the function \code{\link{asy_var_vus}}. In particular, the standard error of the FULL estimate is computed by bootstrap resampling method or by Jackknife approach proposed in Guangming et al. (2013). For the bias-corrected estimates, the standard errors are computed by using asymptotic theory (with respect to FI, MSI, IPW and SPE estimator) or bootstrap resampling method (with respect to KNN estimator). A confidence interval for VUS also is given. A logit transformation is also applied for obtaining the confidence interval.
#'
#' The default value of the number of bootstrap replicates is 250.
#'
#' Note that, before apply the functions \code{vus_mar}, the use of \code{\link{pre_data}} might be needed to check the monotone ordering disease classes and to create the matrix format for disease status.
#'
#' @return \code{vus_mar} returns an object of class inheriting from "vus_mar" class.
#'
#' The function \code{\link{print.vus_mar}} can be used to print a summary of the results.
#'
#' An object of class "vus_mar" is a list containing at least the following components:
#'
#' \item{vus_fit}{the estimate of VUS.}
#' \item{std}{the standard error, obtained by using asymptotic theory or bootstrap resampling method.}
#' \item{call}{the matched call.}
#' \item{t_stat}{t-statistic.}
#' \item{p_val_norm}{p-value correspond to normal-test.}
#' \item{ci_norm}{the confidence interval of VUS by using normal approximation.}
#' \item{ci_logit}{the confidence interval of VUS via logit transform.}
#' \item{ci_level}{the confidence level used.}
#' \item{boot}{the value of \code{boot}.}
#' \item{n_boot}{the number of bootstrap replicates used.}
#'
#' In addition, the name of method used to estimate VUS also is given as the attribute of \code{vus_fit}.
#'
#' @references
#' To Duc, K., Chiogna, M. and Adimari, G. (2020)
#' Nonparametric estimation of ROC surfaces in presence of verification bias.
#' \emph{REVSTAT-Statistical Journal}, \bold{18}, 5, 697–720.
#'
#' To Duc, K., Chiogna, M. and Adimari, G. (2016)
#' Bias-corrected methods for estimating the receiver operating characteristic surface of continuous diagnostic tests.
#' \emph{Electronic Journal of Statistics}, \bold{10}, 3063-3113.
#'
#' Guangming, P., Xiping, W. and Wang, Z. (2013)
#' Non-parameteric statistical inference for $P(X < Y < Z)$.
#' \emph{Sankhya A}, \bold{75}, 1, 118-138.
#'
#' @examples
#' data(EOC)
#' head(EOC)
#'
#'
#' \dontrun{
#' # FULL data estimator
#' dise_full <- pre_data(EOC$D.full, EOC$CA125)
#' dise_vec_full <- dise_full$dise_vec
#' vus_mar("full", diag_test = EOC$CA125, dise_vec = dise_vec_full)
#' }
#'
#' \dontrun{
#' # Preparing the missing disease status
#' dise_na <- pre_data(EOC$D, EOC$CA125)
#' dise_vec_na <- dise_na$dise_vec
#' dise_fact_na <- dise_na$dise
#' # FI estimator
#' rho_out <- rho_mlogit(dise_fact_na ~ CA125 + CA153 + Age, data = EOC,
#'                       test = TRUE)
#'
#' vus_mar("fi", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'         veri_stat = EOC$V, rho_est = rho_out)
#'
#' # MSI estimator
#' vus_mar("msi", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'         veri_stat = EOC$V, rho_est = rho_out)
#'
#' # IPW estimator
#' pi_out <- psglm(V ~ CA125 + CA153 + Age, data = EOC, test = TRUE)
#' vus_mar("ipw", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'         veri_stat = EOC$V, pi_est = pi_out)
#'
#' # SPE estimator
#' vus_mar("spe", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'         veri_stat = EOC$V, rho_est = rho_out, pi_est = pi_out)
#'
#' # KNN estimator, K = 1, Mahalanobis distance
#' x_mat <- cbind(EOC$CA125, EOC$CA153, EOC$Age)
#' rho_maha_1nn <- rho_knn(x_mat = x_mat, dise_vec = dise_vec_na,
#'                         veri_stat = EOC$V, k = 1, type = "mahala")
#' vus_mar("knn", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'         veri_stat = EOC$V, rho_est = rho_maha_1nn)
#' }
#'

##
#' @import utils
#' @import parallel
#' @importFrom Rcpp evalCpp
#' @useDynLib bcROCsurface, .registration = TRUE
#' @export
vus_mar <- function(method = "full", diag_test, dise_vec, veri_stat,
                    rho_est = NULL, pi_est = NULL, ci = TRUE,
                    ci_level = ifelse(ci, 0.95, NULL), boot = FALSE,
                    n_boot = ifelse(ci, 250, NULL), parallel = FALSE,
                    ncpus = ifelse(parallel, detectCores() / 2, NULL),
                    trace = TRUE) {
  ## return the function call
  call <- match.call()
  ## checking the method
  method_temp <- substitute(me, list(me = method))
  ok_method <- c("full", "fi", "msi", "ipw", "spe", "knn")
  if (length(method_temp) > 1)
    stop(gettextf("Please, choose one method from %s",
                  paste(sQuote(ok_method), collapse = ", ")), domain = NA)
  if (!is.character(method_temp)) method_temp <- deparse(method_temp)
  if (!is.element(method_temp, ok_method)) {
    stop(gettextf("the required method should be one of %s", method_temp,
                  paste(sQuote(ok_method), collapse = ", ")),
         domain = NA)
  }
  ## checking the argument diag_test
  if (missing(diag_test)) stop("argument \"diag_test\" is missing \n")
  if (!inherits(diag_test, "numeric") || any(is.na(diag_test)))
    stop("\"diag_test\" must be a numeric vector and not include NA values")
  name_diagnostic <- substitute(diag_test)
  if (!is.character(name_diagnostic)) {
    name_diagnostic <- deparse(name_diagnostic)
  }
  name_diagnostic <- unlist(strsplit(name_diagnostic, NULL))
  if (any(name_diagnostic %in% c("$"))) {
    id_name <- which(name_diagnostic %in% c("$"))
    name_diagnostic <- paste(
      name_diagnostic[(id_name[1] + 1) : length(name_diagnostic)],
      collapse = "")
  } else {
    name_diagnostic <- paste(name_diagnostic, collapse = "")
  }
  method_name <- toupper(method)
  ## checking dise_vec
  if (missing(dise_vec)) stop("argument \"dise_vec\" is missing \n")
  if (!inherits(dise_vec, "matrix") || ncol(dise_vec) != 3 ||
      !all(is.element(na.omit(dise_vec), c(0, 1))))
    stop("variable \"dise_vec\" must be a binary matrix with 3 columns")
  if (length(diag_test) != nrow(dise_vec))
    stop(gettextf("arguments imply differing number of observation: %d",
                  length(diag_test)), gettextf(", %d", nrow(dise_vec)),
         domain = NA)
  dise_vec_flag <- any(is.na(dise_vec))
  ##
  if (!is.null(rho_est)) {
    if (!is.element(class(rho_est), c("prob_dise", "prob_dise_knn")))
      stop("\"rho_est\" not a result of rhoMlogit or rho_knn \n")
    if (is.element(class(rho_est), c("prob_dise_knn"))) {
      if (is.null(rho_est$K))
        stop("Please, choose one value of K \n")
    }
  }
  if (!is.null(pi_est)) {
    if (!is.element(class(pi_est), c("prob_veri")))
      stop("\"pi_est\" not a result of psglm \n")
  }
  ## checking veri_stat
  if (missing(veri_stat)) {
    if (method_temp != "full" || dise_vec_flag)
      stop("\"veri_stat\" is missing")
    if (trace) {
      cat("Hmm, look likes the full data.\n")
      cat("The verification status is not available.\n")
      cat("You are working on FULL or Complete Case approach.\n")
      cat("Number of observation:", length(diag_test), "\n")
      cat("The diagnostic test:", name_diagnostic, "\n")
      cat("Processing .... \n")
      flush.console()
    }
    veri_stat <- NULL
  } else {
    if (all(veri_stat == 0) || !all(is.element(veri_stat, c(0, 1))))
      stop("\"veri_stat\" must have values as 0 and 1.")
    if (nrow(dise_vec) != length(veri_stat) ||
        length(diag_test) != length(veri_stat))
      stop(gettextf("arguments imply differing number of observation: %d",
                    length(diag_test)),
           gettextf(", %d", nrow(dise_vec)), domain = NA)
    if (all(veri_stat == 1)) {
      if (method_temp != "full" || dise_vec_flag)
        stop("all \"veri_stat\" are 1 but \"dise_vec\" has NA!\n")
      if (trace) {
        cat("Hmm, look likes the full data\n")
        cat("Number of observation:", length(diag_test), "\n")
        cat("All subjects underwent the verification process\n")
        cat("You are working on FULL or Complete Case approach\n")
        cat("The diagnostic test:", name_diagnostic, "\n")
        cat("Processing .... \n")
        flush.console()
      }
    } else {
      rv <- mean(veri_stat)
      if (!dise_vec_flag) {
        if (trace) {
          cat("Warning: There are no NA values in variable dise_vec, while",
              paste(round(rv * 100), "%", sep = ""),
              "of the subjects receive disease verification. \n")
          cat("BE CAREFULL OF YOUR INPUT AND RESULTS \n")
          cat("Number of observation:", length(diag_test), "\n")
          cat("You required estimate VUS using", method_name, "approach \n")
          cat("The diagnostic test:", name_diagnostic, "\n")
          cat("Processing .... \n")
          flush.console()
        }
      } else{
        if (trace) {
          cat("Hmm, look likes the incomplete data\n")
          cat("Number of observation:", length(diag_test), "\n")
          cat(paste(round(rv * 100), "%", sep = ""),
              "of the subjects receive disease verification. \n")
          cat("You required estimate VUS using", method_name, "approach \n")
          cat("The diagnostic test:", name_diagnostic, "\n")
          cat("Processing .... \n")
          flush.console()
        }
      }
    }
  }
  ## Main body
  if (method_temp == "full") {
    if (dise_vec_flag) {
      cat("Opp! Look likes wrong method \n")
      ques <- readline("Do you want use Complete Case (CC) approach? [y/n]: ")
      if (ques %in% c("y", "n")) {
        if (ques == "y") {
          if (trace) {
            cat("Number of observation:", length(diag_test), "\n")
            cat("We are estimating VUS by using CC method \n")
            cat("BE CAREFULL OF THE RESULTS. THIS CAN MAKE THE DISTORTED INFERENCE IN VUS \n")
            cat("Processing .... \n")
            flush.console()
          }
          dise_vec <- dise_vec[veri_stat == 1, ]
          diag_test <- diag_test[veri_stat == 1]
          ans <- vus_c(diag_test, dise_vec)
        } else if (ques == "n") {
          stop("The FULL method cannot be performed if disease status have NA value(s).")
        }
      } else {
        stop("The answer was wrong. Please, choose y for Yes, or n for No!")
      }
    }
    ans <- vus_c(diag_test, dise_vec)
  } else if (method_temp == "fi") {
    if (is.null(rho_est))
      stop("The input of argument \"rho_est\" is needed for ",
           method_name, " estimator")
    if (length(diag_test) != nrow(rho_est$values))
      stop(gettextf("arguments imply differing number of observation: %d",
                    length(diag_test), ", %d", nrow(rho_est$values)),
           domain = NA)
    ans <- vus_c(diag_test, rho_est$values)
  } else if (method_temp %in% c("msi", "knn")) {
    if (is.null(rho_est))
      stop("argument \"rho_est\" is needed for", method_name, "estimator")
    dise_vec_temp <- dise_vec
    if (dise_vec_flag) dise_vec_temp[is.na(dise_vec)] <- 99
    dise_msi <- dise_vec_temp * veri_stat + (1 - veri_stat) * rho_est$values
    ans <- vus_c(diag_test, dise_msi)
  } else if (method_temp == "ipw") {
    if (is.null(pi_est))
      stop("argument \"pi_est\" is needed for", method_name, "estimator")
    dise_vec_temp <- dise_vec
    if (dise_vec_flag) dise_vec_temp[is.na(dise_vec)] <- 99
    dise_ipw <- dise_vec_temp * veri_stat / pi_est$values
    ans <- vus_c(diag_test, dise_ipw)
  } else{
    if (is.null(rho_est) || is.null(pi_est))
      stop("arguments \"rho_est\" and \"pi_est\" are needed for",
           method_name, "estimator")
    dise_vec_temp <- dise_vec
    if (dise_vec_flag) dise_vec_temp[is.na(dise_vec)] <- 99
    dise_spe <- dise_vec_temp * veri_stat / pi_est$values -
      (veri_stat / pi_est$values - 1) * rho_est$values
    ans <- vus_c(diag_test, dise_spe)
  }
  attr(ans, "name") <- method_name
  res <- list(vus_fit = ans, call = call)
  class(res) <- "vus_mar"
  if (ci) {
    var_ans <- asy_var_vus(res, diag_test = diag_test, dise_vec = dise_vec,
                           veri_stat = veri_stat, rho_est = rho_est,
                           pi_est = pi_est, boot = boot, n_boot = n_boot,
                           parallel = parallel, ncpus = ncpus)
    cf_ans <- ans + c(-1, 1) * qnorm((1 + ci_level) / 2) * sqrt(var_ans)
    sd_log <- sqrt(var_ans) / (ans * (1 - ans))
    log_cf <- log(ans / (1 - ans)) +
      c(-1, 1) * qnorm((1 + ci_level) / 2) * sd_log
    cf_ans_tran <- exp(log_cf) / (1 + exp(log_cf))
    w <- (ans - 1 / 6) / sqrt(var_ans)
    p_val_norm <- 1 - pnorm(w)
    res <- list(vus_fit = ans, std = sqrt(var_ans), t_stat = w,
                p_val_norm = p_val_norm, ci_norm = cf_ans,
                ci_logit = cf_ans_tran, call = call,
                ci_level = ci_level, boot = boot, n_boot = n_boot)
    class(res) <- "vus_mar"
  }
  if (trace) cat("DONE\n")
  res
}


## The function print.vus_mar
#' @title Print summary results of VUS
#'
#' @description \code{print.vus_mar} prints the results for the output of function \code{\link{vus_mar}}.
#'
#' @method print vus_mar
#' @param x an object of class "vus_mar", a result of a call to \code{\link{vus_mar}}.
#' @param digits minimal number of significant digits, see \code{\link{print.default}}.
#' @param ... further arguments passed to \code{\link{print}} method.
#'
#' @details \code{print.vus_mar} shows a nice format of the summary table for the VUS estimate results. Some information on the diagnostic test, the fitted values of VUS, and confidence intervals are shown.
#'
#' @seealso \code{\link{vus_mar}}
#'
#' @export
print.vus_mar <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\n")
  cat("CALL: ",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n \n", sep = "")
  if (!is.null(x$std)) {
    res_tab <- c(x$ci_level * 100, x$ci_norm, x$ci_logit)
    res_tab <- format(round(res_tab, digits = digits))
    res_tab[1L] <- paste("\n", x$ci_level * 100, "% ", sep = "")
    res_tab[2 * (1L:2L)] <- paste(" (", res_tab[2 * (1L:2L)], ",", sep = "")
    res_tab[2 * (1L:2L) + 1L] <- paste(res_tab[2 * (1L:2L) + 1L], ") ")
    p_val <- x$p_val_norm
    stat <- x$t_stat
    test_tab <- cbind(stat, p_val)
    colnames(test_tab) <- c("Test Statistic", "P-value")
    rownames(test_tab) <- "Normal-test"
    ci_name <- c("       Normal        ", "       Logit        ")
    cat("Estimate of VUS:", format(round(x$vus_fit, digits = digits)), "\n")
    cat("Standard error:", format(round(x$std, digits = digits)), "\n")
    cat("\nIntervals:")
    cat("\nLevel", ci_name)
    cat(res_tab)
    if (attr(x$vus, "name") %in% c("FULL", "KNN")) {
      if (x$boot)
        cat("\nEstimation of Standard Error and Intervals are based on Bootstrap with",
            x$n_boot, "replicates\n")
      else
        cat("\nEstimation of Standard Error and Intervals are based on Jackknife approach \n")
    } else {
      cat("\nEstimation of Standard Error and Intervals are based on Asymptotic Theory \n")
    }
    cat("\n")
    cat("Testing the null hypothesis H0: VUS = 1/6 \n")
    printCoefmat(test_tab, has.Pvalue = TRUE)
  } else {
    cat("Estimate of VUS:", round(x$vus_fit, digits = digits), "\n")
  }
  invisible(x)
}
