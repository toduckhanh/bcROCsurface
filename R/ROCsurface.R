####========================================================================####
## The R code for the true class fractions.                                   ##
##                                                                            ##
####========================================================================####
#' @name ROCsurface
#'
#' @title Receiver operating characteristics surface for a continuous diagnostic test
#' @aliases rocs
#' @aliases rocs.tcf
#'
#'
#' @description
#' \code{rocs.tcf} is used to obtain bias-corrected estimates of the true class fractions (TCFs) for evaluating the accuracy of a continuous diagnostic test for a given cut point \eqn{(c_1, c_2)}, with \eqn{c_1 < c_2}.
#'
#' \code{rocs} provides bias-corrected estimates of the ROC surfaces of the continuous diagnostic test by using TCF.
#'
#'
#' @param method  a estimation method to be used for estimating the true class fractions in presence of verification bias. See 'Details'.
#' @param diag_test  a numeric vector containing the diagnostic test values. \code{NA} values are not allowed.
#' @param dise_vec  a n * 3  binary matrix with the three columns, corresponding to three classes of the disease status. In row i, 1 in column j indicates that the i-th subject belongs to class j, with j = 1, 2, 3. A row of \code{NA} values indicates a non-verified subject.
#' @param veri_stat  a binary vector containing the verification status (1 verified, 0 not verified).
#' @param rho_est  a result of a call to \code{\link{rho_mlogit}} of \code{\link{rho_knn}} to fit the disease model.
#' @param pi_est  a result of a call to \code{\link{psglm}} to fit the verification model.
#' @param ncp  the dimension of cut point grid. It is used to determine the cut points (see 'Details'). Default 100.
#' @param plot  if \code{TRUE}(the default), a 3D plot of ROC surface is produced.
#' @param ellipsoid  a logical value. If TRUE, adds an ellipsoidal confidence region for TCFs at a specified cut point to current plot of ROC surface.
#' @param cpst  a specified cut point, which used to construct the ellipsoid confidence region. If \code{m} ellipsoid confidence regions are required, \code{cpst} must be matrix with \code{m} rows and 2 columns. Default \code{NULL}.
#' @param ci_level  an confidence level to be used for constructing the ellipsoid confidence region; default 0.95.
#' @param surf_col  color to be used for plotting ROC surface and ellipsoid.
#' @param ...  optional arguments to be passed to \code{\link[rgl]{plot3d}}, \code{\link[rgl]{surface3d}}.
#' @param cps  a cut point \eqn{(c_1, c_2)}, with \eqn{c_1 < c_2}, which used to estimate TCFs. If \code{m} estimates of TCFs are required, \code{cps} must be matrix with \code{m} rows and 2 columns.
#' @param boot a logical value. Default = \code{FALSE}. If set to \code{TRUE}, a bootstrap resampling is employed to estimate the asymptotic variance-covariance matrix of TCFs at the cut point \code{cpst}. See more details in \code{\link{asy_cov_tcf}}.
#' @param n_boot  the number of bootstrap replicates, which is used for FULL estimator, or option \code{boot = TRUE}. Usually this will be a single positive integer. Default 250.
#' @param parallel  a logical value. If \code{TRUE}, a parallel computing is employed to the bootstrap resampling process.
#' @param ncpus  number of processes to be used in parallel computing. Default is half of of available cores.
#'
#'
#' @details In a three-class diagnostic problem, quantities used to evaluate the accuracy of a diagnostic test are the true class fractions (TCFs). For a given pair of cut points \eqn{(c_1, c_2)} such that \eqn{c_1 < c_2}, subjects are classified into class 1 (\eqn{D_1}) if \eqn{T < c_1}; class 2 (\eqn{D_2}) if \eqn{c_1 \le T < c_2}; class 3 (\eqn{D_3}) otherwise. The true class fractions of the test \eqn{T} at \eqn{(c_1, c_2)} are defined as
#' \deqn{TCF_1(c_1) = P(T < c_1| D_1 = 1) = 1 - P(T \ge c_1| D_1 = 1),}
#' \deqn{TCF_2(c_1, c_2) = P(c_1 \le T < c_2| D_2 = 1) = P(T \ge c_1| D_2 = 1) - P(T \ge c_2| D_2 = 1),}
#' \deqn{TCF_3(c_2) = P(T > c_2| D_3 = 1) = P(T \ge c_2| D_3 = 1). }
#'
#' The receiver operating characteristic (ROC) surface is the plot of \eqn{TCF_1}, \eqn{TCF_2} and \eqn{TCF_3} by varying the cut point \eqn{(c_1, c_2)} in the domain of the diagnostic test. The cut points \eqn{(c_1, c_2)} are produced by designing a cut point grid with \code{ncp} dimension. In this grid, the points satisfying \eqn{c_1 < c_2} are selected as the cut points. The number of the cut points are obtained as \eqn{ncp(ncp - 1)/2}, for example, the default is 4950.
#'
#' These functions implement the bias-corrected estimators in To Duc et al (2016, 2020) for estimating TCF of a three-class continuous diagnostic test in presence of verification bias. The estimators work under MAR assumption. Five methods are provided, namely:
#' \itemize{
#'    \item Full imputation (FI): uses the fitted values of the disease model to replace the true disease status (both of missing and non-missing values).
#'    \item Mean score imputation (MSI): replaces only the missing values by the fitted values of the disease model.
#'    \item Inverse probability weighted (IPW): weights each observation in the verification sample by the inverse of the sampling fraction (i.e. the probability that the subject was selected for verification).
#'    \item Semiparametric efficient (SPE): replaces the true disease status by the double robust estimates.
#'    \item K nearest-neighbor (KNN): uses K nearest-neighbor imputation to obtain the missing values of the true disease status.
#' }
#'
#' The argument \code{method} must be selected from the collection of the bias-corrected methods, i.e., \code{"full"}, \code{"fi"}, \code{"msi"}, \code{"ipw"}, \code{"spe"} and \code{"knn"}.
#'
#' The ellipsoidal confidence region of TCFs at a given cut point can be constructed by using a normal approximation and plotted in the ROC surface space. The confidence level (default) is 0.95.
#'
#' Note that, before using the functions \code{rocs} and \code{rocs.tcf}, the use of \code{\link{pre_data}} might be needed to check the monotone ordering disease classes and to create the matrix format for disease status.
#'
#' @return \code{rocs} returns a list, with the following components:
#' \item{vals}{the estimates of TCFs at all cut points.}
#' \item{cpoint}{the cut points are used to construct the ROC surface.}
#' \item{ncp}{dimension of the cut point grid.}
#' \item{cpst}{the cut points are used to construct the ellipsoidal confidence regions.}
#' \item{tcf}{the estimates of TCFs at the cut points \code{cpst}.}
#' \item{message}{an integer code or vector. 1 indicates the ellipsoidal confidence region is available.}
#'
#' \code{rocs.tcf} returns a vector having estimates of TCFs at a cut point when \code{cps} is a vector with two elements, or a list of estimates of TCFs at \code{m} cut points when \code{cps} is a \code{m*2} matrix. In addition, some attributes called \code{theta}, \code{beta}, \code{cp} and \code{name} are given. Here, \code{theta} is a probability vector, with 3 element, corresponding to the disease prevalence rates of three classes. \code{beta} is also a probability vector having 4 components, which are used to compute TCFs, see To Duc el al. (2016, 2020) for more details. \code{cp} is the specified cut point that is used to estimate TCFs. \code{name} indicates the method used to estimate TCFs. These attributes are required to compute the asymptotic variance-covariance matrix of TCFs at the given cut point.
#'
#' @references
#'
#' To Duc, K., Chiogna, M. and Adimari, G. (2016)
#' Bias-corrected methods for estimating the receiver operating characteristic surface of continuous diagnostic tests.
#' \emph{Electronic Journal of Statistics}, \bold{10}, 3063-3113.
#'
#' To Duc, K., Chiogna, M. and Adimari, G. (2020)
#' Nonparametric estimation of ROC surfaces in presence of verification bias.
#' \emph{REVSTAT-Statistical Journal}, \bold{18}, 5, 697–720.
#'
#'
#' @seealso \code{\link{psglm}}, \code{\link{rho_mlogit}}, \code{\link[rgl]{plot3d}}.
#'
#' @examples
#' data(EOC)
#' head(EOC)
#'
#' \dontrun{
#' # FULL data estimator
#' dise_full <- pre_data(EOC$D.full, EOC$CA125)
#' dise_vec_full <- dise_full$dise_vec
#' if(requireNamespace("webshot2", quietly = TRUE)){
#'    rocs("full", diag_test = EOC$CA125, dise_vec = dise_vec_full, ncp = 30,
#'         ellipsoid = TRUE, cpst = c(-0.56, 2.31))
#' }
#' }
#'
#' \dontrun{
#' # Preparing the missing disease status
#' dise_na <- pre_data(EOC$D, EOC$CA125)
#' dise_vec_na <- dise_na$dise_vec
#' dise_fact_na <- dise_na$dise
#'
#' # FI estimator
#' rho_out <- rho_mlogit(dise_fact_na ~ CA125 + CA153 + Age, data = EOC,
#'                       test = TRUE)
#' if (requireNamespace("webshot2", quietly = TRUE)) {
#'    rocs("fi", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'         veri_stat = EOC$V, rho_est = rho_out, ncp = 30)
#' }
#'
#' # Plot ROC surface and add ellipsoid confidence region
#' if (requireNamespace("webshot2", quietly = TRUE)) {
#'    rocs("fi", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'         veri_stat = EOC$V, rho_est = rho_out, ncp = 30,
#'         ellipsoid = TRUE, cpst = c(-0.56, 2.31))
#' }
#'
#' # MSI estimator
#' if (requireNamespace("webshot2", quietly = TRUE)) {
#'    rocs("msi", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'         veri_stat = EOC$V, rho_est = rho_out, ncp = 30,
#'         ellipsoid = TRUE, cpst = c(-0.56, 2.31))
#' }
#'
#' # IPW estimator
#' pi_out <- psglm(V ~ CA125 + CA153 + Age, data = EOC, test = TRUE)
#' if (requireNamespace("webshot2", quietly = TRUE)) {
#'    rocs("ipw", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'         veri_stat = EOC$V, pi_est = pi_out, ncp = 30,
#'         ellipsoid = TRUE, cpst = c(-0.56, 2.31))
#' }
#'
#' # SPE estimator
#' if (requireNamespace("webshot2", quietly = TRUE)) {
#'    rocs("spe", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'         veri_stat = EOC$V, rho_est = rho_out, ncp = 30,
#'         pi_est = pi_out, ellipsoid = TRUE, cpst = c(-0.56, 2.31))
#' }
#'
#' # NN estimator
#' x_mat <- cbind(EOC$CA125, EOC$CA153, EOC$Age)
#' k_opt <- cv_knn(x_mat = x_mat, dise_vec = dise_vec_na, veri_stat = EOC$V,
#'                 type = "mahala", plot = TRUE)
#' rho_k_opt <- rho_knn(x_mat = x_mat, dise_vec = dise_vec_na,
#'                      veri_stat = EOC$V, k = k_opt, type = "mahala")
#' if (requireNamespace("webshot2", quietly = TRUE)) {
#'    rocs("knn", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'         veri_stat = EOC$V, rho_est = rho_k_opt, ncp = 30,
#'         ellipsoid = TRUE, cpst = c(-0.56, 2.31))
#' }
#'
#' ## Compute TCFs at three cut points
#' cutps <- rbind(c(0, 0.5), c(0, 1), c(0.5, 1))
#' rocs.tcf("spe", diag_test = EOC$CA125, dise_vec = dise_vec_na,
#'          veri_stat = EOC$V, rho_est = rho_out, ncp = 30,
#'          pi_est = pi_out, cps = cutps)
#'}
#'

##
#' @rdname rocs
#' @export
rocs.tcf <- function(method = "full", diag_test, dise_vec, veri_stat = NULL,
                     rho_est = NULL, pi_est = NULL, cps) {
  ## checking the method
  method_temp <- substitute(me, list(me = method))
  ok_method <- c("full", "fi", "msi", "ipw", "spe", "knn")
  if (length(method_temp) > 1)
    stop(gettextf("Please, choose one method from %s",
                  paste(sQuote(ok_method), collapse = ", ")), domain = NA)
  if (!is.character(method_temp)) method_temp <- deparse(method_temp)
  if (!is.element(method_temp, ok_method)) {
    stop(gettextf("the required method should be one of %s",
                  method_temp, paste(sQuote(ok_method), collapse = ", ")),
         domain = NA)
  }
  ## checking the argument diag_test
  if (missing(diag_test)) stop("argument \"diag_test\" is missing \n")
  if (!inherits(diag_test, "numeric") || any(is.na(diag_test)))
    stop("\"diag_test\" must be a numeric vector and not include NA values")
  method_name <- toupper(method_temp)
  ## checking dise_vec
  if (missing(dise_vec)) stop("argument \"dise_vec\" is missing \n")
  if (!inherits(dise_vec, "matrix") || ncol(dise_vec) != 3 ||
      !all(is.element(na.omit(dise_vec), c(0, 1))))
    stop("variable \"dise_vec\" must be a binary matrix with 3 columns")
  if (length(diag_test) != nrow(dise_vec))
    stop(gettextf("arguments imply differing number of observation: %d",
                  length(diag_test)), gettextf(", %d", nrow(dise_vec)),
         domain = NA)
  ##
  if (!is.null(rho_est)) {
    if (!is.element(class(rho_est), c("prob_dise", "prob_dise_knn")))
      stop("\"rho_est\" not a result of rho_mlogit or rho_knn \n")
    if (is.element(class(rho_est), c("prob_dise_knn"))) {
      if (is.null(rho_est$K))
        stop("Please, choose one value of K \n")
    }
  }
  if (!is.null(pi_est)) {
    if (!is.element(class(pi_est), c("prob_veri")))
      stop("\"pi_est\" not a result of psglm \n")
  }
  ## Main body
  tcf_cal <- function(tt, dd, cp, weight, method) {
    thet_hat <- colSums(dd) / weight
    names(thet_hat) <- c()
    beta_func <- function(tt, dd, ctp, weight) {
      bet <- numeric(4)
      bet[1] <- sum((tt >= ctp[1]) * dd[, 1]) / weight
      bet[2] <- sum((tt >= ctp[1]) * dd[, 2]) / weight
      bet[3] <- sum((tt >= ctp[2]) * dd[, 2]) / weight
      bet[4] <- sum((tt >= ctp[2]) * dd[, 3]) / weight
      return(bet)
    }
    if (!is.matrix(cp)) {
      bet_hat <- beta_func(tt = tt, dd = dd, ctp = cp, weight = weight)
      res <- numeric(3)
      res[1] <- 1 - bet_hat[1] / thet_hat[1]
      res[2] <- (bet_hat[2] - bet_hat[3]) / thet_hat[2]
      res[3] <- bet_hat[4] / thet_hat[3]
      attr(res, "beta") <- bet_hat
      attr(res, "theta") <- thet_hat
      attr(res, "cp") <- cp
      attr(res, "name") <- method
      names(res) <- c("TCF1", "TCF2", "TCF3")
      class(res) <- "tcfs"
    } else {
      res <- list()
      bet_hat <- t(apply(cp, 1, beta_func, tt = tt, dd = dd, weight = weight))
      temp <- matrix(0, nrow = nrow(cp), ncol = 3)
      temp[, 1] <- 1 - bet_hat[, 1] / thet_hat[1]
      temp[, 2] <- (bet_hat[, 2] - bet_hat[, 3]) / thet_hat[2]
      temp[, 3] <- bet_hat[, 4] / thet_hat[3]
      colnames(temp) <- c("TCF1", "TCF2", "TCF3")
      for (i in seq_len(nrow(cp))) {
        temp1 <- temp[i, ]
        attr(temp1, "beta") <- bet_hat[i, ]
        attr(temp1, "theta") <- thet_hat
        attr(temp1, "cp") <- cp[i, ]
        attr(temp1, "name") <- method
        res[[i]] <- temp1
        class(res[[i]]) <- "tcfs"
      }
    }
    return(res)
  }
  nsize <- length(diag_test)
  if (method == "full") {
    if (any(is.na(dise_vec))) {
      stop("The", method_name,
           "method can not access the NA values of disease status")
    } else {
      ans <- tcf_cal(diag_test, dise_vec, cps, weight = nsize,
                     method = method_name)
    }
  }
  if (method == "fi") {
    if (is.null(rho_est))
      stop("argument \"rho_est\" is needed for ", method_name, "estimator")
    else
      ans <- tcf_cal(diag_test, rho_est$values, cps, weight = nsize,
                     method = method_name)
  }
  if (method %in% c("msi", "knn")) {
    if (is.null(rho_est) || is.null(veri_stat)) {
      stop("argument \"rho_est\" is needed for ", method_name, "estimator")
    } else {
      dise_vec_temp <- dise_vec
      if (any(is.na(dise_vec))) dise_vec_temp[is.na(dise_vec)] <- 99
      dise_msi <- dise_vec_temp * veri_stat + (1 - veri_stat) * rho_est$values
      ans <- tcf_cal(diag_test, dise_msi, cps, weight = nsize,
                     method = method_name)
    }
  }
  if (method == "ipw") {
    if (is.null(pi_est) || is.null(veri_stat)) {
      stop("argument \"pi_est\" is needed for ", method_name, "estimator")
    } else {
      dise_vec_temp <- dise_vec
      if (any(is.na(dise_vec))) dise_vec_temp[is.na(dise_vec)] <- 99
      dise_ipw <- dise_vec_temp * veri_stat / pi_est$values
      ans <- tcf_cal(diag_test, dise_ipw, cps,
                     weight = sum(veri_stat / pi_est$values),
                     method = method_name)
    }
  }
  if (method == "spe") {
    if (is.null(rho_est) || is.null(veri_stat) || is.null(pi_est)) {
      stop("arguments \"rho_est\" and \"pi_est\" are needed for ",
           method_name, " estimator")
    } else {
      dise_vec_temp <- dise_vec
      if (any(is.na(dise_vec))) dise_vec_temp[is.na(dise_vec)] <- 99
      dise_spe <- dise_vec_temp * veri_stat / pi_est$values -
        (veri_stat / pi_est$values - 1) * rho_est$values
      ans <- tcf_cal(diag_test, dise_spe, cps, weight = nsize,
                     method = method_name)
    }
  }
  return(ans)
}

## rocs function
#' @rdname rocs
#' @import utils
#' @import rgl
#' @import parallel
#' @export
rocs <- function(method = "full", diag_test, dise_vec, veri_stat,
                 rho_est = NULL, pi_est = NULL, ncp = 100, plot = TRUE,
                 ellipsoid = FALSE, cpst = NULL, ci_level = 0.95,
                 surf_col = c("gray40", "green"), boot = FALSE,
                 n_boot = 250, parallel = FALSE,
                 ncpus = ifelse(parallel, detectCores() / 2, NULL), ...) {
  ## checking the method
  method_temp <- substitute(me, list(me = method))
  ok_method <- c("full", "fi", "msi", "ipw", "spe", "knn")
  if (length(method_temp) > 1)
    stop(gettextf("Please, choose one method from %s",
                  paste(sQuote(ok_method), collapse = ", ")), domain = NA)
  if (!is.character(method_temp)) method_temp <- deparse(method_temp)
  if (!is.element(method_temp, ok_method))
    stop(gettextf("the required method \"%s\" should be one of %s", method_temp,
                  paste(sQuote(ok_method), collapse = ", ")),
         domain = NA)
  ## checking the argument diag_test
  if (missing(diag_test)) stop("argument \"diag_test\" is missing \n")
  if (!inherits(diag_test, "numeric") || any(is.na(diag_test)))
    stop("\"diag_test\" must be a numeric vector and not include NA values")
  name_diagnostic <- substitute(diag_test)
  if (!is.character(name_diagnostic))
    name_diagnostic <- deparse(name_diagnostic)
  name_diagnostic <- unlist(strsplit(name_diagnostic, NULL))
  if (any(name_diagnostic %in% c("$"))) {
    id_name <- which(name_diagnostic %in% c("$"))
    name_diagnostic <- paste(
      name_diagnostic[(id_name[1] + 1) : length(name_diagnostic)],
      collapse = "")
  } else {
    name_diagnostic <- paste(name_diagnostic, collapse = "")
  }
  method_name <- toupper(method_temp)
  ## checking dise_vec
  if (missing(dise_vec)) stop("argument \"dise_vec\" is missing \n")
  if (isFALSE("matrix" %in% class(dise_vec)) || ncol(dise_vec) != 3 ||
      !all(is.element(na.omit(dise_vec), c(0, 1))))
    stop("\"dise_vec\" must be a binary matrix with 3 columns")
  if (length(diag_test) != nrow(dise_vec))
    stop(gettextf("arguments imply differing number of observation: %d",
                  length(diag_test)),
         gettextf(", %d", nrow(dise_vec)), domain = NA)
  dise_vec_flag <- any(is.na(dise_vec))
  ##
  if (ellipsoid && is.null(cpst))
    stop("\"cpst\" is needed to build the ellipsoidal CR for TCFs")
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
  ## checking veri_stat
  if (missing(veri_stat)) {
    if (method_temp != "full" || dise_vec_flag)
      stop("\"veri_stat\" is missing")
    cat("Hmm, look likes the full data\n")
    cat("Number of observation:", length(diag_test), "\n")
    cat("The verification status is not available\n")
    cat("You are working on FULL or Complete Case approach\n")
    cat("The diagnostic test:", name_diagnostic, "\n")
    cat("Processing .... \n")
    flush.console()
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
      cat("Hmm, look likes the full data\n")
      cat("Number of observation:", length(diag_test), "\n")
      cat("All subjects underwent the verification process\n")
      cat("You are working on FULL or Complete Case approach\n")
      cat("The diagnostic test:", name_diagnostic, "\n")
      cat("Processing .... \n")
      flush.console()
    } else {
      rv <- mean(veri_stat)
      if (!dise_vec_flag) {
        cat("Warning: There are no NA values in variable dise_vec, while",
            paste(round(rv * 100), "%", sep = ""),
            "of the subjects receive disease verification. \n")
        cat("BE CAREFULL OF YOUR INPUT AND RESULTS \n")
        cat("Number of observation:", length(diag_test), "\n")
        cat("You required estimate ROC surface using", method_name, "approach \n")
        cat("The diagnostic test:", name_diagnostic, "\n")
        cat("Processing .... \n")
        flush.console()
      } else {
        cat("Hmm, look likes the incomplete data\n")
        cat("Number of observation:", length(diag_test), "\n")
        cat(paste(round(rv * 100), "%", sep = ""),
            "of the subjects receive disease verification. \n")
        cat("You required estimate ROC surface using", method_name,
            "approach \n")
        cat("The diagnostic test:", name_diagnostic, "\n")
        if (ellipsoid) {
          cat("The ellipsoidal CR for TCFs are also constructed \n")
          if (!is.matrix(cpst) && !is.vector(cpst))
            stop("The cut points should be a vector or matrix.")
          if (is.vector(cpst)) {
            if (length(cpst) != 2 || cpst[1] > cpst[2])
              stop("The first cut point must be less than the second one.")
          }
          if (is.matrix(cpst)) {
            if (ncol(cpst) != 2 || any(cpst[, 1] > cpst[, 2]))
              stop("the first column must be less than or equal to the second.")
          }
        }
        cat("Processing .... \n")
        flush.console()
      }
    }
  }
  cp <- c(-Inf, seq(min(diag_test), max(diag_test), length.out = ncp - 2), Inf)
  cp1 <- rep(cp, seq(ncp - 1, 0, by = -1))
  cp2 <- c()
  for (i in 1:(ncp - 1)) {
    cp2 <- c(cp2, cp[-c(1:i)])
  }
  cpoint <- cbind(cp1, cp2)
  roc_point <- rocs.tcf(method = method_temp, diag_test = diag_test,
                        dise_vec = dise_vec, veri_stat = veri_stat,
                        rho_est = rho_est, pi_est = pi_est, cps = cpoint)
  roc_point <- matrix(unlist(roc_point), ncol = 3, byrow = TRUE)
  colnames(roc_point) <- c("TCF1", "TCF2", "TCF3")
  rownames(roc_point) <- paste("(", round(cp1, 3)," , " , round(cp2, 3), ")",
                              sep = "")
  ct1 <- numeric(ncp - 1)
  for (i in 1:(ncp - 1)) {
    ct1[i] <- i * ncp - i * (i + 1) / 2
  }
  tcf1 <- matrix(roc_point[ct1, 1], ncp - 1, ncp - 1, byrow = FALSE)
  tcf3 <- matrix(roc_point[1:(ncp - 1), 3], ncp - 1, ncp - 1, byrow = TRUE)
  tcf2 <- matrix(0, nrow = ncp - 1, ncol = ncp - 1)
  tcf2[lower.tri(tcf2, diag = TRUE)] <- roc_point[, 2]
  tcf2 <- t(tcf2)
  res <- list()
  res$vals <- roc_point
  res$cpoint <- cpoint
  res$ncp <- ncp
  if (plot) {
    open3d()
    my_user_matrix <- rbind(c(-0.8370321, -0.5446390, -0.0523976, 0),
                            c(0.1272045, -0.2868422, 0.9494949, 0),
                            c(-0.5321618, 0.7880925, 0.3093767, 0),
                            c(0, 0, 0, 1))
    par3d(windowRect = 50 + c(0, 0, 640, 640), userMatrix = my_user_matrix)
    plot3d(0, 0, 0, type = "n", box = FALSE, xlab = " ", ylab = " ", zlab = " ",
           xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1), axes = FALSE, ...)
    axes3d(edges = c("x--", "y--", "z--"))
    mtext3d("TCF 1", "x--", line = 2, at = 0.35)
    mtext3d("TCF 2", "z--", line = 4, at = 0.55)
    mtext3d("TCF 3", "y--", line = 4, at = 0.15, level = 2)
    if (any(c("full","fi", "msi", "ipw", "spe") %in% method_temp)) {
      title_plot <- paste(toupper(method_temp), "estimator")
    }
    if ("knn" %in% method_temp) {
      title_plot <- paste(rho_est$K, "NN-",
                          switch(rho_est$type, "eucli" = "Euclidean",
                                 "manha" = "Manhattan", "canber" = "Canberra",
                                 "lagran" = "Lagrange",
                                 "mahala" = "Mahalanobis"),
                          sep = "")
    }
    bgplot3d({
      plot.new()
      title(main = title_plot, line = 1)
    })
    surface3d(tcf1, tcf3, tcf2, col = surf_col[1], alpha = 0.5, ...)
  }
  if (ellipsoid) {
    shade_ellips <- function(orgi, sig, lev){
      t1 <- sig[2, 2]
      sig[2, 2] <- sig[3, 3]
      sig[3, 3] <- t1
      t1 <- sig[1, 2]
      sig[1, 2] <- sig[1, 3]
      sig[1, 3] <- t1
      sig[lower.tri(sig)] <- sig[upper.tri(sig)]
      ellips <- ellipse3d(sig, centre = orgi[c(1,3,2)],
                          t = sqrt(qchisq(lev, 3)))
      return(ellips)
    }
    res$cpst <- cpst
    if (is.vector(cpst)) {
      tcf_orgi <- rocs.tcf(method = method_temp, diag_test = diag_test,
                           dise_vec = dise_vec, veri_stat = veri_stat,
                           rho_est = rho_est, pi_est = pi_est, cps = cpst)
      tcf_sig <- asy_cov_tcf(tcf_orgi, diag_test = diag_test,
                             dise_vec = dise_vec, veri_stat = veri_stat,
                             rho_est = rho_est, pi_est = pi_est, boot = boot,
                             n_boot = n_boot, parallel = parallel,
                             ncpus = ncpus)
      sig_test <- rcond(tcf_sig)
      res$tcf <- tcf_orgi
      if (sig_test < .Machine$double.eps) {
        cat("The asymptotic variance-covariance matrix of TCFs at",
            paste("(", cpst[1],", ", cpst[2],")", sep = ""),
            "is not singular!\n")
        cat("The ellipsoidal confidence region is not available!\n")
        if (!boot && method_temp != "full")
          cat("Try again with bootstrap process!\n")
        plot3d(tcf_orgi[1], tcf_orgi[3], tcf_orgi[2], type = "s", col = "red",
               radius = 0.01, add = TRUE, ...)
        res$message <- 0
      } else {
        res$message <- 1
        ellip_tcf <- shade_ellips(tcf_orgi, tcf_sig, ci_level)
        plot3d(ellip_tcf, box = FALSE, col = surf_col[2], alpha = 0.5,
               xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1), xlab = " ",
               ylab = " ", zlab = " ", add = TRUE, ...)
        plot3d(tcf_orgi[1], tcf_orgi[3], tcf_orgi[2], type = "s", col = "red",
               radius = 0.01, add = TRUE, ...)
      }
    }
    if (is.matrix(cpst)) {
      res$tcf <- matrix(0, nrow = nrow(cpst), ncol = 3,
                        dimnames = list(paste("(", round(cpst[, 1], 3), " , " ,
                                              round(cpst[, 2], 3), ")",
                                              sep = ""),
                                        c("TCF1", "TCF2", "TCF3")))
      res$message <- numeric(nrow(cpst))
      for (i in seq_len(nrow(cpst))) {
        tcf_orgi <- rocs.tcf(method = method_temp, diag_test = diag_test,
                             dise_vec = dise_vec, veri_stat = veri_stat,
                             rho_est = rho_est, pi_est = pi_est,
                             cps = cpst[i, ])
        tcf_sig <- asy_cov_tcf(tcf_orgi, diag_test = diag_test,
                               dise_vec = dise_vec, veri_stat = veri_stat,
                               rho_est = rho_est, pi_est = pi_est, boot = boot,
                               n_boot = n_boot, parallel = parallel,
                               ncpus = ncpus)
        sig_test <- rcond(tcf_sig)
        res$tcf[i, ] <- tcf_orgi
        if (sig_test < .Machine$double.eps) {
          cat("The asymptotic variance-covariance matrix of TCFs at ",
              paste("(",cpst[i, 1],", ",cpst[i, 2],")", sep = ""),
              "is not singular!\n")
          cat("The ellipsoidal confidence region is not available!\n")
          if (!boot && method_temp != "full")
            cat("Try again with bootstrap process!\n")
          plot3d(tcf_orgi[1], tcf_orgi[3], tcf_orgi[2], type = "s", col = "red",
                 radius = 0.01, add = TRUE, ...)
          res$message[i] <- 0
        } else {
          res$message[i] <- 1
          ellip_tcf <- shade_ellips(tcf_orgi, tcf_sig, ci_level)
          plot3d(ellip_tcf, box = FALSE, col = surf_col[2], alpha = 0.5,
                 xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1), xlab = " ",
                 ylab = " ", zlab = " ", add = TRUE, ...)
          plot3d(tcf_orgi[1], tcf_orgi[3], tcf_orgi[2], type = "s", col = "red",
                 radius = 0.01, add = TRUE, ...)
        }
      }
    }
  }
  cat("DONE\n")
  cat("===============================================================\n")
  selec_ord <- order(rowSums(roc_point) - 1, decreasing = TRUE)
  cat("Some values of TCFs:\n")
  print(roc_point[selec_ord[1:6], ], 3)
  cat("\n")
  if (ellipsoid) {
    cat("Some information for Ellipsoidal Confidence Region(s):\n")
    cat("Confidence level:", ci_level, "\n")
    if (is.vector(cpst)) {
      cat("TCFs at", paste("(", cpst[1], ", ", cpst[2], ")", sep = ""),"are:\n")
      temp <- tcf_orgi
      attributes(temp) <- NULL
      names(temp) <- c("TCF1", "TCF2", "TCF3")
      print(temp, 3)
    }
    if (is.matrix(cpst)) {
      colnames(res$tcf) <- c("TCF1", "TCF2", "TCF3")
      rownames(res$tcf) <- paste("(",round(cpst[,1], 3)," , ",
                                 round(cpst[,2], 3),")",
                                 sep = "")
      cat("TCFs at", rownames(res$tcf), "are:\n")
      print(res$tcf, 3)
    }
  }
  cat("===============================================================\n")
  invisible(res)
}
