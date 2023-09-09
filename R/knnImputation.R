####==========================================================================####
## The R code for fitting the disease models with K nearest-neighbor regression.##
##                                                                              ##
####==========================================================================####
#' @title K nearest-neighbor (KNN) regression
#'
#' @description  \code{rho_knn} uses the KNN approach to estimate the probabilities of the disease status in case of three categories.
#'
#' @param x_mat  a numeric design matrix.
#' @param dise_vec   a n * 3  binary matrix with three columns, corresponding to the three classes of the disease status. In row i, 1 in column j indicates that the i-th subject belongs to class j, with j = 1, 2, 3. A row of \code{NA} values indicates a non-verified subject.
#' @param veri_stat  a binary vector containing the verification status (1 verified, 0 not verified).
#' @param k  an integer value/vector, which indicates the number of nearest neighbors. It should be less than the number of the verification subjects.
#' @param type  a distance measure.
#' @param trace switch for tracing estimation process. Default \code{FALSE}.
#'
#' @details \code{type} should be selected as one of \code{"eucli"}, \code{"manha"}, \code{"canber"}, \code{"lagran"}, \code{"mahala"} corresponding to Euclidean, Manhattan, Canberra, Lagrange and Mahalanobis distance. In practice, the selection of a suitable distance is typically dictated by features of the data and possible subjective evaluations. For example, if the covariates are heterogeneous with respect to their variances (which is particularly true when the variables are measured on heterogeneous scales), the choice of the Mahalanobis distance may be a good choice.
#'
#' For the number of nearest neighbors, a small value of \code{k}, within the range 1-3, may be a good choice. In general, the choice of \code{k} may depend on the dimension of the feature space, and propose to use cross--validation to find \code{k} in case of high--dimensional covariate. See \code{\link{cv_knn}}.
#'
#'
#' @return \code{rho_knn} returns a list containing the following components:
#'  \item{values}{estimates of the probabilities.}
#'  \item{X}{a design model matrix.}
#'  \item{K}{the number of nearest neighbors.}
#'  \item{type}{the chosen distance.}
#'
#' @references
#' To Duc, K., Chiogna, M. and Adimari, G. (2020)
#' Nonparametric estimation of ROC surfaces in presence of verification bias.
#' \emph{REVSTAT-Statistical Journal}. \bold{18}, 5, 697–720.
#'
#' @examples
#' data(EOC)
#' x_mat <- cbind(EOC$CA125, EOC$CA153, EOC$Age)
#' dise_na <- pre_data(EOC$D, EOC$CA125)
#' dise_vec_na <- dise_na$dise_vec
#'
#' ## Euclidean distance, k = 1
#' out_ecul_1nn <- rho_knn(x_mat, dise_vec_na, EOC$V, k = 1, type = "eucli")
#'
#' ## Manhattan distance, k = 1
#' out_manh_1nn <- rho_knn(x_mat, dise_vec_na, EOC$V, k = 1, type = "manha")
#'
#' ## Canberra distance, k = 3
#' out_canb_1nn <- rho_knn(x_mat, dise_vec_na, EOC$V, k = 3, type = "canber")
#'
#' ## Lagrange distance, k = 3
#' out_lagr_1nn <- rho_knn(x_mat, dise_vec_na, EOC$V, k = 3, type = "lagran")
#'
#' ## Mahalanobis distance, k = c(1,3)
#' out_maha_13nn <- rho_knn(x_mat, dise_vec_na, EOC$V, k = c(1, 3),
#'                          type = "mahala")
#'
#' @export
rho_knn <- function(x_mat, dise_vec, veri_stat, k,
                    type = c("eucli", "manha", "canber", "lagran", "mahala"),
                    trace = FALSE) {
  if (!inherits(x_mat, "matrix")) stop("\"x_mat\" not a matrix \n")
  if (!inherits(dise_vec, "matrix") || ncol(dise_vec) != 3 ||
      !all(is.element(na.omit(dise_vec), c(0, 1))))
    stop("\"dise_vec\" must be a binary matrix with 3 columns")
  if (nrow(x_mat) != nrow(dise_vec))
    stop(gettextf("arguments imply differing number of observation: %d",
                  nrow(x_mat)), gettextf(", %d", nrow(dise_vec)), domain = NA)
  if (missing(veri_stat)) stop("object \"veri_stat\" is missing \n")
  if (missing(k)) stop("object \"k\" is missing \n")
  nk <- length(k)
  if (nk == 0) stop("what is the value of k!? \n")
  is_wholenumber <- function(x) all(x == round(x)) & all(x > 0)
  if (!is_wholenumber(k))
    stop("\"k\" not a positive integer number or vector \n")
  n <- nrow(x_mat)
  s_inv <- solve(cov(x_mat))
  type <- match.arg(type)
  dista <- switch(type,
                  eucli = function(xi, xx, sx_inv) {
                    ss <- colSums((xi - t(xx)) ^ 2)
                    return(sqrt(ss))
                  },
                  manha = function(xi, xx, sx_inv) {
                    ss <- colSums(abs(xi - t(xx)))
                    return(ss)
                  },
                  canber = function(xi, xx, sx_inv) {
                    ss <- colSums(abs(xi - t(xx)) / (abs(xi) + t(abs(xx))))
                    return(ss)
                  },
                  lagran = function(xi, xx, sx_inv) {
                    ss <- apply(abs(xi - t(xx)), 2, max)
                    return(ss)
                  },
                  mahala = function(xi, xx, sx_inv) {
                    ss1 <- (xi - t(xx))
                    ss <- diag(t(ss1) %*% sx_inv %*% ss1)
                    return(sqrt(ss))
                  }
  )
  if (trace) {
    cat("Fitting the disease model by using K nearest-neighbor imputation.\n")
    cat("K:", k, "\n")
    cat("Distance measure:", switch(type,
                                    eucli = "Euclidean",
                                    manha = "Manhattan",
                                    canber = "Canberra",
                                    lagran = "Lagrange",
                                    mahala = "Mahalanobis"), "\n")
    cat("\n")
  }
  dise_vec_veri <- dise_vec[veri_stat == 1, ]
  x_mat_veri <- x_mat[veri_stat == 1, ]
  id_ms <- which(veri_stat == 0)
  id_veri <- which(veri_stat == 1)
  if (nk == 1) {
    res <- matrix(0, nrow = n, ncol = 3)
    for (i in id_ms) {
      dist_tem <- dista(x_mat[i, ], x_mat_veri, s_inv)
      id_tem <- order(dist_tem)
      y_temp <- dise_vec_veri[id_tem, ]
      res[i, ] <- c(mean(y_temp[1:k, 1]), mean(y_temp[1:k, 2]),
                    mean(y_temp[1:k, 3]))
    }
    for (i in id_veri) {
      dist_tem <- dista(x_mat[i, ], x_mat_veri, s_inv)
      id_tem <- order(dist_tem)[-1]
      y_temp <- dise_vec_veri[id_tem, ]
      res[i, ] <- c(mean(y_temp[1:k, 1]), mean(y_temp[1:k, 2]),
                    mean(y_temp[1:k, 3]))
    }
    fit <- list(values = res, X = x_mat, K = k, type = type)
    class(fit) <- "prob_dise_knn"
  } else {
    res <- list()
    for (j in 1:nk) {
      res[[j]] <- matrix(0, nrow = n, ncol = 3)
    }
    for (i in id_ms) {
      dist_tem <- dista(x_mat[i, ], x_mat_veri, s_inv)
      id_tem <- order(dist_tem)
      y_temp <- dise_vec_veri[id_tem, ]
      for (j in 1:nk) {
        res[[j]][i, ] <- c(mean(y_temp[1:k[j], 1]), mean(y_temp[1:k[j], 2]),
                           mean(y_temp[1:k[j], 3]))
      }
    }
    for (i in id_veri) {
      dist_tem <- dista(x_mat[i, ], x_mat_veri, s_inv)
      id_tem <- order(dist_tem)[-1]
      y_temp <- dise_vec_veri[id_tem, ]
      for (j in 1:nk) {
        res[[j]][i, ] <- c(mean(y_temp[1:k[j], 1]), mean(y_temp[1:k[j], 2]),
                           mean(y_temp[1:k[j], 3]))
      }
    }
    fit <- list()
    for (j in 1:nk) {
      fit[[j]] <- list(values = res[[j]], X = x_mat, K = k[j], type = type)
      class(fit[[j]]) <- "prob_dise_knn"
    }
  }
  invisible(fit)
}
