####========================================================================####
## The R code for the cross-validation (KNN estimator).                       ##
##                                                                            ##
####========================================================================####
#' @title Cross-validation for K nearest-neighbor regression
#'
#' @description This function calculates the estimated cross-validation prediction error for K nearest-neighbor regression and returns a suitable choice for K.
#'
#' @param x_mat  a numeric design matrix, which used in \code{\link{rho_knn}} to estimate probabilities of the disease status.
#' @param dise_vec  a n * 3  binary matrix with three columns, corresponding to the three classes of the disease status. In row i, 1 in column j indicates that the i-th subject belongs to class j, with j = 1, 2, 3. A row of \code{NA} values indicates a non-verified subject.
#' @param veri_stat  a binary vector containing the verification status (1 verified, 0 not verified).
#' @param k_list  a list of candidate values for K. If \code{NULL}(the default), the set \eqn{\{1, 2, ..., n.ver\}}{{1, 2, ..., n.ver}} is employed, where, \eqn{n.ver} is the number of verified subjects.
#' @param type  a type of distance, see \code{\link{rho_knn}} for more details. Default \code{"eucli"}.
#' @param plot  if \code{TRUE}, a plot of cross-validation prediction error is produced.
#'
#' @details Data are divided into two groups, the first contains the data corresponding to veri_stat = 1, whereas the second contains the data corresponding to veri_stat = 0. In the first group, the discrepancy between the true disease status and the KNN estimates of the probabilities of the disease status is computed by varying \code{k} from 1 to the number of verification subjects, see To Duc et al. (2020). The optimal value of \code{k} is the value that corresponds to the smallest value of the discrepancy.
#'
#' @return A suitable choice for k is returned.
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
#' cv_knn(x_mat, dise_vec_na, EOC$V, type = "mahala", plot = TRUE)
#'
#' @import graphics
#' @export
cv_knn <- function(x_mat, dise_vec, veri_stat, k_list = NULL, type = "eucli",
                   plot = FALSE) {
  if (!inherits(x_mat, "matrix")) stop("\"x_mat\" not a matrix \n")
  if (!inherits(dise_vec, "matrix") || ncol(dise_vec) != 3 ||
      !all(is.element(na.omit(dise_vec), c(0, 1))))
    stop("variable \"dise_vec\" must be a binary matrix with 3 columns")
  if (nrow(x_mat) != nrow(dise_vec))
    stop(gettextf("arguments imply differing number of observation: %d",
                  nrow(x_mat)), gettextf(", %d", nrow(dise_vec)), domain = NA)
  if (missing(veri_stat)) stop("object \"veri_stat\" is missing \n")
  x_ver <- x_mat[veri_stat == 1, ]
  y_ver <- dise_vec[veri_stat == 1, ]
  n_ver <- nrow(x_ver)
  if (!is.null(k_list)) {
    if (length(k_list) >= n_ver) stop("k needs to less than n_ver - 1 !! \n")
    else kk <- k_list
  } else {
    kk <- seq(1, ceiling(n_ver / 2))
  }
  nkk <- length(kk)
  cv_cri <- numeric(nkk)
  rho_ver <- rho_knn(x_ver, y_ver, rep(1, n_ver), k = kk, type = type)
  for (i in 1:nkk) {
    cv_cri[i] <- sum(colSums(abs(y_ver[, -3] - rho_ver[[i]]$values[, -3])))
  }
  cv_cri <- cv_cri / (n_ver * (ncol(dise_vec) - 1))
  ans <- which.min(cv_cri)
  name_type <- switch(type,
                      eucli = "Euclidean", manha = "Manhattan",
                      canber = "Canberra", lagran = "Lagrange",
                      mahala = "Mahalanobis")
  if (plot) {
    plot(1 / kk, cv_cri, type = "b", xlab = "1/K", ylab = "The discrepancy",
         main = paste("The CV of KNN regression with", name_type, "Distance"))
    points(1 / ans, min(cv_cri), pch = 16)
    grid()
    legend("topright", "The good choice for k", pch = 16, cex = 0.7)
  }
  return(ans)
}
