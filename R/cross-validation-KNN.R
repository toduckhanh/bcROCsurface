####==========================================================================####
## The R code for the cross-validation (KNN estimator).                         ##
##                                                                              ##
####==========================================================================####
#' @title Cross-validation for K nearest-neighbor regression
#'
#' @description This function calculates the estimated cross-validation prediction error for K nearest-neighbor regression and returns a good choice for K.
#'
#' @param X  a numeric design matrix, which used in \code{\link{rhoKNN}} to estimate probabilities of the disease status.
#' @param Dvec  a n * 3  binary matrix with three columns, corresponding to the three classes of the disease status. In row i, 1 indicates, the i-th subject belongs to class j, with j = 1, 2, 3. A row of \code{NA} values indicates a non-verified subject.
#' @param V  a binary vector containing the verification status (1 verified, 0 not verified).
#' @param k.list  a list of candidate values for finding the optimal value of K. If \code{NULL}(the default) then a set of \eqn{\{1, 2, ..., n.ver\}} is employed. Here, \eqn{n.ver} is the number of verified subjects.
#' @param type  a type of distance, see \code{\link{rhoKNN}} for more details.
#' @param plot  if \code{TRUE}, a plot of cross-validation prediction error is produced.
#'
#' @details The data are divided into 2 groups, verified and non verified. In the verified group, the discrepancy between the true disease status and the KNN estimates of the probabilities of the disease status is computed by varying \code{k} from 1 to the number of verification subjects, see To Duc et al (2016). The optimal value of \code{k} is the point that corresponds to the smallest value of the discrepancy.
#'
#' @return A good choice for K is returned.
#'
#' @references
#'
#' To Duc, K., Chiogna, M., Adimari, G. (2016): Nonparametric Estimation of ROC Surfaces Under Verification Bias. \url{https://arxiv.org/abs/1604.04656v1}. Submitted.
#'
#' @examples
#' data(EOC)
#' XX <- cbind(EOC$CA125, EOC$CA153, EOC$Age)
#' Dna <- preDATA(EOC$D, EOC$CA125)
#' Dvec.na <- Dna$Dvec
#' CVknn(XX, Dvec.na, EOC$V, type = "mahala", plot = TRUE)
#'
#' @import graphics
#' @export
CVknn <- function(X, Dvec, V, k.list = NULL, type, plot = FALSE){
  x.ver <- X[V == 1, ]
  y.ver <- Dvec[V == 1, ]
  n.ver <- nrow(x.ver)
  if(!is.null(k.list)){
    if(length(k.list) >= n.ver) stop("k needs to less than n.ver - 1 !! \n")
    else kk <- k.list
  }
  else kk <- seq(1, ceiling(n.ver/2))
  nkk <- length(kk)
  cv_cri <- numeric(nkk)
  rho.ver <- rhoKNN(x.ver, y.ver, rep(1, n.ver), k = kk, type)
  for(i in 1:nkk){
    cv_cri[i] <- sum(colSums(abs(y.ver[, -3] - rho.ver[[i]]$values[, -3])))
  }
  cv_cri <- cv_cri/(n.ver*(ncol(Dvec) - 1))
  ans <- which.min(cv_cri)
  name.type <- switch(type,
                      eucli = "Euclidean", manha = "Manhattan",
                      canber = "Canberra", lagran = "Lagrange",
                      mahala = "Mahalanobis")
  if(plot){
    plot(1/kk, cv_cri, type = "b", xlab = "1/K", ylab = "The discrepancy",
         main = paste("The CV of KNN regression with", name.type, "Distance"))
    points(1/ans, min(cv_cri), pch = 16)
    grid()
    legend("topright", "The good choice for K", pch = 16, cex = 0.7)
  }
  return(ans)
}

