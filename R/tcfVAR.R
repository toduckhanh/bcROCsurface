####==========================================================================####
## The R code for the asymptotic covariance matrix of true class fractions.     ##
##                                                                              ##
##                                                                              ##
####==========================================================================####
#'
#' @title Asymptotic variance-covariance estimation for True Class Fractions (TCFs) at the cut point \eqn{(c_1, c_2)}
#'
#' @description
#' \code{asyCovTCF} computes the asymptotic variance-covariance matrix of full data (FULL) and bias-corrected estimates (i.e. full imputation, mean score imputation, inverse probabilites weighted, semiparametric efficient and K nearest neighbor) of TCFs.
#'
#'
#' @param obj_tcf  a result of a call to \code{\link{ROCs.tcf}}.
#' @param T  a numeric vector containing the diagnostic test values. \code{NA} values of \code{T} are not accepted.
#' @param Dvec a n * 3  binary matrix with three columns, corresponding to the three classes of the disease status. In row i, 1 indicates, the i-th subject belongs to class j, with j = 1, 2, 3. A row of \code{NA} values indicates a non-verified subject.
#' @param V  a binary vector containing the verification status (1 verified, 0 not verified).
#' @param rhoEst  a result of a call to \code{\link{rhoMLogit}} of \code{\link{rhoKNN}} to fit the disease model.
#' @param piEst  a result of a call to \code{\link{psglm}} to fit the verification model.
#' @param cp  cut point \eqn{(c_1, c_2)} is used for estimating TCFs.
#' @param BOOT  a logical value. Default = \code{FALSE}. If set to \code{TRUE}, a bootstrap resampling is employed to estimate the asymptotic variance-covariance matrix of bias-corrected TCFs.
#' @param nR  the number of bootstrap replicates, which used for \code{BOOT = TRUE} or FULL estimate. Usually this will be a single positive integer. Default 2500.
#' @param parallel  a logical value. If \code{TRUE}, a parallel computing is employed to the bootstrap resampling process.
#' @param ncpus  number of processes to be used in parallel computing. Default is half of available cores.
#'
#' @details For bias-corrected estimates of TCFs, the asymptotic variance-covariance matrix at a fixed cut point is estimated by using the Delta method. The function \code{asyCovTCF} implements the explicit forms, which presented in To Duc et al. (2016a, 2016b). In addition, the bootstrap procedure is also implemented.
#'
#' For FULL estimate, the asymptotic variance-covariance matrix is computed by bootstrap procedure.
#'
#' @return This function returns an estimated asymptotic variance-covariance matrix of FULL estimate and bias-corrected estimates of TCFs at a fixed cut point.
#'
#' @references
#'
#' To Duc, K., Chiogna, M. and Adimari, G. (2016a)
#' Bias-corrected methods for estimating the receiver operating characteristic surface of continuous diagnostic tests.
#' \emph{Electronic Journal of Statistics}. In press.
#'
#' To Duc, K., Chiogna, M. and Adimari, G. (2016b)
#' Nonparametric Estimation of ROC Surfaces Under Verification Bias.
#' \url{https://arxiv.org/abs/1604.04656v1}. Submitted.
#'
#'
#' @examples
#' data(EOC)
#' attach(EOC)
#'
#' # FULL data estimator
#' Dfull <- preDATA(D.full, CA125)
#' Dvec.full <- Dfull$Dvec
#'
#' full.tcf <- ROCs.tcf("full", T = CA125, Dvec = Dvec.full, cps = c(2, 4))
#' full.var <- asyCovTCF(full.tcf, T = CA125, Dvec = Dvec.full, cp = c(2, 4))
#'
#' # Preparing the missing disease status
#' Dna <- preDATA(D, CA125)
#' Dvec.na <- Dna$Dvec
#'
#' rho.out <- rhoMLogit(Dna$D ~ CA125 + CA153 + Age, data = EOC, test = TRUE)
#'
#' ## FI estimator
#' fi.tcf <- ROCs.tcf("fi", T = CA125, Dvec = Dvec.na, V = V,
#'                    rhoEst = rho.out, cps = c(2,4))
#' fi.var <- asyCovTCF(fi.tcf, T = CA125, Dvec = Dvec.na, V = V,
#'                     rhoEst = rho.out, cp = c(2,4))
#'
#' ## MSI estimator
#' msi.tcf <- ROCs.tcf("msi", T = CA125, Dvec = Dvec.na, V = V,
#'                     rhoEst = rho.out, cps = c(2,4))
#' msi.var <- asyCovTCF(msi.tcf, T = CA125, Dvec = Dvec.na, V = V,
#'                          rhoEst = rho.out, cp = c(2,4))
#'
#' ## IPW estimator
#' pi.out <- psglm(V ~ CA125 + CA153 + Age, data = EOC, test = TRUE)
#'
#' ipw.tcf <- ROCs.tcf("ipw", T = CA125, Dvec = Dvec.na, V = V,
#'                     piEst = pi.out, cps = c(2,4))
#' ipw.var <- asyCovTCF(ipw.tcf, T = CA125, Dvec = Dvec.na, V = V,
#'                           piEst = pi.out, cp = c(2,4))
#'
#' ## SPE estimator
#' spe.tcf <- ROCs.tcf("spe", T = CA125, Dvec = Dvec.na, V = V,
#'                    rhoEst = rho.out, piEst = pi.out, cps = c(2,4))
#' spe.var <- asyCovTCF(spe.tcf, T = CA125, Dvec = Dvec.na, V = V,
#'                      rhoEst = rho.out, piEst = pi.out, cp = c(2,4))
#'
#' ## KNN estimators
#' XX <- cbind(CA125, CA153, Age)
#' rho.1nn <- rhoKNN(XX, Dvec.na, V, k = 1, type = "mahala")
#' knn.tcf <- ROCs.tcf("knn", T = CA125, Dvec = Dvec.na, V = V,
#'                    rhoEst = rho.1nn, cps = c(2,4))
#' knn.var <- asyCovTCF(knn.tcf, T = CA125, Dvec = Dvec.na, V = V,
#'                         rhoEst = rho.1nn, cp = c(2,4))
#' @import parallel
#' @export
asyCovTCF <- function(obj_tcf, T, Dvec, V, rhoEst, piEst, cp, BOOT = FALSE,
                      nR = 2500, parallel = FALSE,
                      ncpus = ifelse(parallel, detectCores()/2, NULL)){
  method <- tolower(attr(obj_tcf, "name"))
  Dvec.temp <- Dvec
  if(any(is.na(Dvec))) Dvec.temp[is.na(Dvec)] <- 99
  tcf.orgi <- obj_tcf
  tcf.thet <- attr(tcf.orgi, "theta")
  tcf.bet <- attr(tcf.orgi, "beta")
  if(method == "full"){
    bst.full <- function(dt, inds, cp){
      dat <- dt[inds, ]
      ROCs.tcf(method = "full", T = dat[,1], Dvec = as.matrix(dat[, c(2:4)]),
               cps = cp)
    }
    data <- data.frame(T, Dvec.temp)
    if(parallel){
      res.bst <- boot(data, bst.full, R = nR, cp = cp, parallel = "snow",
                      ncpus = ncpus)
    }
    else res.bst <- boot(data, bst.full, R = nR, cp = cp)
    ans <- var(res.bst$t)
  }
  else if(method == "knn"){
    if(!BOOT){
      if(is.null(rhoEst$K)) stop("what is the value of k!? \n")
      XX <- rhoEst$X
      rho.temp <- rhoKNN(XX, Dvec, V, k = 2, type = rhoEst$type)$values
      pi.est <- psknn(XX, V, rhoEst$type)
      ans <- asy.Cov.KNN(T, tcf.thet, tcf.bet, cp, pi.est, rho.temp, rhoEst$K)
    }
    else{
      bst.knn <- function(dt, inds, cp, k, type){
        dat <- dt[inds,]
        XX <- as.matrix(dat[, -c(1:5)])
        rho.knn <- rhoKNN(XX, dat[, c(2:4)], dat[, 5], k = k, type = type)
        ROCs.tcf(method = "knn", T = dat[,1], Dvec = as.matrix(dat[, c(2:4)]),
                 V = dat[, 5], rhoEst = rho.knn, cps = cp)
      }
      data <- data.frame(T, Dvec.temp, V, rhoEst$X)
      if(parallel){
        res.bst <- boot(data, bst.knn, R = nR, cp = cp, k = rhoEst$K,
                        type = rhoEst$type, parallel = "snow", ncpus = ncpus)
      }
      else res.bst <- boot(data, bst.knn, R = nR, cp = cp, k = rhoEst$K,
                           type = rhoEst$type)
      ans <- var(res.bst$t)
    }
  }
  else if(method == "fi"){
    if(!BOOT){
      term1 <- EstFuncIE_deriv(V, T, rhoEst$X, cp, rhoEst$coeff, rhoEst$values,
                               rhoEst$Hess, tcf.thet, tcf.bet, m = 0)
      term2 <- EstFuncIE(Dvec.temp, V, T, rhoEst$X, cp, rhoEst$coeff,
                         rhoEst$values, tcf.thet, tcf.bet, m = 0)
      term1_inv <- solve(term1)
      Sig <- term1_inv %*% (t(term2) %*% term2) %*% t(term1_inv)
      hh <- h_deriv(tcf.thet, tcf.bet, length(rhoEst$coeff))
      ans <- hh %*% Sig %*% t(hh)
    }
    else{
      bst.fi <- function(dt, inds, formula, cp){
        dat <- dt[inds,]
        out <- rhoMLogit(formula, data = dat)
        ROCs.tcf(method = "fi", T = dat[,1],
                 Dvec = as.matrix(dat[, c("D1", "D2", "D3")]),
                 rhoEst = out, cps = cp)
      }
      name.var <- all.vars(rhoEst$formula)
      data <- data.frame(T, rhoEst$D, rhoEst$X[,name.var[-1]], Dvec.temp)
      names(data) <- c("T", name.var, "D1", "D2", "D3")
      form <- as.formula(paste(paste(name.var[1],"~"), paste(name.var[-1],
                                                             collapse = " + ")))
      if(parallel){
        res.bst <- boot(data, bst.fi, R = nR, formula = form, cp = cp,
                        parallel = "snow", ncpus = ncpus)
      }
      else res.bst <- boot(data, bst.fi, R = nR, formula = form, cp = cp)
      ans <- var(res.bst$t)
    }
  }
  else if(method == "msi"){
    if(!BOOT){
      term1 <- EstFuncIE_deriv(V, T, rhoEst$X, cp, rhoEst$coeff, rhoEst$values,
                               rhoEst$Hess, tcf.thet, tcf.bet, m = 1)
      term2 <- EstFuncIE(Dvec.temp, V, T, rhoEst$X, cp, rhoEst$coeff,
                         rhoEst$values, tcf.thet, tcf.bet, m = 1)
      term1_inv <- solve(term1)
      Sig <- term1_inv %*% (t(term2) %*% term2) %*% t(term1_inv)
      hh <- h_deriv(tcf.thet, tcf.bet, length(rhoEst$coeff))
      ans <- hh %*% Sig %*% t(hh)
    }
    else{
      bst.msi <- function(dt, inds, formula, cp){
        dat <- dt[inds,]
        out <- rhoMLogit(formula, data = dat)
        ROCs.tcf(method = "msi", T = dat[,1],
                 Dvec = as.matrix(dat[, c("D1", "D2", "D3")]), V = dat[, "V"],
                 rhoEst = out, cps = cp)
      }
      name.var <- all.vars(rhoEst$formula)
      data <- data.frame(T, rhoEst$D, rhoEst$X[,name.var[-1]], V, Dvec.temp)
      names(data) <- c("T", name.var, "V", "D1", "D2", "D3")
      form <- as.formula(paste(paste(name.var[1],"~"), paste(name.var[-1],
                                                             collapse = " + ")))
      if(parallel){
        res.bst <- boot(data, bst.msi, R = nR, formula = form, cp = cp,
                        parallel = "snow", ncpus = ncpus)
      }
      else res.bst <- boot(data, bst.msi, R = nR, formula = form, cp = cp)
      ans <- var(res.bst$t)
    }
  }
  else if(method == "ipw"){
    if(!BOOT){
      term1 <- EstFuncIPW_deriv(Dvec.temp, V, T, piEst$X, cp, piEst$coeff,
                                piEst$values, piEst$Hess, tcf.thet, tcf.bet,
                                piEst$model)
      term2 <- EstFuncIPW(Dvec.temp, V, T, piEst$X, cp, piEst$coeff,
                          piEst$values, tcf.thet, tcf.bet, piEst$model)
      term1_inv <- solve(term1)
      Sig <- term1_inv %*% (t(term2) %*% term2) %*% t(term1_inv)
      hh <- h_deriv(tcf.thet, tcf.bet, length(piEst$coeff))
      ans <- hh %*% Sig %*% t(hh)
    }
    else{
      bst.ipw <- function(dt, inds, formula, cp){
        dat <- dt[inds,]
        out <- psglm(formula, data = dat, test = FALSE, trace = FALSE)
        ROCs.tcf(method = "ipw", T = dat[,1],
                 Dvec = as.matrix(dat[, c("D1", "D2", "D3")]), V = dat[, 2],
                 piEst = out, cps = cp)
      }
      name.var <- all.vars(piEst$formula)
      data <- data.frame(T, V, piEst$X[,name.var[-1]], Dvec.temp)
      names(data) <- c("T", name.var, "D1", "D2", "D3")
      form <- as.formula(paste(paste(name.var[1],"~"), paste(name.var[-1],
                                                             collapse = " + ")))
      if(parallel){
        res.bst <- boot(data, bst.ipw, R = nR, formula = form, cp = cp,
                        parallel = "snow", ncpus = ncpus)
      }
      else res.bst <- boot(data, bst.ipw, R = nR, formula = form, cp = cp)
      ans <- var(res.bst$t)
    }
  }
  else{
    if(!BOOT){
      term1 <- EstFuncSPE_deriv(Dvec.temp, V, T, rhoEst$X, piEst$X, cp,
                                rhoEst$coeff, rhoEst$values, rhoEst$Hess,
                                piEst$coeff, piEst$values, piEst$Hess, tcf.thet,
                                tcf.bet, piEst$model)
      term2 <- EstFuncSPE(Dvec.temp, V, T, rhoEst$X, piEst$X, cp, rhoEst$coeff,
                          rhoEst$values, piEst$coeff, piEst$values, tcf.thet,
                          tcf.bet, piEst$model)
      term1_inv <- solve(term1)
      Sig <- term1_inv %*% (t(term2) %*% term2) %*% t(term1_inv)
      hh <- h_deriv(tcf.thet, tcf.bet, length(rhoEst$coeff) + length(piEst$coeff))
      ans <- hh %*% Sig %*% t(hh)
    }
    else{
      bst.spe <- function(dt, inds, formula.rho, formula.pi, cp){
        dat <- dt[inds,]
        out.rho <- rhoMLogit(formula.rho, data = dat)
        out.pi <- psglm(formula.pi, data = dat, test = FALSE, trace = FALSE)
        ROCs.tcf(method = "spe", T = dat[,1],
                 Dvec = as.matrix(dat[, c("D1", "D2", "D3")]), V = dat[, "V"],
                 rhoEst = out.pi, piEst = out.pi, cps = cp)
      }
      name.var.rho <- all.vars(rhoEst$formula)
      name.var.pi <- all.vars(piEst$formula)
      data.rho <- data.frame(rhoEst$X[,name.var.rho[-1]])
      data.pi <- data.frame(piEst$X[,name.var.pi[-1]])
      data <- data.frame(T, V, rhoEst$D, intersect(data.rho, data.pi), Dvec.temp)
      names(data) <- c("T", name.var.pi[1], name.var.rho[1],
                       intersect(name.var.rho, name.var.pi), "D1", "D2", "D3")
      form.rho <- as.formula(paste(paste(name.var.rho[1],"~"),
                                   paste(name.var.rho[-1], collapse = " + ")))
      form.pi <- as.formula(paste(paste(name.var.pi[1],"~"),
                                  paste(name.var.pi[-1], collapse = " + ")))
      if(parallel){
        res.bst <- boot(data, bst.spe, R = nR, formula.rho = form.rho,
                        formula.pi = form.pi, cp = cp, parallel = "snow",
                        ncpus = ncpus)
      }
      else res.bst <- boot(data, bst.spe, R = nR, formula.rho = form.rho,
                           formula.pi = form.pi, cp = cp)
      ans <- var(res.bst$t)
    }
  }
  rownames(ans) <- colnames(ans) <- paste("TCF", c(1:3), sep = "")
  return(ans)
}

