####==========================================================================####
## This file consists of some functions that are related to the computation     ##
## asymptotic variance of VUS                                                   ##
## Date: 05/07/2016																															##
####==========================================================================####
##
#' @title Asymptotic variance estimation for VUS
#'
#' @description \code{asyVarVUS} computes the asymptotic variance of full data (FULL) and bias-corrected estimates (i.e. full imputation, mean score imputation, inverse probabilites weighted, semiparametric efficient and K nearest neighbor) of VUS.
#'
#' @param obj_vus  a result of a call to \code{\link{vus}}.
#' @param T  a numeric vector containing the diagnostic test values. \code{NA} values of \code{T} are not accepted.
#' @param Dvec  a n * 3  binary matrix with three columns, corresponding to the three classes of the disease status. In row i, 1 indicates, the i-th subject belongs to class j, with j = 1, 2, 3. A row of \code{NA} values indicates a non-verified subject.
#' @param V  a binary vector containing the verification status (1 verified, 0 not verified).
#' @param rhoEst  a result of a call to \code{\link{rhoMLogit}} of \code{\link{rhoKNN}} to fit the disease model.
#' @param piEst  a result of a call to \code{\link{psglm}} to fit the verification model.
#' @param BOOT a logical value. Default = \code{FALSE}. If set to \code{TRUE}, a bootstrap resampling is employed to estimate the asymptotic variance of bias-corrected VUS estimates except K nearest neighbor.
#' @param nR  the number of bootstrap replicates, which is used for FULL or KNN estimator, or option \code{BOOT = TRUE}. The defaut is 250.
#' @param parallel  a logical value. If \code{TRUE}, a parallel computing is employed to the bootstrap resampling process.
#' @param ncpus  number of processes to be used in parallel computing. Default is half of available cores.
#'
#'
#' @details
#' For the FULL and K nearest neighbor estimates, a bootstrap resampling process is employed to obtain asymptotic variances.
#'
#' For the full imputation, mean score imputation, inverse probabilites weighted and semiparametric efficient estimates of VUS, the asymptotic variances are computed by using the explicit form. Furthermore, the bootstrap procedure is also considered in case of small sample size.
#'
#' @return \code{asyVarVUS} returns a estimated value of the asymptotic variance.
#'
#' @references
#' To Duc, K., Chiogna, M. and Adimari, G. (2016)
#' Bias-corrected methods for estimating the receiver operating characteristic surface of continuous diagnostic tests.
#' \emph{Electronic Journal of Statistics}. In press.
#'
#'
#' @examples
#' data(EOC)
#' attach(EOC)
#' # Preparing the missing disease status
#' Dna <- preDATA(D, CA125)
#' Dvec.na <- Dna$Dvec
#'
#' rho.out <- rhoMLogit(Dna$D ~ CA125 + CA153 + Age, data = EOC, test = TRUE)
#' pi.out <- psglm(V ~ CA125 + CA153 + Age, data = EOC, test = TRUE)
#' vus.spe <- vus("spe", T = CA125, Dvec = Dvec.na, V = V, rhoEst = rho.out,
#'                piEst = pi.out, ci = FALSE)
#' var.spe <- asyVarVUS(vus.spe, T = CA125, Dvec = Dvec.na, V = V, rhoEst = rho.out,
#'                     piEst = pi.out)
#' \dontrun{
#' var.bst.spe <- asyVarVUS(vus.spe, T = CA125, Dvec = Dvec.na, V = V,
#'                          rhoEst = rho.out, piEst = pi.out, BOOT = TRUE,
#'                          parallel = TRUE)
#' }
#' @importFrom Rcpp evalCpp
#' @import parallel
#' @import boot
#' @export
asyVarVUS <- function(obj_vus, T, Dvec, V = NULL, rhoEst = NULL, piEst = NULL,
                      BOOT = FALSE, nR = 250, parallel = FALSE,
                      ncpus = ifelse(parallel, detectCores()/2, NULL)){
  VUS_obj <- obj_vus$vus.fit
  meth <- tolower(attr(VUS_obj, "name"))
  n <- length(T)
  Dvec.temp <- Dvec
  if(any(is.na(Dvec))) Dvec.temp[is.na(Dvec)] <- 99
  if(meth == "fi"){
    if(!BOOT){
      score <- mlogEstFunc(Dvec.temp, V, rhoEst$X, rhoEst$values)
      hess <- solve(rhoEst$Hess)
      der_rho1 <- t(rho_deriv(rhoEst$X, rhoEst$values, ref_level = "1"))
      der_rho2 <- t(rho_deriv(rhoEst$X, rhoEst$values, ref_level = "2"))
      der_rho3 <- -(der_rho1 + der_rho2)
      Q_Fi <- asyVarVUS_C(T, rhoEst$values, VUS_obj, score, hess, der_rho1,
                          der_rho2, der_rho3)
      ans_var <- n^4*sum(Q_Fi^2)/prod(colSums(rhoEst$values)^2)
    }
    else{
      bst.fi <- function(dt, inds, formula){
        dat <- dt[inds,]
        out <- rhoMLogit(formula, data = dat)
        vusC(dat[,1], out$values)
      }
      name.var <- all.vars(rhoEst$formula)
      data <- data.frame(T, rhoEst$D, rhoEst$X[,name.var[-1]])
      names(data) <- c("T", name.var)
      form <- as.formula(paste(paste(name.var[1],"~"), paste(name.var[-1],
                                                             collapse = " + ")))
      if(parallel){
        res.bst <- boot(data, bst.fi, R = nR, formula = form, parallel = "snow",
                        ncpus = ncpus)
      }
      else res.bst <- boot(data, bst.fi, R = nR, formula = form)
      ans_var <- as.numeric(var(res.bst$t))
    }
  }
  else if(meth == "msi"){
    if(!BOOT){
      D_MSI <- as.matrix(V*Dvec.temp + (1 - V)*rhoEst$values)
      score <- mlogEstFunc(Dvec.temp, V, rhoEst$X, rhoEst$values)
      hess <- solve(rhoEst$Hess)
      der_D_MSI1 <- t((1 - V)*rho_deriv(rhoEst$X, rhoEst$values, ref_level = "1"))
      der_D_MSI2 <- t((1 - V)*rho_deriv(rhoEst$X, rhoEst$values, ref_level = "2"))
      der_D_MSI3 <- -(der_D_MSI1 + der_D_MSI2)
      Q_Msi <- asyVarVUS_C(T, D_MSI, VUS_obj, score, hess, der_D_MSI1,
                           der_D_MSI2, der_D_MSI3)
      ans_var <- n^4*sum(Q_Msi^2)/prod(colSums(D_MSI)^2)
    }
    else{
      bst.msi <- function(dt, inds, formula){
        dat <- dt[inds,]
        out <- rhoMLogit(formula, data = dat)
        Dmsi <- as.matrix(dat[, c("D1", "D2", "D3")]*dat[, "V"] +
                            (1 - dat[, "V"])*out$values)
        vusC(dat[,1], Dmsi)
      }
      name.var <- all.vars(rhoEst$formula)
      data <- data.frame(T, rhoEst$D, rhoEst$X[,name.var[-1]], V, Dvec.temp)
      names(data) <- c("T", name.var, "V", "D1", "D2", "D3")
      form <- as.formula(paste(paste(name.var[1],"~"), paste(name.var[-1],
                                                             collapse = " + ")))
      if(parallel){
        res.bst <- boot(data, bst.msi, R = nR, formula = form, parallel = "snow",
                        ncpus = ncpus)
      }
      else res.bst <- boot(data, bst.msi, R = nR, formula = form)
      ans_var <- as.numeric(var(res.bst$t))
    }
  }
  else if(meth == "ipw"){
    if(!BOOT){
      D_IPW <- as.matrix(V*Dvec.temp/piEst$values)
      score <- piEstFunc(V, piEst$X, piEst$values, piEst$coeff, piEst$model)
      hess <- solve(piEst$Hess)
      der_pi_inv <- pi_inv_deriv(piEst$X, piEst$values, piEst$coeff, piEst$model)
      der_D_IPW1 <- t(V*Dvec.temp[,1]*der_pi_inv)
      der_D_IPW2 <- t(V*Dvec.temp[,2]*der_pi_inv)
      der_D_IPW3 <- t(V*Dvec.temp[,3]*der_pi_inv)
      Q_Ipw <- asyVarVUS_C(T, D_IPW, VUS_obj, score, hess, der_D_IPW1,
                           der_D_IPW2, der_D_IPW3)
      ans_var <- sum(Q_Ipw^2)/prod((colSums(D_IPW)/sum(V/piEst$values))^2)/n^2
    }
    else{
      bst.ipw <- function(dt, inds, formula){
        dat <- dt[inds,]
        out <- psglm(formula, data = dat, test = FALSE, trace = FALSE)
        Dipw <- as.matrix(dat[, c("D1", "D2", "D3")]*dat[, 2]/out$values)
        vusC(dat[,1], Dipw)
      }
      name.var <- all.vars(piEst$formula)
      data <- data.frame(T, V, piEst$X[,name.var[-1]], Dvec.temp)
      names(data) <- c("T", name.var, "D1", "D2", "D3")
      form <- as.formula(paste(paste(name.var[1],"~"), paste(name.var[-1],
                                                             collapse = " + ")))
      if(parallel){
        res.bst <- boot(data, bst.ipw, R = nR, formula = form, parallel = "snow",
                        ncpus = ncpus)
      }
      else res.bst <- boot(data, bst.ipw, R = nR, formula = form)
      ans_var <- as.numeric(var(res.bst$t))
    }
  }
  else if(meth == "spe"){
    if(!BOOT){
      temp1 <- V/piEst$values
      D_SPE <- as.matrix(temp1*Dvec.temp - (temp1 - 1)*rhoEst$values)
      score1 <- mlogEstFunc(Dvec.temp, V, rhoEst$X, rhoEst$values)
      score2 <- piEstFunc(V, piEst$X, piEst$values, piEst$coeff, piEst$model)
      score <- cbind(score1, score2)
      hess1 <- cbind(solve(rhoEst$Hess), matrix(0, nrow = nrow(rhoEst$Hess),
                                                ncol = ncol(piEst$Hess)))
      hess2 <- cbind(matrix(0, nrow = nrow(piEst$Hess), ncol = ncol(rhoEst$Hess)),
                     solve(piEst$Hess))
      hess <- rbind(hess1, hess2)
      der_rho1 <- rho_deriv(rhoEst$X, rhoEst$values, ref_level = "1")
      der_rho2 <- rho_deriv(rhoEst$X, rhoEst$values, ref_level = "2")
      der_pi_inv <- pi_inv_deriv(piEst$X, piEst$values, piEst$coeff, piEst$model)
      der_D_SPE1 <- t(cbind(-(temp1 - 1)*der_rho1,
                            V*(Dvec.temp[,1] - rhoEst$values[,1])*der_pi_inv))
      der_D_SPE2 <- t(cbind(-(temp1 - 1)*der_rho2,
                            V*(Dvec.temp[,2] - rhoEst$values[,2])*der_pi_inv))
      der_D_SPE3 <- -(der_D_SPE1 + der_D_SPE2)
      Q_Spe <- asyVarVUS_C(T, D_SPE, VUS_obj, score, hess, der_D_SPE1,
                           der_D_SPE2, der_D_SPE3)
      ans_var <- n^4*sum(Q_Spe^2)/prod(colSums(D_SPE)^2)
    }
    else{
      bst.spe <- function(dt, inds, formula.rho, formula.pi){
        dat <- dt[inds,]
        out.rho <- rhoMLogit(formula.rho, data = dat)
        out.pi <- psglm(formula.pi, data = dat, test = FALSE, trace = FALSE)
        Dspe <- as.matrix(dat[, c("D1", "D2", "D3")]*dat[, "V"]/out.pi$values -
                            (dat[, "V"]/out.pi$values - 1)*out.rho$values)
        vusC(dat[,1], Dspe)
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
                        formula.pi = form.pi, parallel = "snow", ncpus = ncpus)
      }
      else res.bst <- boot(data, bst.spe, R = nR, formula.rho = form.rho,
                           formula.pi = form.pi)
      ans_var <- as.numeric(var(res.bst$t))
    }
  }
  else if(meth == "knn"){
    bst.knn <- function(dt, inds, k, type){
      dat <- dt[inds,]
      XX <- as.matrix(dat[, -c(1:5)])
      rho.knn <- rhoKNN(XX, dat[, c(2:4)], dat[, 5], k = k, type = type)
      Dknn <- as.matrix(dat[, c(2:4)]*dat[, 5] + (1 - dat[, 5])*rho.knn$values)
      vusC(dat[,1], Dknn)
    }
    data <- data.frame(T, Dvec.temp, V, rhoEst$X)
    if(parallel){
      res.bst <- boot(data, bst.knn, R = nR, k = rhoEst$K, type = rhoEst$type,
                      parallel = "snow", ncpus = ncpus)
    }
    else res.bst <- boot(data, bst.knn, R = nR, k = rhoEst$K, type = rhoEst$type)
    ans_var <- as.numeric(var(res.bst$t))
  }
  else{
    bst.full <- function(dt, inds){
      dat <- dt[inds,]
      vusC(dat[,1], as.matrix(dat[,c(2:4)]))
    }
    data <- data.frame(T, Dvec.temp)
    if(parallel){
      res.bst <- boot(data, bst.full, R = nR, parallel = "snow", ncpus = ncpus)
    }
    else res.bst <- boot(data, bst.full, R = nR)
    ans_var <- as.numeric(var(res.bst$t))
  }
  return(ans_var)
}


