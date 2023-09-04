####========================================================================####
## The R code for some specials functions, which are using in the computation ##
## for the asymptotic covariances of TCF's at a specified cut point.          ##
##                                                                            ##
####========================================================================####

## the derivative of function h(beta, theta, tau)
h_deriv <- function(thet, bet, n_nuis) {
  res <- matrix(0, nrow = (6 + n_nuis), ncol = 3)
  res[, 1] <- c(bet[1] / thet[1] ^ 2, 0, -1 / thet[1], 0, 0, 0, rep(0, n_nuis))
  res[, 2] <- c(0, -(bet[2] - bet[3]) / thet[2] ^ 2, 0, 1 / thet[2],
                -1 / thet[2], 0, rep(0, n_nuis))
  res[, 3] <- c(bet[4] / (1 - thet[1] - thet[2]) ^ 2,
                bet[4] / (1 - thet[1] - thet[2]) ^ 2,
                0, 0, 0, 1 / (1 - thet[1] - thet[2]), rep(0, n_nuis))
  return(t(res))
}

## the first derivatives of rho.
rho_deriv <- function(u_matrix, rho_hat, ref_level = c("1", "2")) {
  # u_matrix: the design matrix is used to fit the disease model.
  del_2 <- -u_matrix * rho_hat[, 1] * rho_hat[, 2]
  ref_level <- match.arg(ref_level)
  if (ref_level == "1") {
    del_1 <- u_matrix * rho_hat[, 1] * (1 - rho_hat[, 1])
    res <- cbind(del_1, del_2)
  } else if (ref_level == "2") {
    del_1 <- u_matrix * rho_hat[, 2] * (1 - rho_hat[, 2])
    res <- cbind(del_2, del_1)
  } else {
    stop("Choose reference level in set {`1`, `2`}.")
  }
  return(res)
}

## the first derivatives of 1/pi.
pi_inv_deriv <- function(u_matrix, pi_hat, pi_cov,
                         model = c("logit", "probit", "threshold")) {
  # u_matrix: the design matrix is used to fit the vefication model.
  model <- match.arg(model)
  res <- switch(EXPR = model,
                logit = -u_matrix * (1 - pi_hat) / pi_hat,
                probit = -u_matrix *
                  as.numeric(dnorm(u_matrix %*% pi_cov)) / pi_hat,
                threshold = -u_matrix / pi_hat ^ 2)
  return(res)
}

## Estimating function for rho_k = Pr(D_k = 1 | T, A).
mlog_est_func <- function(dd, vv, u_matrix, rho_hat) {
  temp <- u_matrix * vv
  b1 <- (dd[, 1] - rho_hat[, 1]) * temp
  b2 <- (dd[, 2] - rho_hat[, 2]) * temp
  return(cbind(b1, b2))
}

## Estimating function for pi = Pr(V = 1| T, A).
pi_est_func <- function(vv, u_matrix, pi_hat, pi_cov,
                        model = c("logit", "probit", "threshold")) {
  model <- match.arg(model)
  res <- switch(EXPR = model,
                logit = (vv - pi_hat) * u_matrix,
                probit = u_matrix * as.numeric(dnorm(u_matrix %*% pi_cov)) *
                  (vv - pi_hat) / (pi_hat * (1 - pi_hat)),
                threshold = (vv / pi_hat - 1) * u_matrix / (1 - pi_hat))
  return(res)
}

## Estimating function of Imputation Estimator (FI and MSI).
est_func_ie <- function(dd, vv, tt, u_matrix, cc, cov_rho, rho_hat,
                        thet, bet, m) {
  n_par <- 6 + length(cov_rho)
  res <- matrix(0, nrow = nrow(u_matrix), ncol = n_par)
  res[, 1] <- vv * (m * dd[, 1] - thet[1] + (1 - m) * rho_hat[, 1]) +
    (1 - vv) * (rho_hat[, 1] - thet[1])
  res[, 2] <- vv * (m * dd[, 2] - thet[2] + (1 - m) * rho_hat[, 2]) +
    (1 - vv) * (rho_hat[, 2] - thet[2])
  res[, 3] <- vv * ((tt >= cc[1]) * (m * dd[, 1] + (1 - m) * rho_hat[, 1]) -
                      bet[1]) + (1 - vv) * ((tt >= cc[1]) * rho_hat[, 1] -
                                              bet[1])
  res[, 4] <- vv * ((tt >= cc[1]) * (m * dd[, 2] + (1 - m) * rho_hat[, 2]) -
                      bet[2]) + (1 - vv) * ((tt >= cc[1]) * rho_hat[, 2] -
                                              bet[2])
  res[, 5] <- vv * ((tt >= cc[2]) * (m * dd[, 2] + (1 - m) * rho_hat[, 2]) -
                      bet[3]) + (1 - vv) * ((tt >= cc[2]) * rho_hat[, 2] -
                                              bet[3])
  res[, 6] <- vv * ((tt >= cc[2]) * (m * dd[, 3] + (1 - m) * rho_hat[, 3]) -
                      bet[4]) + (1 - vv) * ((tt >= cc[2]) * rho_hat[, 3] -
                                              bet[4])
  res[, 7:n_par] <- mlog_est_func(dd, vv, u_matrix, rho_hat)
  return(res)
}

## Derivative of Estimating function of Imputation Estimator (FI and MSI).
est_func_ie_deriv <- function(vv, tt, u_matrix, cc, cov_rho, rho_hat, rho_hes,
                              thet, bet, m) {
  der_rho_1 <- rho_deriv(u_matrix, rho_hat, ref_level = "1")
  der_rho_2 <- rho_deriv(u_matrix, rho_hat, ref_level = "2")
  der_rho_3 <- -(der_rho_1 + der_rho_2)
  temp1 <- (1 - m * vv) * der_rho_1
  temp2 <- (1 - m * vv) * der_rho_2
  a <- rbind(colSums(temp1), colSums(temp2))
  b1 <- colSums((tt >= cc[1]) * temp1)
  b2 <- colSums((tt >= cc[1]) * temp2)
  b3 <- colSums((tt >= cc[2]) * temp2)
  b4 <- colSums((tt >= cc[2]) * (1 - m * vv) * der_rho_3)
  b <- rbind(b1, b2, b3, b4)
  res_left <- rbind(diag(-length(tt), 6, 6),
                    matrix(0, nrow = length(cov_rho), ncol = 6))
  res_right <- rbind(a, b, -rho_hes)
  res <- cbind(res_left, res_right)
  return(res)
}

## Estimating function of IPW Estimator
est_func_ipw <- function(dd, vv, tt, u_matrix, cc, cov_pi, pi_hat, thet, bet,
                         model) {
  n_par <- 6 + length(cov_pi)
  res <- matrix(0, nrow = length(tt), ncol = n_par)
  res[, 1] <- vv * (dd[, 1] - thet[1]) / pi_hat
  res[, 2] <- vv * (dd[, 2] - thet[2]) / pi_hat
  res[, 3] <- vv * ((tt >= cc[1]) * dd[, 1] - bet[1]) / pi_hat
  res[, 4] <- vv * ((tt >= cc[1]) * dd[, 2] - bet[2]) / pi_hat
  res[, 5] <- vv * ((tt >= cc[2]) * dd[, 2] - bet[3]) / pi_hat
  res[, 6] <- vv * ((tt >= cc[2]) * dd[, 3] - bet[4]) / pi_hat
  res[, 7:n_par] <- pi_est_func(vv, u_matrix, pi_hat, cov_pi, model)
  return(res)
}

## Derivative of Estimating function of IPW Estimator
est_func_ipw_deriv <- function(dd, vv, tt, u_matrix, cc, cov_pi, pi_hat, pi_hes,
                               thet, bet, model) {
  der_pi_inv <- pi_inv_deriv(u_matrix, pi_hat, cov_pi, model)
  term <- vv * der_pi_inv
  a1 <- colSums(term * (dd[, 1] - thet[1]))
  a2 <- colSums(term * (dd[, 2] - thet[2]))
  b11 <- colSums(term * ((tt >= cc[1]) * dd[, 1] - bet[1]))
  b12 <- colSums(term * ((tt >= cc[1]) * dd[, 2] - bet[2]))
  b22 <- colSums(term * ((tt >= cc[2]) * dd[, 2] - bet[3]))
  b23 <- colSums(term * ((tt >= cc[2]) * dd[, 3] - bet[4]))
  res_left <- rbind(diag(-sum(vv / pi_hat), 6, 6),
                    matrix(0, nrow = length(cov_pi), ncol = 6))
  res_right <- rbind(a1, a2, b11, b12, b22, b23, -pi_hes)
  res <- cbind(res_left, res_right)
  return(res)
}

## Estimating function of SPE Estimator
est_func_spe <- function(dd, vv, tt, u_matrix_rho, u_matrix_pi, cc, cov_rho,
                         rho_hat, cov_pi, pi_hat, thet, bet, model) {
  n_par <- 6 + length(cov_rho) + length(cov_pi)
  res <- matrix(0, nrow = length(tt), ncol = n_par)
  res[, 1] <- vv * (dd[, 1] - thet[1]) / pi_hat -
    (rho_hat[, 1] - thet[1]) * (vv - pi_hat) / pi_hat
  res[, 2] <- vv * (dd[, 2] - thet[2]) / pi_hat -
    (rho_hat[, 2] - thet[2]) * (vv - pi_hat) / pi_hat
  res[, 3] <- vv * ((tt >= cc[1]) * dd[, 1] - bet[1]) / pi_hat -
    (vv - pi_hat) * ((tt >= cc[1]) * rho_hat[, 1] - bet[1]) / pi_hat
  res[, 4] <- vv * ((tt >= cc[1]) * dd[, 2] - bet[2]) / pi_hat -
    (vv - pi_hat) * ((tt >= cc[1]) * rho_hat[, 2] - bet[2]) / pi_hat
  res[, 5] <- vv * ((tt >= cc[2]) * dd[, 2] - bet[3]) / pi_hat -
    (vv - pi_hat) * ((tt >= cc[2]) * rho_hat[, 2] - bet[3]) / pi_hat
  res[, 6] <- vv * ((tt >= cc[2]) * dd[, 3] - bet[4]) / pi_hat -
    (vv - pi_hat) * ((tt >= cc[2]) * rho_hat[, 3] - bet[4]) / pi_hat
  res[, 7:(6 + length(cov_rho))] <- mlog_est_func(dd, vv, u_matrix_rho, rho_hat)
  res[, (7 + length(cov_rho)):n_par] <- pi_est_func(vv, u_matrix_pi, pi_hat,
                                                    cov_pi, model)
  return(res)
}

## Derivative of Estimating function of SPE Estimator
est_func_spe_deriv <- function(dd, vv, tt, u_matrix_rho, u_matrix_pi, cc,
                               cov_rho, rho_hat, rho_hes, cov_pi, pi_hat,
                               pi_hes, thet, bet, model) {
  der_rho_1 <- rho_deriv(u_matrix_rho, rho_hat, ref_level = "1")
  der_rho_2 <- rho_deriv(u_matrix_rho, rho_hat, ref_level = "2")
  der_rho_3 <- -(der_rho_1 + der_rho_2)
  temp_v <- (1 - vv / pi_hat)
  temp1 <- temp_v * der_rho_1
  temp2 <- temp_v * der_rho_2
  der_pi_inv <- pi_inv_deriv(u_matrix_pi, pi_hat, cov_pi, model)
  temp3 <- vv * der_pi_inv
  h <- rbind(colSums(temp1), colSums(temp2))
  g1 <- colSums((tt >= cc[1]) * temp1)
  g2 <- colSums((tt >= cc[1]) * temp2)
  g3 <- colSums((tt >= cc[2]) * temp2)
  g4 <- colSums((tt >= cc[2]) * temp_v * der_rho_3)
  g <- rbind(g1, g2, g3, g4)
  d1 <- colSums(temp3 * (rho_hat[, 1] - dd[, 1]))
  d2 <- colSums(temp3 * (rho_hat[, 2] - dd[, 2]))
  e11 <- colSums(temp3 * ((tt >= cc[1]) * (rho_hat[, 1] - dd[, 1])))
  e12 <- colSums(temp3 * ((tt >= cc[1]) * (rho_hat[, 2] - dd[, 2])))
  e22 <- colSums(temp3 * ((tt >= cc[2]) * (rho_hat[, 2] - dd[, 2])))
  e23 <- colSums(temp3 * ((tt >= cc[2]) * (rho_hat[, 3] - dd[, 3])))
  res_left <- rbind(diag(-length(tt), 6, 6),
                    matrix(0, nrow = length(cov_rho) + length(cov_pi),
                           ncol = 6))
  res_right_rho <- rbind(h, g, -rho_hes,
                         diag(0, nrow = length(cov_pi), ncol = length(cov_rho)))
  res_right_pi <- rbind(d1, d2, e11, e12, e22, e23,
                        diag(0, nrow = length(cov_rho), ncol = length(cov_pi)),
                        -pi_hes)
  res <- cbind(res_left, res_right_rho, res_right_pi)
  return(res)
}
