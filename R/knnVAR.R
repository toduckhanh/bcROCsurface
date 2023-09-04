####========================================================================####
## The R code for estimating asymptotic covariance matrix of KNN estimators   ##
##                                                                            ##
##                                                                            ##
####========================================================================####

psknn <- function(x, y,
                  type = c("eucli", "manha", "canber", "lagran", "mahala")) {
  n <- nrow(x)
  s_inv <- solve(cov(x))
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
  res <- numeric(n)
  k1 <- rep(0, n)
  for (i in 1:n) {
    dist_tem <- dista(x[i, ], x, s_inv)
    id_tem <- order(dist_tem)[-1]
    y_temp <- y[id_tem]
    if (y[i] == 0) {
      if (y_temp[1] == 0) {
        k1[i] <- which((1 - y_temp) == 0)[1]
      } else {
        k1[i] <- which((1 - y_temp) == 1)[1]
      }
    } else {
      if (y_temp[1] == 1) {
        k1[i] <- which((1 - y_temp) == 1)[1]
      } else {
        k1[i] <- which((1 - y_temp) == 0)[1]
      }
    }
    res[i] <- mean(y_temp[1:k1[i]])
  }
  return(res)
}

# ### Computing the variances
asy_cov_knn <- function(diag_test, thet, bet, cc, pi_hat, rho_hat, k) {
  n <- length(diag_test)
  res <- matrix(0, 3, 3)
  bet[5] <- bet[2] - bet[3]
  ### Computing the variances
  omega1 <- mean((rho_hat[, 1] * (1 - rho_hat[, 1]) * (1 - pi_hat)) *
                  ((k + 1) / k + (1 - pi_hat) / pi_hat))
  omega2 <- mean((rho_hat[, 2] * (1 - rho_hat[, 2]) * (1 - pi_hat)) *
                  ((k + 1) / k + (1 - pi_hat) / pi_hat))
  omega3 <- mean((rho_hat[, 3] * (1 - rho_hat[, 3]) * (1 - pi_hat)) *
                  ((k + 1) / k + (1 - pi_hat) / pi_hat))
  sig_sq <- (thet * (1 - thet) + c(omega1, omega2, omega3))
  omega11 <- mean(((diag_test >= cc[1]) * rho_hat[, 1] * (1 - rho_hat[, 1]) *
                    (1 - pi_hat)) * ((k + 1) / k + (1 - pi_hat) / pi_hat))
  omega12 <- mean((diag_test >= cc[1]) * (rho_hat[, 2] * (1 - rho_hat[, 2]) *
                    (1 - pi_hat)) * ((k + 1) / k + (1 - pi_hat) / pi_hat))
  omega22 <- mean((diag_test >= cc[2]) * (rho_hat[, 2] * (1 - rho_hat[, 2]) *
                    (1 - pi_hat)) * ((k + 1) / k + (1 - pi_hat) / pi_hat))
  omega23 <- mean((diag_test >= cc[2]) * (rho_hat[, 3] * (1 - rho_hat[, 3]) *
                    (1 - pi_hat)) * ((k + 1) / k + (1 - pi_hat) / pi_hat))
  omega <- omega12 - omega22
  sig_sq_bet <- (bet * (1 - bet) + c(omega11, omega12, omega22, omega23, omega))
  gam <- sig_sq_gam <- numeric(4)
  gam[1] <- thet[1] - bet[1]
  gam[2] <- thet[2] - bet[2]
  gam[3] <- thet[2] - bet[3]
  gam[4] <- thet[3] - bet[4]
  sig_sq_gam[1] <- (gam[1] * (1 - gam[1]) + omega1 - omega11)
  sig_sq_gam[2] <- (gam[2] * (1 - gam[2]) + omega2 - omega12)
  sig_sq_gam[3] <- (gam[3] * (1 - gam[3]) + omega2 - omega22)
  sig_sq_gam[4] <- (gam[4] * (1 - gam[4]) + omega3 - omega23)
  sig111 <- (sig_sq[1] + sig_sq_bet[1] - sig_sq_gam[1]) / 2
  sig212 <- (sig_sq[2] + sig_sq_bet[2] - sig_sq_gam[2]) / 2
  sig222 <- (sig_sq[2] + sig_sq_bet[3] - sig_sq_gam[3]) / 2
  sig323 <- (sig_sq[3] + sig_sq_bet[4] - sig_sq_gam[4]) / 2
  res[1, 1] <- bet[1] ^ 2 * sig_sq[1] / thet[1] ^ 4 +
    sig_sq_bet[1] / thet[1] ^ 2 - 2 * bet[1] * sig111 / thet[1] ^ 3
  res[2, 2] <- bet[5] ^ 2 * sig_sq[2] / thet[2] ^ 4 +
    sig_sq_bet[5] / thet[2] ^ 2 - 2 * bet[5] * (sig212 - sig222) / thet[2] ^ 3
  res[3, 3] <- bet[4] ^ 2 * sig_sq[3] / thet[3] ^ 4 +
    sig_sq_bet[4] / thet[3] ^ 2 - 2 * bet[4] * sig323 / thet[3] ^ 3
  ### Estimating xi_12
  omeg_sig12 <- mean((rho_hat[, 1] * rho_hat[, 2] * (1 - pi_hat)) *
                      ((k + 1) / k + (1 - pi_hat) / pi_hat))
  omeg1112 <- mean((diag_test >= cc[1]) * (rho_hat[, 1] * rho_hat[, 2] *
                      (1 - pi_hat)) * ((k + 1) / k + (1 - pi_hat) / pi_hat))
  omeg1122 <- mean((diag_test >= cc[2]) * (rho_hat[, 1] * rho_hat[, 2] *
                      (1 - pi_hat)) * ((k + 1) / k + (1 - pi_hat) / pi_hat))
  omeg11_12_22 <- omeg1112 - omeg1122
  sig_11_12_22 <- -(bet[1] * bet[5] + omeg11_12_22)
  sig_1_12_22 <- -(thet[1] * bet[5] + omeg11_12_22)
  sig_2_11 <- -(thet[2] * bet[1] + omeg1112)
  sig_thet_12 <- -(thet[1] * thet[2] + omeg_sig12)
  res[1, 2] <- -sig_11_12_22 / (thet[1] * thet[2]) + bet[1] * sig_1_12_22 /
                  (thet[1] ^ 2 * thet[2]) -
                  bet[5] * (bet[1] * sig_thet_12 / thet[1] ^ 2 -
                  sig_2_11 / thet[1]) / thet[2] ^ 2
  ###Estimating xi_13
  omeg1123 <- mean((diag_test >= cc[2]) * (rho_hat[, 1] * rho_hat[, 3] *
                   (1 - pi_hat)) * ((k + 1) / k + (1 - pi_hat) / pi_hat))
  omeg311 <- mean((diag_test >= cc[1]) * (rho_hat[, 1] * rho_hat[, 3] *
                   (1 - pi_hat)) * ((k + 1) / k + (1 - pi_hat) / pi_hat))
  sig_1_23 <- -(thet[1] * bet[4] + omeg1123)
  sig_11_23 <- -(bet[1] * bet[4] + omeg1123)
  sig_3_11 <- -(thet[3] * bet[1] + omeg311)
  res[1, 3] <- (bet[1] * sig_1_23 / thet[1] - sig_11_23) / (thet[1] * thet[3]) +
                bet[4] * (bet[1] * (sig_sq[1] + sig_thet_12) / thet[1] +
                            sig_3_11) / (thet[1] * thet[3] ^ 2)
  ###Estimating xi_23
  omeg312 <- mean((diag_test >= cc[1]) * (rho_hat[, 2] * rho_hat[, 3] *
                  (1 - pi_hat)) * ((k + 1) / k + (1 - pi_hat) / pi_hat))
  omeg322 <- mean((diag_test >= cc[2]) * (rho_hat[, 2] * rho_hat[, 3] *
                  (1 - pi_hat)) * ((k + 1) / k + (1 - pi_hat) / pi_hat))
  omeg3_12_22 <- omeg312 - omeg322
  sig_23_12_22 <- -bet[4] * bet[5]
  sig_2_23 <- -(thet[2] * bet[4] + omeg322)
  sig_3_12_22 <- -(thet[3] * bet[5] + omeg3_12_22)
  res[2, 3] <- (sig_23_12_22 - bet[5] * sig_2_23 / thet[2]) /
    (thet[3] * thet[2]) + bet[4] * (-sig_3_12_22 - bet[5] *
            (sig_sq[2] + sig_thet_12) / thet[2]) / (thet[2] * thet[3] ^ 2)
  res[lower.tri(res)] <- res[upper.tri(res)]
  return(res / n)
}
