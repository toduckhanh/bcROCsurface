####========================================================================####
## The R code for the disease models.                                         ##
## The default model is multinomial logistic.                                 ##
####========================================================================####
#' @title Fitting disease models via multinomial logistic models
#'
#' @description \code{rho_mlogit} is used to fit multinomial logistic models to the disease process in the verified subjects.
#'
#' @param formula  an object of class "formula": a symbolic description of the model to be fitted.
#' @param data  an optional data frame containing the variables in the model.
#' @param test  a logical value indicating whether p-values of the regression coefficients should be returned. Default \code{FALSE}.
#' @param maxit  maximum number of iterations. Default 500.
#' @param trace  switch for tracing estimation process. Default \code{FALSE}.
#'
#' @details In the formula, the response must be a result of \code{\link{pre_data}}, a factor with three levels, say 1, 2, 3. These levels correspond to three classes of disease status, e.g., non-dieseased, intermediate, diseased. The last class (class 3) is considered as the reference level in multinomal logistic model. In presence of verification bias, the missing (\code{NA}) values correspond to non verified subjects.
#'
#' @return \code{rho_mlogit} returns a list containing the following components:
#'  \item{coeff}{a vector of estimated coefficients.}
#'  \item{values}{fitted values of the model.}
#'  \item{Hess}{the Hessian of the measure of fit at the estimated coefficients.}
#'  \item{D}{the disease status vector used.}
#'  \item{X}{a design model matrix.}
#'  \item{formula}{the fomular supplied.}
#'
#' @seealso \code{\link[nnet]{multinom}}, \code{\link[nnet]{nnet}}
#'
#' @references
#' To Duc, K., Chiogna, M. and Adimari, G. (2016)
#' Bias-corrected methods for estimating the receiver operating characteristic surface of continuous diagnostic tests.
#' \emph{Electronic Journal of Statistics}, \bold{10}, 3063-3113.
#'
#' @examples
#' data(EOC)
#' dise_na <- pre_data(EOC$D, EOC$CA125)
#' dise_fact_na <- dise_na$D
#' out <- rho_mlogit(dise_fact_na ~ CA125 + CA153 + Age, data = EOC,
#'                   test = TRUE, trace = TRUE)
#'
#' @import nnet
#' @export
rho_mlogit <- function(formula, data, test = FALSE, maxit = 500,
                       trace = FALSE) {
  if (missing(data)) {
    data <- environment(formula)
    cat("Warning: the data input is missing, the global variables are used!\n")
  }
  md <- model.frame(formula, data, na.action = NULL)
  x_design <- model.matrix(formula, md)
  dise2 <- model.response(md)
  dise2 <- relevel(as.factor(dise2), ref = "3")
  formula_new <- update.formula(formula, dise2 ~ .)
  data_temp <- data
  data_temp$dise2 <- dise2
  tem1_out <- multinom(formula_new, data = data_temp, na.action = na.omit,
                       maxit = maxit, trace = trace, Hess = TRUE)
  if (trace) {
    cat("Fitting the disease model by using multinomial logistic model via nnet package.\n")
    cat("FORMULAR:", deparse(update.formula(formula, Disease  ~ .)), "\n")
    cat("\n")
  }
  res_coef <- t(coef(tem1_out))
  colnames(res_coef) <- c("1", "2")
  res_pr <- predict(tem1_out, newdata = data_temp, type = "probs")
  res_pr <- res_pr[, c(2, 3, 1)]
  colnames(res_pr) <- c("1", "2", "3")
  hess <- tem1_out$Hessian
  rownames(hess) <- colnames(hess) <- c()
  if (test) {
    z <- res_coef / t(summary(tem1_out)$standard.errors)
    p <- (1 - pnorm(abs(z), 0, 1)) * 2
    cat("====================================================================\n")
    cat("The p-value calculation for the regression coefficients:\n")
    print(p, 4L)
    cat("====================================================================\n")
  }
  fit <- list(coeff = res_coef, values = res_pr, Hess = hess,
              D = model.response(md), X = x_design, formula = formula)
  class(fit) <- "prob_dise"
  invisible(fit)
}
