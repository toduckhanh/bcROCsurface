####==========================================================================####
## The R code for the verification models.                                      ##
## The suggestion regression models are logit, probit and threshold. The        ##
## default model is logit.                                                      ##
####==========================================================================####
#' @title Fitting verification models
#'
#' @description \code{psglm} is used to fit generalized linear models to the verification process. This function requires a symbolic formula of the linear predictor, and a specified regression model.
#'
#' @param formula  an object of class "formula": a symbolic description of the model to be fitted.
#' @param data  an optional data frame containing the variables in the model.
#' @param model  a specified model to be used in the fitting. The suggestion regression models are logit, probit and threshold. If \code{model} is ignored, then \code{psglm} use a default model as logit.
#' @param test  a logical value indicating whether p-values of the regression coefficients should be returned.
#' @param trace switch for tracing estimation process. Default \code{TRUE}.
#' @param ... optional arguments to be passed to \code{glm}.
#'
#' @details \code{psglm} estimates the verification probabilities of the patients. The suggestion model is designed as a list containing: logit, probit and threshold.
#'
#' @return \code{psglm} returns a list containing the following components:
#'  \item{coeff}{a vector of estimated coefficients.}
#'  \item{values}{fitted values of the model.}
#'  \item{Hess}{the Hessian of the measure of fit at the estimated coefficients.}
#'  \item{x}{a design model matrix.}
#'  \item{formula}{the formula supplied.}
#'  \item{model}{the model object used.}
#'
#' @seealso \code{\link{glm}}
#'
#' @examples
#' data(EOC)
#' out <- psglm(V ~ CA125 + CA153 + Age, data = EOC, test = TRUE)
#'
#'
#' @import stats
#' @export
psglm <- function(formula, data, model = "logit", test = FALSE,
                  trace = TRUE, ...) {
  sugg_model <- c("logit", "probit", "threshold")
  model_temp <- substitute(mod, list(mod = model))
  if (!is.character(model_temp)) model_temp <- deparse(model_temp)
  if (missing(data)) {
    data <- environment(formula)
    cat("Warning: the data input is missing, the global variables are used!\n")
  }
  if (model_temp %in% sugg_model) {
    if (model == "threshold") {
      md_temp <- glm(formula, data = data, family = gaussian, x = TRUE, ...)
    } else {
      md_temp <- glm(formula, data = data, family = binomial(link = model),
                     x = TRUE, ...)
    }
    if (trace) {
      cat("Fitting the verification model by using", deparse(model),
          "regression.\n")
      cat("FORMULAR:", deparse(update.formula(formula, Verification ~ .)), "\n")
      cat("\n")
    }
    res_coef <- coef(md_temp)
    res_est <- predict(md_temp, type = "response")
    res_hess <- solve(summary(md_temp)$cov.scaled)
    if (test) {
      cat("=================================================================\n")
      cat("The p-value calculation for the regression coefficients:\n")
      print(summary(md_temp)$coefficients[, c(3, 4)])
      cat("=================================================================\n")
    }
    fit <- list(coeff = res_coef, values = res_est, hess = res_hess,
                X = md_temp$x, formula = formula, model = model)
    class(fit) <- "prob_veri"
  } else {
    stop(gettextf("model \"%s\" must be one of %s", model_temp,
                  paste(sQuote(sugg_model), collapse = ", ")),
         domain = NA)
  }
  invisible(fit)
}
