####========================================================================####
## This file consists of R function pre_disease.                              ##
## This is used to make a suitable version of disease status,                 ##
## which version is passed to the functions to                                ##
## obtain bias-corrected estimates of ROC surface and VUS                     ##
##                                                                            ##
####========================================================================####
##
#' @title Preparing monotone ordered disease classes
#'
#' @description \code{pre_data} is used to check and make a suitable monotone increasing ordering of the disease classes. In addition, this function also creates a binary matrix format of disease status to pass into the functions of \code{bcROCsurface} package.
#'
#' @param dise  a numeric vector/factor containing the disease status.
#' @param diag_test  a numeric vector containing the diagnostic test values. \code{NA} values are not admitted.
#' @param plot  if TRUE (the default) then a boxplot of diagnostic test based on three ordered groups is produced.
#'
#' @details The ROC surface analysis implemented in the package is coherent when the ordering of the diagnostic classes is monotone increasing. That is, for a diagnostic test \eqn{T} and three disease classes, 1, 2 and 3, the monotone increasing ordering of interest is \eqn{T_1 < T_2 < T_3}. Here, \eqn{T_1}, \eqn{T_2} and \eqn{T_3} are the measurements of diagnostic test \eqn{T} corresponding to class 1, 2 and 3, respectively. Note that, if an umbrella or tree ordering is of interest, then the results of ROC surface analysis is not reliable.
#'
#' In order to find out the monotone ordering, we compute the medians of \eqn{T_1}, \eqn{T_2} and \eqn{T_3}, and then sort the three medians in ascending order. After that, the three disease classes are reordered corresponding to the order of medians.
#'
#' To be used in the functions of package \code{bcROCsurface}, the vector of disease status must be presented as a n * 3 binary matrix with the three columns, corresponding to the three classes.
#'
#' With real data, the application of this function is the first step in the use of ROC analysis. Note that, if the user is sure that the disease classes follow a monotone increasing ordering and the disease matrix is available, then the use of \code{pre_data} is not necessary.
#'
#' @return This function returns a list containting a factor \code{dise} of ordered disease status and a binary matrix \code{dise_vec} of the disease status and a vector \code{order} containing the sequence of class labels.
#'
#' @examples
#' data(EOC)
#' dise_full <- pre_data(EOC$D.full, EOC$CA125)
#'
#'
#' @export
pre_data <- function(dise, diag_test, plot = TRUE) {
  if (!is.factor(dise)) dise <- as.factor(dise)
  if (length(levels(dise)) != 3) stop("\"dise\" do not have three classes!")
  if (length(diag_test) != length(dise))
    stop(gettextf("arguments imply differing number of observation: %d",
                  length(diag_test)), gettextf(", %d", length(dise)),
         domain = NA)
  classes <- sort(as.character(unique(na.omit(dise))))
  dise_temp <- as.numeric(factor(dise, levels = classes, labels = c(1, 2, 3)))
  if (any(is.na(dise_temp))) cat("There are missing disease status.\n")
  diag_test1 <- na.omit(diag_test[dise_temp == 1])
  diag_test2 <- na.omit(diag_test[dise_temp == 2])
  diag_test3 <- na.omit(diag_test[dise_temp == 3])
  mean_temp <- c(mean(diag_test1), mean(diag_test2), mean(diag_test3))
  median_temp <- c(median(diag_test1), median(diag_test2), median(diag_test3))
  id_median <- order(median_temp)
  cat("The sample means of diagostic test based on three classes.\n")
  for (i in 1:3) {
    cat("(", i, ")", classes[i], ":", round(mean_temp[i], 3), "\n")
  }
  cat("The sample median of diagostic test based on three classes.\n")
  for (i in 1:3) {
    cat("(", i, ")", classes[i], ":", round(median_temp[i], 3), "\n")
  }
  cat("The ordering based on median: ")
  cat(paste(classes[id_median], c(rep("<", 2), ""), sep = " "), "\n")
  rel_diff <- function(x, y) return((x - y) / max(abs(c(x, y))))
  p1 <- rel_diff(median_temp[1], median_temp[2])
  p2 <- rel_diff(median_temp[1], median_temp[3])
  p3 <- rel_diff(median_temp[2], median_temp[3])
  if (abs(p1) <= 0.01) {
    cat("The compararison between the sample medians are employed!\n")
    cat("Warning:\n")
    cat(" (+) The ROC surface analysis may be not reliable.\n")
    cat(" (+) The true ordering may be umbrella! In fact, ")
    if (p3 > 0) cat(paste(classes[c(1, 3, 2)], c(">", "<", ""),
                          sep = " "), "\n")
    if (p3 < 0) cat(paste(classes[c(1, 3, 2)], c("<", ">", ""),
                          sep = " "), "\n")
  }
  if (abs(p2) <= 0.01) {
    cat("The compararison between the sample medians are employed!\n")
    cat("Warning:\n")
    cat(" (+) The ROC surface analysis may be not reliable.\n")
    cat(" (+) The true ordering may be umbrella! In fact, ")
    if (p1 > 0) cat(paste(classes[c(1, 2, 3)], c(">", "<", ""),
                          sep = " "), "\n")
    if (p1 < 0) cat(paste(classes[c(1, 2, 3)], c("<", ">", ""),
                          sep = " "), "\n")
  }
  if (abs(p3) <= 0.01) {
    cat("The compararison between the sample medians are employed!\n")
    cat("Warning:\n")
    cat(" (+) The ROC surface analysis may be not reliable.\n")
    cat(" (+) The true ordering may be umbrella! In fact, ")
    if (p2 > 0) cat(paste(classes[c(2, 1, 3)], c(">", "<", ""),
                          sep = " "), "\n")
    if (p2 < 0) cat(paste(classes[c(2, 1, 3)], c("<", ">", ""),
                          sep = " "), "\n")
  }
  d <- numeric(length(dise_temp))
  d[dise_temp == 1] <- which(id_median == 1)
  d[dise_temp == 2] <- which(id_median == 2)
  d[dise_temp == 3] <- which(id_median == 3)
  d[is.na(dise_temp)] <- NA
  d1 <- d2 <- d3 <- rep(0, length(d))
  d1[d == 1] <- 1
  d2[d == 2] <- 1
  d3[d == 3] <- 1
  if (any(is.na(d))) {
    d1[is.na(d)] <- d2[is.na(d)] <- d3[is.na(d)] <- NA
  }
  dise_vec <- cbind(d1, d2, d3)
  colnames(dise_vec) <- c("D1", "D2", "D3")
  if (plot) {
    boxplot(diag_test ~ d, names = classes[id_median], na.action = na.omit,
            col = c("gray", "cornflowerblue", "red"), ylab = "Diagnostic Test",
            xlab = "Three Ordered Groups")
  }
  res <- list(dise = as.factor(d), dise_vec = dise_vec, order = id_median)
  invisible(res)
}
