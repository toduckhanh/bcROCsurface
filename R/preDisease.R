####==========================================================================####
## This file consists of R function preDisease. This is used to make a suitble  ##
## version of disease status, which version is passed to the functions to       ##
## obtain bias-corrected estimates of ROC surface and VUS                       ##
##                                                                              ##
## Date: 21/07/2016																															##
####==========================================================================####
##
#' @title Preparing monotone ordered disease classes
#'
#' @description \code{preDATA} is used to check and make a suitable monotone increasing ordering of the disease classes. In addition, this function also creates a binary matrix format of disease status to pass into the functions of \code{bcROCsurface} package.
#'
#' @param D  a numeric vector/factor containing the disease status.
#' @param T  a numeric vector containing the diagnostic test values. \code{NA} values are not admitted.
#' @param plot  if TRUE (the default) then a boxplot of diagnostic test T based on three ordered groups is produced.
#'
#' @details The ROC surface and VUS are only performed when the ordering of interest is monotone increasing. That is, for diagnostic test \eqn{T} and three classes 1, 2 and 3 of disease status, the monotone increasing ordering of interest is \eqn{T_1 < T_2 < T_3}. Here, \eqn{T_1}, \eqn{T_2} and \eqn{T_3} are the measurements of diagnostic test \eqn{T} corresponding to class 1, 2 and 3, respectively. Note that, if an umbrella or tree ordering is of interest, then the results of ROC surface analysis is not reliable.
#'
#' In order to find out the ordering of the disease classes, we employ the explorator analysis based on the primary ordering of the disease status (e.g. non-diseased, intermediate, diseased). After that, we re-describe the ordering such that it corresponds to the ordering for diagnostic test values.
#'
#' To pass into the functions of package \code{bcROCsurface}, the vector of disease status must be presented as a n * 3 binary matrix with the three columns, corresponding to three classes of the disease status.
#'
#' For any real data, the application of this function is the first step in the use of ROC analysis. Note that, if the user be sure that the disease classes are monotone increasing ordering and the matrix format is available, then the use of \code{preDATA} can be ignorble.
#'
#' @return This function returns a list containting a vector \code{D} of ordered disease status, a binary matrix \code{Dvec} of the disease status and a vector \code{order} containing the sequence of class labels.
#'
#' @examples
#' data(EOC)
#' attach(EOC)
#' Dfull <- preDATA(D.full, CA125)
#'
#'
#' @export
preDATA <- function(D, T, plot = TRUE){
  if(!is.factor(D)) D <- as.factor(D)
  if (length(levels(D)) != 3) stop("\"D\" do not have three classes!")
  classes <- sort(as.character(unique(na.omit(D))))
  D.temp <- as.numeric(factor(D, levels = classes, labels = c(1, 2, 3)))
  if(any(is.na(D.temp))) cat("There are missing disease status.\n")
  T1 <- na.omit(T[D.temp == 1])
  T2 <- na.omit(T[D.temp == 2])
  T3 <- na.omit(T[D.temp == 3])
  mean.temp <- c(mean(T1), mean(T2), mean(T3))
  median.temp <- c(median(T1), median(T2), median(T3))
  id.median <- order(median.temp)
  cat("The sample means of diagostic test based on three classes.\n")
  for(i in 1:3){
    cat("(",i,")", classes[i], ":", round(mean.temp[i], 3), "\n")
  }
  cat("The sample median of diagostic test based on three classes.\n")
  for(i in 1:3){
    cat("(",i,")", classes[i], ":", round(median.temp[i], 3), "\n")
  }
  cat("The ordering based on median: ")
  cat(paste(classes[id.median], c(rep("<",2), ""), sep = " "), "\n")
  rel.diff <- function(x ,y) return((x - y)/max(abs( c(x,y) )))
  p1 <- rel.diff(median.temp[1], median.temp[2])
  p2 <- rel.diff(median.temp[1], median.temp[3])
  p3 <- rel.diff(median.temp[2], median.temp[3])
  if(abs(p1) <= 0.01){
    cat("The compararison between the sample medians are employed!\n")
    cat("Warning:\n")
    cat(" (+) The ROC surface analysis may be not reliable.\n")
    cat(" (+) The true ordering may be umbrella! In fact, ")
    if(p3 > 0) cat(paste(classes[c(1,3,2)], c(">", "<", ""), sep = " "),"\n")
    if(p3 < 0) cat(paste(classes[c(1,3,2)], c("<", ">", ""), sep = " "),"\n")
  }
  if(abs(p2) <= 0.01){
    cat("The compararison between the sample medians are employed!\n")
    cat("Warning:\n")
    cat(" (+) The ROC surface analysis may be not reliable.\n")
    cat(" (+) The true ordering may be umbrella! In fact, ")
    if(p1 > 0) cat(paste(classes[c(1,2,3)], c(">", "<", ""), sep = " "),"\n")
    if(p1 < 0) cat(paste(classes[c(1,2,3)], c("<", ">", ""), sep = " "),"\n")
  }
  if(abs(p3) <= 0.01){
    cat("The compararison between the sample medians are employed!\n")
    cat("Warning:\n")
    cat(" (+) The ROC surface analysis may be not reliable.\n")
    cat(" (+) The true ordering may be umbrella! In fact, ")
    if(p2 > 0) cat(paste(classes[c(2,1,3)], c(">", "<", ""), sep = " "),"\n")
    if(p2 < 0) cat(paste(classes[c(2,1,3)], c("<", ">", ""), sep = " "),"\n")
  }
  d <- numeric(length(D.temp))
  d[D.temp == 1] <- which(id.median == 1)
  d[D.temp == 2] <- which(id.median == 2)
  d[D.temp == 3] <- which(id.median == 3)
  d[is.na(D.temp)] <- NA
  d1 <- d2 <- d3 <- rep(0, length(d))
  d1[d == 1] <- 1
  d2[d == 2] <- 1
  d3[d == 3] <- 1
  if (any(is.na(d))) {
    d1[is.na(d)] <- d2[is.na(d)] <- d3[is.na(d)] <- NA
  }
  Dvec <- cbind(d1, d2, d3)
  colnames(Dvec) <- c("D1", "D2", "D3")
  if(plot){
    boxplot(T ~ d, names = classes[id.median], na.action = na.omit,
            col = c("gray", "blue", "red"), ylab = "Diagnostic Test",
            xlab = "Three Ordered Groups")
  }
  return(list(D = d, Dvec = Dvec, order = id.median))
}
