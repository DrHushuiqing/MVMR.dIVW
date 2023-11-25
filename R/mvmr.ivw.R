#' Perform inverse-variance weighted (IVW) estimator for two-sample summary-data multivariable Mendelian randomization
#'
#' @param beta.exposure A data.frame or matrix. Each row contains the estimated marginal effect of a SNP on K exposures, usually obtained from a GWAS
#' @param se.exposure A data.frame or matrix of estimated standard errors of beta.exposure
#' @param beta.outcome A vector of the estimated marginal effect of a SNP on outcome, usually obtained from a GWAS
#' @param se.outcome A vector of estimated standard errors of beta.outcome
#' @param P A K-by-K matrix for the estimated shared correlation matrix, where K is the number of exposure
#'
#' @return A list with elements
#' \item{beta.hat}{Estimated direct effects of each exposure on the outcome}
#' \item{beta.se}{Estimated standard errors of beta.hat}
#' \item{iv_strength_parameter}{The minimum eigenvalue of the sample IV strength matrix, which quantifies the IV strength in the sample}
#' @import MVMR
#' @export
#'
#' @examples
#' library(MVMR)
#' data("rawdat_mvmr")
#' beta.exposure <- rawdat_mvmr[,c("LDL_beta","HDL_beta","Trg_beta")]
#' se.exposure <- rawdat_mvmr[,c("LDL_se","HDL_se","Trg_se")]
#' beta.outcome <- rawdat_mvmr$SBP_beta
#' se.outcome <- rawdat_mvmr$SBP_se
#' P <- matrix(0.3, nrow = 3, ncol = 3)
#' diag(P) <- 1
#' mvmr.ivw(beta.exposure = beta.exposure,
#' se.exposure = se.exposure,
#' beta.outcome = beta.outcome,
#' se.outcome = se.outcome,
#' P = P)
#'
mvmr.ivw <- function(beta.exposure, se.exposure, beta.outcome, se.outcome, P) {
  if (ncol(beta.exposure) <= 1 | ncol(se.exposure) <= 1) {stop("this function is developed for multivariable MR")}
  if (ncol(P) != ncol(beta.exposure)) {stop("The shared correlation matrix has a different number of columns than the input beta.exposure")}
  if (nrow(beta.exposure) != length(beta.outcome)) {stop("The number of SNPs in beta.exposure and beta.outcome is different")}
  beta.exposure <- as.matrix(beta.exposure)
  se.exposure <- as.matrix(se.exposure)
  # the number of SNPs
  p <- nrow(beta.exposure)
  # the number of exposures
  K <- ncol(beta.exposure)
  # diagonal W matrix
  W<- diag(se.outcome^(-2))
  # create a list of Sigma Xj matrices
  Vj <- lapply(1:p, function(j) diag(se.exposure[j,]) %*% P %*% diag(se.exposure[j,]))
  # calculate square root inverse of Vj
  Vj_root_inv <- lapply(Vj, function(x) {
    x_eigen <- eigen(x)
    x_eigen$vectors %*% diag(1/sqrt(x_eigen$values)) %*% t(x_eigen$vectors)
  })
  # calcualte IV strenght parameter
  IV_strength_matrix <- Reduce("+",lapply(1:p, function(j) {
    beta.exposure.V <- Vj_root_inv[[j]] %*% beta.exposure[j,]
    beta.exposure.V %*% t(beta.exposure.V)})) - p*diag(K)
  iv_strength_parameter <- min(eigen(IV_strength_matrix/sqrt(p))$values)
  # get V matrix
  V <- Reduce("+",lapply(1:p, function(j) {Vj[[j]] * (se.outcome[j]^(-2))}))
  # get M matrix
  M <- t(beta.exposure)%*%W%*%beta.exposure
  # multivariable IVW estimator
  mvmr.IVW <- solve(M) %*% t(beta.exposure) %*% W %*% beta.outcome
  # calculate variance covariance matrix
  IVW_Vt <- Reduce("+", lapply(1:p, function(j) {
    m <- beta.exposure[j,] %*% t(beta.exposure[j,]) * (se.outcome[j]^(-2))
    v <- Vj[[j]] * (se.outcome[j]^(-2))
    bvb <- as.numeric(t(mvmr.IVW) %*% v %*% mvmr.IVW)
    vbbv <- v %*% mvmr.IVW %*% t(mvmr.IVW) %*% v
    m*(1+bvb) + vbbv
  }))
  mvmr.IVW.se <- sqrt(diag(solve(M)%*%IVW_Vt%*%solve(M)))
  return(list(beta.hat = mvmr.IVW,
              beta.se = mvmr.IVW.se,
              iv_strength_parameter = iv_strength_parameter))
}
