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