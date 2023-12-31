mvmr.divw <- function(beta.exposure, se.exposure, beta.outcome, se.outcome, P, phi_cand=0, overlap = FALSE) {
  if (ncol(beta.exposure) <= 1 | ncol(se.exposure) <= 1) {stop("this function is developed for multivariable MR")}
  if ((!overlap) & ncol(P) != ncol(beta.exposure)) {stop("With independent exposure and outcome datasets, please provide an K-by-K estimated shared correlation matrix, where K is the number of exposures")}
  if (nrow(beta.exposure) != length(beta.outcome)) {stop("The number of SNPs in beta.exposure and beta.outcome is different")}
  if (overlap & (ncol(P) == ncol(beta.exposure))) {stop("With overlapping exposure and outcome datasets, please provide an (K+1)-by-(K+1) estimated shared correlation matrix, where K is the number of exposures")}
  beta.exposure <- as.matrix(beta.exposure)
  se.exposure <- as.matrix(se.exposure)
  # the number of SNPs
  p <- nrow(beta.exposure)
  # the number of exposures
  K <- ncol(beta.exposure)
  # diagonal W matrix
  W<- diag(se.outcome^(-2))
  # create a list of Sigma Xj matrices
  if (!overlap) {
    Vj <- lapply(1:p, function(j) diag(se.exposure[j,]) %*% P %*% diag(se.exposure[j,]))
  } else {
    Vj <- lapply(1:p, function(j) diag(c(se.exposure[j,],se.outcome[j])) %*% P %*% diag(c(se.exposure[j,],se.outcome[j])))
  }
  # calculate square root inverse of Vj
  Vj_root_inv <- lapply(Vj, function(x) {
    x_eigen <- eigen(x[1:K,1:K])
    x_eigen$vectors %*% diag(1/sqrt(x_eigen$values)) %*% t(x_eigen$vectors)
  })
  # calcualte IV strenght parameter
  IV_strength_matrix <- Reduce("+",lapply(1:p, function(j) {
    beta.exposure.V <- Vj_root_inv[[j]] %*% beta.exposure[j,]
    beta.exposure.V %*% t(beta.exposure.V)})) - p*diag(K)
  iv_strength_parameter <- min(eigen(IV_strength_matrix/sqrt(p))$values)
  # get V matrix
  V <- Reduce("+",lapply(1:p, function(j) {Vj[[j]][1:K,1:K] * (se.outcome[j]^(-2))}))
  # get M matrix
  M <- t(beta.exposure)%*%W%*%beta.exposure
  # get M-V matrix
  MV <- M-V
  # perform eigen-decomposition on MV
  MV_eigvalues <- eigen(MV)$values
  MV_eigen <- eigen(MV)
  # multivariable (a)dIVW estimator
  if (is.null(phi_cand)) {
  phi_cand <- c(0, exp(seq(0, 15, by = 0.5) - min(iv_strength_parameter)/2))
  }
  phi_length <- length(phi_cand)
  MV.l.inv.long <- Reduce(rbind, lapply(1:phi_length, function(l) {
    MV_eigen$vectors %*% diag(1/(MV_eigvalues + phi_cand[l]/MV_eigvalues)) %*% t(MV_eigen$vectors)}
  ))
  if (!overlap) {
    # (a)dIVW for independent datasets
    beta.est <- MV.l.inv.long %*% t(beta.exposure) %*% W %*% (beta.outcome)
    prof.lik <- sapply(1:phi_length, function(l) {
      beta.hat <- beta.est[(1+(l-1)*K):(l*K)]
      bvb.test <- sapply(1:p, function(j) t(beta.hat) %*% Vj[[j]][1:K,1:K] %*% beta.hat)
      S <- diag(1/(se.outcome^2 + bvb.test))
      1/p * t(beta.outcome - beta.exposure %*% beta.hat) %*% S %*%
        (beta.outcome - beta.exposure %*% beta.hat)})
    phi_selected <- phi_cand[which.min(prof.lik)]
    MV.l.inv <- MV_eigen$vectors %*% diag(1/(MV_eigvalues + phi_selected/MV_eigvalues)) %*% t(MV_eigen$vectors)
    mvmr.adIVW <- MV.l.inv %*% t(beta.exposure) %*% W %*% beta.outcome
    adIVW_Vt <- Reduce("+",lapply(1:p, function(j) {
      m <- beta.exposure[j,] %*% t(beta.exposure[j,]) * (se.outcome[j]^(-2))
      v <- Vj[[j]][1:K,1:K]*(se.outcome[j]^(-2))
      bvb <- as.numeric(t(mvmr.adIVW) %*% v %*% mvmr.adIVW)
      vbbv <- v %*% mvmr.adIVW %*% t(mvmr.adIVW) %*% v
      m*(1+bvb) + vbbv
    }))
    mvmr.adIVW.se <- sqrt(diag(MV.l.inv%*%adIVW_Vt%*%MV.l.inv))
  } else {
    # (a)dIVW allowing for overlap
    Adj_term <- Reduce("+",lapply(1:p, function(j) {Vj[[j]][1:K,(K+1)] * se.outcome[j]^{-2}}))
    beta.est <- MV.l.inv.long %*% t(beta.exposure) %*% W %*% (beta.outcome) - MV.l.inv.long %*% Adj_term
    prof.lik <- sapply(1:phi_length, function(l) {
      beta.hat <- beta.est[(1+(l-1)*K):(l*K)]
      bvb.test <- sapply(1:p, function(j) t(beta.hat) %*% Vj[[j]][1:K,1:K] %*% beta.hat)
      bvxy <- sapply(1:p, function(j) {
        sigmaxy <- Vj[[j]][1:K,K+1]
        (t(beta.hat) %*% sigmaxy)[1,1]
      })
      S <- diag(1/(se.outcome^2 + bvb.test - 2 * bvxy))
      1/p * t(beta.outcome - beta.exposure %*% beta.hat) %*% S %*%
        (beta.outcome - beta.exposure %*% beta.hat)})
    phi_selected <- phi_cand[which.min(prof.lik)]
    MV.l.inv <- MV_eigen$vectors %*% diag(1/(MV_eigvalues + phi_selected/MV_eigvalues)) %*% t(MV_eigen$vectors)
    mvmr.adIVW <- MV.l.inv %*% t(beta.exposure) %*% W %*% beta.outcome - MV.l.inv %*% Adj_term
    adIVW_Vt_overlap <- Reduce("+",lapply(1:p, function(j) {
      m <- beta.exposure[j,] %*% t(beta.exposure[j,]) * (se.outcome[j]^(-2))
      v <- Vj[[j]][1:K,1:K] * (se.outcome[j]^(-2))
      bvb <- as.numeric(t(mvmr.adIVW) %*% v %*% mvmr.adIVW)
      vbbv <- v %*% mvmr.adIVW %*% t(mvmr.adIVW) %*% v
      sigmaxy <- Vj[[j]][1:K,K+1]
      A1 <- sigmaxy %*% t(sigmaxy) * se.outcome[j]^(-4)
      A2 <- Vj[[j]][1:K,1:K] %*% mvmr.adIVW %*% t(sigmaxy) * se.outcome[j]^(-4)
      A3 <- sigmaxy %*% t(mvmr.adIVW) %*% Vj[[j]][1:K,1:K] * se.outcome[j]^(-4)
      A4 <- (t(mvmr.adIVW) %*% sigmaxy * se.outcome[j]^(-2))[1,1] *  m
      A5 <- (t(mvmr.adIVW) %*% sigmaxy * se.outcome[j]^(-2))[1,1] *  sigmaxy %*% t(sigmaxy) * se.outcome[j]^(-4)
      m*(1+bvb) + vbbv + A1 - A2 - A3 - 2 * A4 + 4 * A5
    }))
    mvmr.adIVW.se <- sqrt(diag(MV.l.inv%*%adIVW_Vt_overlap%*%MV.l.inv))
  }
  return(list(beta.hat = mvmr.adIVW,
              beta.se = mvmr.adIVW.se,
              iv_strength_parameter = iv_strength_parameter,
              phi_selected = phi_selected))
}
