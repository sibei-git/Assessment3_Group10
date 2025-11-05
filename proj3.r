# task1

library(splines)

#probability of each fatal disease duration(pi(j))
d <- 1:80;edur <- 3.151;sdur <- .469
pd <- dlnorm(d,edur,sdur);
pd <- pd/sum(pd)

#create Xbar_X_S
Xbar_X_S <- function(t,K=80,L=80,L0=30){
  t <- as.numeric(t)
  n <- length(t)
  
  # t-axis of f
  tf <- seq(from = t[1]-L0, to = t[n], by = 1)
  nf <- length(tf)
  
  # B-spline basis(number=K) 
  knots_all <- seq(min(tf), max(tf), length.out = K - 2)
  dk <- mean(diff(knots_all))
  knots <- c(rep(min(knots_all) - 2*dk, 2), knots_all,
             rep(max(knots_all) + 2*dk, 2))
  Xtilde <- splineDesign(knots = knots, x = tf, ord = 4, outer.ok = TRUE)
  
  # model matrix of deconvolution
  X <- matrix(0.0, n, K)
  row_idx <- L0 + outer(1:n, 1:L, FUN = "-")
  valid_v <- (row_idx >= 1) & (row_idx <= nf)
  for (j in 1:L) {
    idx <- which(valid_v[, j])
    if (length(idx)) {
      for (ii in seq_along(idx)) {
        X[idx[ii], ] <- X[idx[ii], ] + pd[j] * Xtilde[row_idx[idx[ii], j], ]
      }
    }
  }
  # Penalty matrix S
  S <- crossprod(diff(diag(K), diff = 2))
  list(Xbar = Xbar, X = X, S = S, tf = tf, knots = knots, pi = pd)
}
# data <- read.table('/Users/youchao/Desktop/Extended Statistical Programming for lecture/assignment 3/engcov.txt', header = TRUE)
# t_obs <- data$julian
# y <- data$nhs
# out<-Xbar_X_S(t_obs)

