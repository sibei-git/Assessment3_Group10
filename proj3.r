# task1
par(mfrow=c(1,2))
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
  knots <- c(
    min(knots_all) - 3*dk, min(knots_all) - 2*dk, min(knots_all) - 1*dk,
    knots_all,
    max(knots_all) + 1*dk, max(knots_all) + 2*dk, max(knots_all) + 3*dk
  )
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
  list(Xtilde = Xtilde, X = X, S = S, tf = tf, knots = knots, pi = pd)
}

##data
data <- read.table('/Users/kangxinwei/sds/statistical programming/Assessment3_Group10/engcov .txt', header = TRUE)
t_obs <- data$julian
y <- data$nhs
out<-Xbar_X_S(t_obs)

#task 2
negloglik_penalty <- function(gamma, X, y, S, lambda) {
  beta <- exp(gamma)              
  mu <- as.vector(X %*% beta)       
  mu <- pmax(mu, 1e-10)         
  
  loglik <- sum(y * log(mu) - mu)  # no log(yi!),after derivation, it eauals to 0
  penalty <- 0.5 * lambda * t(beta) %*% S %*% beta
  
  return(-loglik + penalty)     
}

grad_negloglik_penalty <- function(gamma, X, y, S, lambda){
  beta <- exp(gamma)
  mu <- as.vector(X %*% beta)
  mu <- pmax(mu, 1e-10)
  grad_beta <- t(X) %*% (1 - y / mu)
  grad_gamma <- as.vector(beta * grad_beta + lambda * (beta * (S %*% beta)))
  return(grad_gamma)
}

# xga <- rgamma(K)
# ga0 <- negloglik_penalty(xga, X, y, S, lambda = 5e-5)
# ga1 <- rep(0,K)
# eps <- 1e-8
# for (i in 1:K){
#   xa_0 <- xga
#   xa_0[i] <- xga[i]+eps
#   ga1[i] <- negloglik_penalty(xa_0, X, y, S, lambda = 5e-5)
# }
# fd <- (ga1-ga0)/eps
# fd1 <- grad_negloglik_penalty(xga, X, y, S, lambda = 5e-5)
# fd1-fd

#task3
K <- ncol(X)
init_gamma <- rep(0, K)
opt <- optim(
  par = init_gamma,
  fn = negloglik_penalty,
  gr = grad_negloglik_penalty,
  method = "BFGS",
  control = list(maxit = 200, trace = TRUE),
  X = X, y = y, S = S, lambda = 5e-5
)

opt$convergence  # 0 表示收敛
gamma_hat <- opt$par
beta_hat <- exp(gamma_hat)
mu_hat <- as.vector(X %*% beta_hat)
f_hat <- as.vector(Xtilde %*% beta_hat)

# 计算所有数据的范围，确保y轴能容纳三条曲线
y_all <- range(c(y, mu_hat, f_hat))

# 绘图
plot(t_obs, y, type = "p", col = "grey",
     main = "Observed, Fitted Deaths and Estimated Infections",
     xlab = "Time (days)", ylab = "Count",
     ylim = y_all)

# 添加拟合死亡曲线
lines(t_obs, mu_hat, col = "blue", lwd = 2)

# 添加感染曲线
lines(out$tf, f_hat, col = "red", lwd = 2)

# 添加图例
legend("topleft",
       legend = c("Observed Deaths", "Fitted Deaths", "Estimated Infections"),
       col = c("grey", "blue", "red"),
       lty = c(NA, 1, 1), pch = c(1, NA, NA),
       bty = "n")  # 无边框图例

# 若想标出观测期起点：
abline(v = t_obs[1], lty = 2, col = "darkgrey")

