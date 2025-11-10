

# task1

library(splines)

#probability of each fatal disease duration(pi(j))
d <- 1:80;edur <- 3.151;sdur <- .469
pd <- dlnorm(d,edur,sdur);
pd <- pd/sum(pd)

#create Xtilde_X_S
Xtilde_X_S <- function(t,K=80,L=80,L0=30){
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


#task3
fit_death <- function(y, X, Xtilde, S, lambda = 5e-5, maxit = 200) {
  # 初始化
  K <- ncol(X)
  init_gamma <- rep(0, K)
  
  # 使用 BFGS 优化
  opt <- optim(
    par = init_gamma,
    fn = negloglik_penalty,
    gr = grad_negloglik_penalty,
    method = "BFGS",
    control = list(maxit = maxit, trace = TRUE),
    X = X, y = y, S = S, lambda = lambda
  )
  
  if (opt$convergence != 0) {
    warning("Optimization did not converge. Try adjusting lambda or maxit.")
  }
  
  # 计算估计值
  gamma_hat <- opt$par
  beta_hat  <- exp(gamma_hat)
  mu_hat    <- as.vector(X %*% beta_hat)
  f_hat     <- as.vector(Xtilde %*% beta_hat)
  
  return(list(
    opt = opt,
    gamma_hat = gamma_hat,
    beta_hat = beta_hat,
    mu_hat = mu_hat,
    f_hat = f_hat,
    lambda = lambda
  ))
}


plot_death <- function(t_obs, y, tf, out, mu_hat, f_hat, lambda = NULL) {
  y_all <- range(c(y, mu_hat, f_hat))
  
  # 绘制实际死亡点
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
  
}

#task4
log_lambda_interval <- seq(-13,-7,length=50)
bic_criterion <- function(X,S,y,log_lambda_interval){
  n <- length(y)
  K <- ncol(X)
  init_gamma <- rep(0, K)
  bic_interval <- numeric(0)
  lambda_interval <- exp((log_lambda_interval))
  for (lam in lambda_interval){
    opt <- optim(
      par = init_gamma,
      fn = negloglik_penalty,
      gr = grad_negloglik_penalty,
      method = "BFGS",
      control = list(maxit = 200, trace = TRUE),
      X = X, y = y, S = S, lambda = lam
    )
    init_gamma <-opt$par
    beta <- exp(opt$par)
    mu <- pmax(as.vector(X %*% beta), 1e-10)
    loglik <- sum(y * log(mu) - mu)
    w <- as.numeric((y/(mu^2)))
    H_lambda <- t(X) %*% (X * w) + lam*S
    H_lambda <- H_lambda+diag(1e-8, ncol(X))
    H_0 <-  t(X) %*% (X * w)
    edf <- sum(diag(solve(H_lambda,H_0)))
    bic <- -2*loglik + log(n)*edf
    bic_interval <- c(bic_interval,bic)
  }
  lambda_hat <- lambda_interval[which.min(bic_interval)]
  return(lambda_hat)
}

#task5
w_negloglik_penalty <- function(gamma, wb, X, y, S, lambda) {
  beta <- exp(gamma)              
  mu <- as.vector(X %*% beta)       
  mu <- pmax(mu, 1e-10)         
  
  loglik <- sum(wb*(y * log(mu) - mu))  # no log(yi!),after derivation, it eauals to 0
  penalty <- 0.5 * lambda * t(beta) %*% S %*% beta
  
  return(-loglik + penalty)     
}

w_grad_negloglik_penalty <- function(gamma, wb, X, y, S, lambda){
  beta <- exp(gamma)
  mu <- as.vector(X %*% beta)
  mu <- pmax(mu, 1e-10)
  grad_beta <- t(X) %*% (wb*(1 - y / mu))
  grad_gamma <- as.vector(beta * grad_beta + lambda * (beta * (S %*% beta)))
  return(grad_gamma)
}

bootstrap <- function(X, Xtilde, S, y, log_lambda_interval){
  lambda_hat <- bic_criterion(X,S,y,log_lambda_interval)
  
  n <- length(y)
  K <- ncol(X)
  
  wb_mat <- replicate(200, tabulate(sample(n, replace = TRUE), n))
  
  fhat_boot <- matrix(NA, nrow = nrow(Xtilde), ncol = 200)
  
  for(i in 1:200){
    wb <- wb_mat[, i]
    opt <- opt <- optim(
      par = rep(0, K),
      fn = w_negloglik_penalty,
      gr = w_grad_negloglik_penalty,
      method = "BFGS",
      control = list(maxit = 200, trace = FALSE),
      X = X, y = y, S = S, lambda = lambda_hat, wb = wb
    )
    beta_i <- exp(opt$par)
    fhat_boot[, i] <- as.vector(Xtilde %*% beta_i)
  }
  return(list(
    fhat_boot = fhat_boot,
    lambda_hat = lambda_hat,
    wb_mat = wb_mat
  ))
}


#task6
plot_bootstrap <- function(tf, fhat_boot) {
  # 计算均值与置信区间
  f_mean <- rowMeans(fhat_boot)
  f_low  <- apply(fhat_boot, 1, quantile, probs = 0.025)
  f_high <- apply(fhat_boot, 1, quantile, probs = 0.975)
  
  # 绘图
  plot(tf, f_mean, type = "l", col = "red", lwd = 2,
       main = "Estimated infection curve with 95% bootstrap CI",
       xlab = "Time (days)", ylab = "f(t)")
  
  lines(tf, f_low, col = "grey60", lty = 2)
  lines(tf, f_high, col = "grey60", lty = 2)
  
  # 灰色置信区间阴影
  polygon(c(tf, rev(tf)), c(f_low, rev(f_high)),
          col = rgb(0.9, 0.9, 0.9, 0.5), border = NA)
}

data <- read.table('D:/studying/Extended statistical programming/Assessment/engcov.txt', header = TRUE)

# Generate X, S
t_obs <- data$julian
y <- data$nhs
out<-Xtilde_X_S(t_obs)
X <- out$X
Xtilde <- out$Xtilde
S <- out$S
tf <- out$tf
par(mfrow = c(2, 2))
# 用固定 λ = 5e-5 拟合并绘图
fit_result <- fit_death(y, X, Xtilde, S, lambda = 5e-5)
plot_death(t_obs, y, tf, out, fit_result$mu_hat, fit_result$f_hat, lambda = fit_result$lambda)

# BIC 搜索最优 λ 
log_lambda_interval <- seq(-13, -7, length = 50)
lambda_hat <- bic_criterion(X, S, y, log_lambda_interval)
cat("Best lambda found by BIC:", lambda_hat, "\n")

# 用最优 λ 重新拟合并绘图
fit_best <- fit_death(y, X, Xtilde, S, lambda = lambda_hat)
plot_death(t_obs, y, tf, out, fit_best$mu_hat, fit_best$f_hat, lambda = lambda_hat)

# Bootstrap 重采样计算置信区间
boot_result <- bootstrap(X, Xtilde, S, y, log_lambda_interval)
fhat_boot <- boot_result$fhat_boot
lambda_hat <- boot_result$lambda_hat
cat("Bootstrap finished, using lambda =", lambda_hat, "\n")

# 绘制感染曲线及 95% 置信带
plot_bootstrap(tf, fhat_boot)
