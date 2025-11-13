# Chao You s2785107; Chenyi Jia s2792046; Kangxin Wei s2817655
# Our team members completed this assignment together, completing each step 
# through discussion.In the early stage, we read the notes together to 
# understand the requirements and learn functions. Later, we program together, 
#find bugs, and solve them together.

# https://github.com/sibei-git/Assessment3_Group10.git
#-------------------------------------------------------------------------------
# Smooth Deconvolution
# This script estimates the underlying infection curve f(t)
# from observed daily deaths y(t) using:
#  a known infection-to-death delay distribution π(j),
#  B-spline basis for smooth representation of f(t),
#  penalized Poisson regression (with smoothness penalty λβᵀSβ),
#  automatic selection of λ using BIC,
#  bootstrap estimation for uncertainty intervals.
#
# Output:
#   1. Fitted deaths vs observed deaths plot.
#   2. Estimated infection curve f(t).
#   3. Infection curve with 95% bootstrap confidence interval.
#-------------------------------------------------------------------------------
library(splines)  # provides splineDesign() for B-spline basis
# Define infection-to-death probability distribution π(j)
# log(d) ~ Normal(3.151, 0.469^2)
d <- 1:80                # possible durations from infection to death
edur <- 3.151            # Mean of log-normal delay (log-scale)
sdur <- .469             # SD of log-normal delay
pd <- dlnorm(d, edur, sdur)   # Log-normal distribution of delay
pd <- pd / sum(pd)       # Normalize so that π(j) sums to 1


# Task 1: Function to build matrices for deconvolution
# Function: Constructs the B-spline basis (Xtilde),
# convolution model matrix (X), and penalty matrix (S).
# Arguments:
#   t: Numeric vector of observed time points.
#   K: Number of spline basis functions (default = 80).
#   L: Maximum infection-to-death delay (default = 80).
#   L0: Offset to extend f(t) grid earlier than first observed death (default = 30)
# Returns:
#   Xtilde: Spline basis matrix for infection curve f(t).
#   X: Convolution model matrix linking infections to deaths.
#   S: Penalty matrix for smoothness (2nd derivative).
#   tf: Time grid for estimated infections.
#   knots: Knots used for spline basis.
#   pi: Probability distribution π(j) for delay.
Xtilde_X_S <- function(t, K = 80, L = 80, L0 = 30) {
  t <- as.numeric(t)
  n <- length(t)  # Number of observed points
  
  # Construct time grid for f(t): start earlier to capture infection onset
  # tf:from 30 days before the first death (day 62) to the last death (day 211)
  tf <- seq(from = t[1] - L0, to = t[n], by = 1)
  nf <- length(tf)
  
  # Define internal knots and extend boundaries for cubic B-spline (ord=4)
  knots_all <- seq(min(tf), max(tf), length.out = K - 2)
  dk <- mean(diff(knots_all))  # Average knot spacing
  
  # Add 3 boundary knots at each side for smooth spline ends
  knots <- c(
    min(knots_all) - 3*dk, min(knots_all) - 2*dk, min(knots_all) - 1*dk,
    knots_all,
    max(knots_all) + 1*dk, max(knots_all) + 2*dk, max(knots_all) + 3*dk
  )
  
  # Build B-spline basis for infection curve f(t)
  # Xtilde: The [nf x K] B-spline basis matrix, 
  # where Xtilde[i, k] is the value of basis function k at time tf[i].
  Xtilde <- splineDesign(knots = knots, x = tf, ord = 4, outer.ok = TRUE)
  
  # Build the [n x K] Deconvolution Matrix X
  # This matrix maps the infection basis (Xtilde) to the observed deaths (mu)
  # by "convolving" the infection function f(t) with the delay probability pd.
  # In order to get mu = X %*% beta
  # Initialize convolution matrix X (deaths = convolution of infections)
  X <- matrix(0.0, n, K)
  # For each observed day, compute indices corresponding to past infection times
  row_idx <- L0 + outer(1:n, 1:L, FUN = "-")
  valid_v <- (row_idx >= 1) & (row_idx <= nf)
  
  # Loop over every possible infection-to-death delay 'j' (from 1 to L)
  # Combine contributions from each delay j weighted by π(j)
  for (j in 1:L) {
    # Find all observed death days 'i' for which a delay 'j' is valid.
    idx <- which(valid_v[, j])
    # If any valid days exist for this delay 'j'
    if (length(idx)) {
      # Loop over each valid death day 'i' (indexed by 'ii')
      for (ii in seq_along(idx)) {
        # X[i, ] = X[i, ] + (probability of j-day delay) * (B-spline basis)
        # Xtilde[row_idx[...], ] gets the K basis values at time (t_i - j).
        X[idx[ii], ] <- X[idx[ii], ] + pd[j] * Xtilde[row_idx[idx[ii], j], ]
      }
    }
  }
  
  # Construct 2nd difference penalty matrix (encourages smoothness)
  S <- crossprod(diff(diag(K), diff = 2))
  
  # Return all components
  list(Xtilde = Xtilde, X = X, S = S, tf = tf, knots = knots, pi = pd)
}

# Task 2: Penalized negative log-likelihood and gradient
# Function: Computes penalized negative log-likelihood for Poisson deaths.
# Arguments:
#   gamma: log(β) coefficients to ensure β > 0.
#   X, y: Model matrix and observed deaths.
#   S: Smoothness penalty matrix.
#   lambda: Smoothing parameter controlling smoothness of f(t).
# Returns:
#   Scalar value of penalized negative log-likelihood.
negloglik_penalty <- function(gamma, X, y, S, lambda) {
  beta <- exp(gamma) # β = exp(γ) ensures positivity
  mu <- as.vector(X %*% beta) # Expected deaths
  mu <- pmax(mu, 1e-10)  # Avoid log(0)
  
  loglik <- sum(y * log(mu) - mu) # Poisson log-likelihood (without constant)
  penalty <- 0.5 * lambda * t(beta) %*% S %*% beta  # Smoothness penalty
  
  return(-loglik + penalty) # Return penalized negative log-likelihood
}

# Function: Gradient of penalized negative log-likelihood w.r.t γ.
# Arguments and Returns same as above.
grad_negloglik_penalty <- function(gamma, X, y, S, lambda) {
  beta <- exp(gamma)
  mu <- as.vector(X %*% beta)
  mu <- pmax(mu, 1e-10)
  
  grad_beta <- t(X) %*% (1 - y / mu) # Gradient w.r.t β
  grad_gamma <- as.vector(beta * grad_beta + lambda * (beta * (S %*% beta)))# Chain rule
  return(grad_gamma)
}

# Task 3: Fit penalized Poisson model for deaths
# Function: Fits penalized Poisson regression using BFGS optimization.
# Arguments:
#   y, X, Xtilde, S: Model inputs.
#   lambda: Smoothing parameter.
#   maxit: Maximum iterations for optimizer.
# Returns:
#   List containing estimated parameters and fitted values:
#   gamma_hat, beta_hat, mu_hat (fitted deaths), f_hat (estimated infections)
fit_death <- function(y, X, Xtilde, S, lambda = 5e-5, maxit = 200) {
  K <- ncol(X)
  init_gamma <- rep(0, K)  # Initialize γ = 0
  
  # Optimize penalized likelihood via BFGS
  opt <- optim(
    par = init_gamma,# par:The initial guess for the parameters (gamma vector)
    fn = negloglik_penalty,# fn:The objective function to minimize
    gr = grad_negloglik_penalty,# gr:The gradient function for 'fn'
    method = "BFGS", # use the BFGS algorithm
    control = list(maxit = maxit, trace = TRUE),
    X = X, y = y, S = S, lambda = lambda
  )
  
  # Warn if optimizer did not converge
  if (opt$convergence != 0) {
    warning("Optimization did not converge. Try adjusting lambda or maxit.")
  }
  
  # Recover fitted parameters and values
  gamma_hat <- opt$par #it stores the optimal parameter vector found by optim
  beta_hat <- exp(gamma_hat)
  # Calculate the fitted death curve (mu_hat)
  # Use the optimal beta_hat and the deconvolution matrix X (mu = X * beta)
  mu_hat <- as.vector(X %*% beta_hat)
  # Calculate the estimated infection curve (f_hat)
  # Use the optimal beta_hat and the Xtilde (f = Xtilde * beta)
  f_hat <- as.vector(Xtilde %*% beta_hat)
  
  return(list(
    opt = opt, gamma_hat = gamma_hat, beta_hat = beta_hat,
    mu_hat = mu_hat, f_hat = f_hat, lambda = lambda
  ))
}


# Function: Plot observed deaths, fitted deaths, and infection curve.
# Arguments:
#   t_obs: Observed time vector.
#   y: Observed deaths.
#   tf: Time grid for infection curve.
#   out: Output from Xtilde_X_S() containing tf and knots.
#   mu_hat: Fitted deaths from model.
#   f_hat: Estimated infections f(t).
#   lambda: (Optional) smoothing parameter for plot title.
# Returns:
#  plots figure directly
plot_death <- function(t_obs, y, tf, out, mu_hat, f_hat, lambda = NULL) {
  y_all <- range(c(y, mu_hat, f_hat))  # Determine y-axis range
  
  # Plot observed deaths (grey points)
  plot(t_obs, y, type = "p", col = "darkgreen",
       main = "Observed, Fitted Deaths and Estimated Infections",
       xlab = "Time (days)", ylab = "Count",
       ylim = y_all)
  
  # Add fitted death curve (blue)
  lines(t_obs, mu_hat, col = "blue", lwd = 2)
  
  # Add estimated infection curve (red)
  lines(out$tf, f_hat, col = "red", lwd = 2)
  
  # Add legend and start marker
  legend("topright",
         legend = c("Observed Deaths", "Fitted Deaths", "Estimated Infections"),
         col = c("darkgreen", "blue", "red"),
         lty = c(NA, 1, 1), pch = c(1, NA, NA),
         bty = "n")
  abline(v = t_obs[1], lty = 2, col = "darkgrey")
}

# Task 4: Select λ using Bayesian Information Criterion (BIC)
# Function: Performs a grid search over a range of λ values to find the
# one that minimizes the BIC.
# Arguments:
#   X, S, y: Model inputs (death matrix, penalty matrix, observed deaths).
#   log_lambda_interval: A sequence of log(λ) values to search over.
# Returns:
#   lambda_hat: The single λ value from the interval that results
#   in the minimum BIC score.
bic_criterion <- function(X, S, y, log_lambda_interval) {
  n <- length(y) # n: Sample size (number of observed days)
  K <- ncol(X) # K: Number of B-spline basis functions
  init_gamma <- rep(0, K) # Initial guess for gamma
  bic_interval <- numeric(0) # Vector to store BIC scores for each λ
  lambda_interval <- exp(log_lambda_interval)#Convert from log-scale to λ values
  
  # Loop over every candidate λ value in the grid
  for (lam in lambda_interval) {
    # Fit model for current λ
    opt <- optim(
      par = init_gamma,
      fn = negloglik_penalty,
      gr = grad_negloglik_penalty,
      method = "BFGS",
      control = list(maxit = 200, trace = TRUE),
      X = X, y = y, S = S, lambda = lam
    )
    # Use the result of this fit as a "warm start" for the next iteration
    # This makes the search much faster
    init_gamma <- opt$par 
    
    # Compute log-likelihood
    beta <- exp(opt$par)
    mu <- pmax(as.vector(X %*% beta), 1e-10)
    loglik <- sum(y * log(mu) - mu)
    
    # Compute effective degrees of freedom (edf)
    # edf = trace(H_lambda_inv * H_0)
    
    # Weight matrix w = diag(y_i / mu_i^2)
    w <- as.numeric(y / (mu^2))
    # Penalized Hessian: H_lambda = X'WX + λS
    H_lambda <- t(X) %*% (X * w) + lam * S
    # Add a small stabilization term (ridge) to ensure H_lambda is invertible
    H_lambda <- H_lambda + diag(1e-8, ncol(X))
    # Unpenalized Hessian: H_0 = X'WX
    H_0 <- t(X) %*% (X * w)
    edf <- sum(diag(solve(H_lambda, H_0))) # Calculate edf
    
    # Compute BIC
    # BIC = -2l + log(n) * EDF
    bic <- -2 * loglik + log(n) * edf
    bic_interval <- c(bic_interval, bic)# store the BIC score for this λ
  }
  # Find the λ that corresponds to the minimum BIC score
  lambda_hat <- lambda_interval[which.min(bic_interval)]#Choose λ with smallest BIC
  return(lambda_hat)
}

# Task 5: Bootstrap resampling to estimate uncertainty in f(t)
# Function:Computes the weighted penalized negative log-likelihood
# for the Poisson model used in bootstrap resampling.
# Arguments:
#   gamma: log(β) coefficients to ensure β > 0.
#   wb: Vector of bootstrap weights indicating how many times 
#            each observation was sampled.
#   X, S, y: Model inputs (death matrix, penalty matrix, observed deaths).
#   lambda: Smoothing parameter controlling penalty strength.
# Returns:
#   A single numeric value representing the penalized negative log-likelihood.
w_negloglik_penalty <- function(gamma, wb, X, y, S, lambda) {
  beta <- exp(gamma) #β = exp(γ) ensures positivity
  mu <- as.vector(X %*% beta) # Compute fitted mean values
  mu <- pmax(mu, 1e-10) # Avoid log(0) or division by zero by bounding μ below 1e-10
  
  # Compute the weighted Poisson log-likelihood (ignoring constant term)
  # Each observation is weighted by its bootstrap weight wb
  loglik <- sum(wb * (y * log(mu) - mu))
  
  # Compute smoothness penalty term
  penalty <- 0.5 * lambda * t(beta) %*% S %*% beta
  
  #Return the penalized negative log-likelihood
  return(-loglik + penalty)
}

# Function:Computes the gradient of the weighted penalized negative
# log-likelihood with respect to γ = log(β).
# Arguments:
#   gamma: log(β) coefficients to ensure β > 0.
#   wb: Numeric vector of bootstrap weights (length n).
#   X, y, S: Model inputs (death matrix, observed deaths, penalty matrix).
# Returns:
#   Numeric vector (length K) — gradient ∂L/∂γ.
w_grad_negloglik_penalty <- function(gamma, wb, X, y, S, lambda) {
  # Convert log-scale parameters gamma to positive coefficients beta
  beta <- exp(gamma)
  mu <- as.vector(X %*% beta) # Compute fitted mean values
  mu <- pmax(mu, 1e-10) # Avoid log(0) or division by zero
  
  # Compute the gradient of the weighted negative log-likelihood with respect to β
  grad_beta <- t(X) %*% (wb * (1 - y / mu))
  # Chain rule
  grad_gamma <- as.vector(beta * grad_beta + lambda * (beta * (S %*% beta)))
  return(grad_gamma)
}

# Function:Performs nonparametric bootstrap to estimate uncertainty in the 
# infection curve f(t). For each resample,fit weighted penalized model using 
# BFGS optimization and store estimated f(t) curve. Repeat this 200 times 
# to obtain bootstrap samples.
# Arguments:
#   X: Model and basis matrices .
#   Xtilde:Spline basis matrix for infection curve f(t).
#   S: Penalty matrix.
#   y: Observed deaths.
#   log_lambda_interval: A sequence of log(λ) values to search over.
# Returns:
#   A list containing:
#     fhat_boot: matrix of bootstrap infection curves.
#     lambda_hat: Optimal λ selected by BIC.
#     wb_mat: matrix of bootstrap weights.
bootstrap <- function(X, Xtilde, S, y, log_lambda_interval) {
  # Select λ using BIC
  lambda_hat <- bic_criterion(X, S, y, log_lambda_interval)
  
  n <- length(y)
  K <- ncol(X)
  
  # Generate 200 bootstrap weight vectors 
  wb_mat <- replicate(200, tabulate(sample(n, replace = TRUE), n))
  fhat_boot <- matrix(NA, nrow = nrow(Xtilde), ncol = 200)# Store infection estimates
  
  # Loop over 200 bootstrap samples: fit weighted model 
  # and store infection estimates
  for (i in 1:200) {
    wb <- wb_mat[, i]
    
    # Fit weighted model using BFGS
    opt <- optim(
      par = rep(0, K),
      fn = w_negloglik_penalty,
      gr = w_grad_negloglik_penalty,
      method = "BFGS",
      control = list(maxit = 200, trace = FALSE),
      X = X, y = y, S = S, lambda = lambda_hat, wb = wb
    )
    
    # Convert log-coefficients to positive beta values
    beta_i <- exp(opt$par)
    # Compute infection estimate for sample i
    fhat_boot[, i] <- as.vector(Xtilde %*% beta_i)
  }
  
  return(list(
    fhat_boot = fhat_boot,
    lambda_hat = lambda_hat,
    wb_mat = wb_mat
  ))
}


# Task 6: Plot infection curve with 95% bootstrap CI
# Function:Visualize the estimated infection curve f(t) along with its
# 95% bootstrap confidence interval (CI). The function computes
# the mean and quantile bounds across bootstrap samples and
# overlays them on a single plot.
# Arguments:
#   t_obs: Observed time points (days).
#   y: Observed deaths.
#   tf: Time grid for infection curve.
#   fit_best: Model object from fit_death() using optimal λ.
#   fhat_boot: Matrix of bootstrap infection curves .
# Returns:
#   None (plots directly).
plot_final <- function(t_obs, y, tf, fit_best, fhat_boot) {
  # Compute mean and 95% CI of infection curve f(t)
  f_low  <- apply(fhat_boot, 1, quantile, probs = 0.025)
  f_high <- apply(fhat_boot, 1, quantile, probs = 0.975)
  
  # Set up empty plot with sensible axis limits
  plot(tf, fit_best$f_hat, type = "n",
       main = "Daily deaths and infection rate with 95% CI",
       xlab = "Time (days)", ylab = "Count",
       ylim = c(0, max(y, f_high, na.rm = TRUE)))
  
  # Add shaded 95% CI for infection curve (dark grey)
  polygon(c(tf, rev(tf)), c(f_low, rev(f_high)),
          col = rgb(0.8, 0.6, 0.9, 0.4), border = NA)
  
  # Add observed and fitted deaths
  points(t_obs, y, col = "darkgreen", pch = 1)           # Observed deaths (grey points)
  lines(t_obs, fit_best$mu_hat, col = "blue", lwd = 2)  # Fitted deaths (blue line)
  
  # Add mean infection curve f(t)
  lines(tf, fit_best$f_hat, col = "red", lwd = 2)

  # Add legend
  legend("topright",
         legend = c("Observed Deaths", "Fitted Deaths",
                    "Estimated Infections", "95% CI"),
         col = c("darkgreen", "blue", "red",rgb(0.8, 0.6, 0.9, 0.6)),
         lty = c(NA, 1, 1, 1), pch = c(1, NA, NA, NA),
         bty = "n")
}

# Run the complete pipeline
data <- read.table('/Users/youchao/Assessment3_Group10/engcov.txt', header = TRUE)
t_obs <- data$julian
y <- data$nhs

# Build design matrices
out <- Xtilde_X_S(t_obs)
X <- out$X
Xtilde <- out$Xtilde
S <- out$S
tf <- out$tf

# Fit, plot, and bootstrap
# Fit with fixed λ = 5e-5 and plot observed deaths, fitted deaths, and infection curve
fit_result <- fit_death(y, X, Xtilde, S, lambda = 5e-5)
plot_death(t_obs, y, tf, out, fit_result$mu_hat, fit_result$f_hat, lambda = fit_result$lambda)

# Select λ via BIC
log_lambda_interval <- seq(-13, -7, length = 50)
lambda_hat <- bic_criterion(X, S, y, log_lambda_interval)
cat("Best lambda found by BIC:", lambda_hat, "\n")

# Refit using best λ
fit_best <- fit_death(y, X, Xtilde, S, lambda = lambda_hat)

# Bootstrap to estimate uncertainty
boot_result <- bootstrap(X, Xtilde, S, y, log_lambda_interval)
fhat_boot <- boot_result$fhat_boot
lambda_hat <- boot_result$lambda_hat
cat("Bootstrap finished, using lambda =", lambda_hat, "\n")

# Plot infection curve with 95% bootstrap CI 
# and the daily death data against day of the year
plot_final(t_obs, y, tf, fit_best, fhat_boot)