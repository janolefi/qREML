## Caracaras case study
##
## HMMs for the activity (log VeDBA) of a caracara, derived from
## accelerometer data. We first fit simple parametric HMMs with 3 and 4 states
## and compare them via AIC/BIC. We then fit a 3-state HMM with nonparametric
## (spline-based) state-dependent densities using qreml(), and finally refit
## the same model by full REML and marginal ML for comparison.

## Packages
# install.packages(c("RTMB", "LaMa", "scales"))
library(LaMa)
library(scales) # for transparent colors
source("./case_studies/utils.R")


## Colors
color <- LaMaColors(3)

## Data
data <- read.csv("./data/caracara_downsampled.csv")
head(data); nrow(data)


# Simple parametric HMMs --------------------------------------------------

## Likelihood function for a simple homogeneous HMM
nll <- function(par){
  getAll(par, dat) # makes all elements of par and dat available by name
  sigma <- exp(log_sigma); REPORT(mu); REPORT(sigma) # sd must be positive
  Gamma <- tpm(eta) # transition probability matrix from unconstrained parameters
  allprobs <- statedens(x, dnorm, mu, sigma) # n x N matrix of state-dependent densities
  -forward("stationary", Gamma, allprobs) # initial distribution is stationary
}

## Fit the 3-state parametric model
par <- list(mu = c(-5, -4, -2), # state-dependent mean
            log_sigma = rep(log(0.3), 3), # state-dependent sd
            eta = rep(-2, 6)) # tpm parameters (off-diagonal, multinomial logit)
dat <- list(x = data$logVDBA, N = 3)

obj3 <- MakeADFun(nll, par, silent = TRUE) # build the AD objective function
opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr) # numerical maximum likelihood

mod3 <- lama(obj3) # collect estimates and reported quantities
summary(mod3)

## Plot the estimated state-dependent distributions
# pdf("./case_studies/figs/caracaras_simple.pdf", width = 8.5, height = 4)
par(mfrow = c(1,2))
plot(mod3, "statedens", legend = FALSE, n = 300, xlab = "log(VeDBA)")
legend("topright", legend = c("Resting", "Low activity", "High activity"),
       col = color, lwd = 2, bty = "n")


## Fit the 4-state parametric model
par <- list(mu = c(-5, -4.5, -3, -1.5), # state-dependent mean
            log_sigma = rep(log(0.3), 4), # state-dependent sd
            eta = rep(-2, 12)) # tpm parameters
dat <- list(x = data$logVDBA, N = 4)

obj4 <- MakeADFun(nll, par, silent = TRUE)
opt4 <- nlminb(obj4$par, obj4$fn, obj4$gr)

mod4 <- lama(obj4)
summary(mod4)

## Plot the estimated state-dependent distributions
color2 <- c("#F0E442", color) # extra color for the additional resting state
plot(mod4, "statedens", legend = FALSE, n = 300, col = color2, xlab = "log(VeDBA)")
legend("topright", legend = c("Resting 1", "Resting 2", "Low activity", "High activity"),
       col = color2, lwd = 2, bty = "n")
# dev.off()


## Information criteria for the parametric models
AIC(mod3, mod4)
BIC(mod3, mod4)


# Model with nonparametric emission distributions -------------------------

## Penalised likelihood function
# each state-dependent density is a weighted sum of B-spline basis densities
pnll <- function(par) {
  getAll(par, dat)
  Gamma <- tpm(eta) # transition probability matrix

  ## Compute mixture weights (multinomial logit)
  alpha <- exp(cbind(beta, rep(0, N))) # last column fixed to zero
  alpha <- alpha / rowSums(alpha); REPORT(alpha)

  ## Compute state-dependent densities
  allprobs <- matrix(1, nrow(Z), N) # missing observations keep probability 1
  ind <- which(!is.na(x))
  allprobs[ind,] <- Z[ind,] %*% t(alpha) # sum over weighted basis functions

  # penalty() adds the smoothness penalty on beta and reports it to qreml()
  -forward("stationary", Gamma, allprobs) + penalty(beta, S, lambda)
}


## Prepare model matrices and initial coefficients for the smooth densities
N <- 3
smoothDens <- smooth_dens_construct(
  data["logVDBA"],
  par = list(logVDBA = list(mean = c(-5, -4, -2.5), sd = rep(1.5, 3))), # initial densities
  k = 30 # number of basis functions
)

Z <- smoothDens$Z$logVDBA # design matrix
S <- smoothDens$S$logVDBA # penalty matrix
Z_p <- smoothDens$Z_predict$logVDBA # prediction design matrix
xseq <- smoothDens$xseq$logVDBA # prediction grid
beta <- smoothDens$coef$logVDBA # initial coefficients for basis functions

par <- list(eta = rep(-3, N*(N-1)), beta = beta)
dat <- list(x = data$logVDBA,
            Z = Z, N = 3, S = S,
            lambda = rep(100, 3)) # initial smoothing parameters, one per state

## Fit the model via qreml
t <- Sys.time()
mod <- qreml(pnll, par, dat, random = "beta") # beta is treated as a random effect
(time_qreml <- Sys.time() - t)

summary(mod)
round(mod$lambda, 2) # estimated smoothing parameters

## Extract parameters and decode states
alpha <- mod$alpha # weight matrix
delta <- mod$delta # stationary distribution
mod$states <- viterbi(mod) # most likely state sequence

## Information criteria, comparing all three models
AIC(mod, mod4, mod3)
BIC(mod, mod4, mod3)


## Visualise results
# pdf("./case_studies/figs/caracara_log.pdf", width = 8.5, height = 4)

Dens <- Z_p %*% t(alpha) # state-dependent densities on the prediction grid
mDens <- rowSums(Dens * matrix(delta, nrow(Dens), N, byrow = TRUE)) # marginal density
par(mfrow = c(1,2), mar = c(5,4,3.5,1) + 0.1)

# state-dependent (weighted) and marginal densities on top of the data histogram
hist(data$logVDBA, prob = TRUE, breaks = 50, bor = "white",
     main = "", xlab = "log(VeDBA)", ylab = "Density")
for(j in 1:N) lines(xseq, delta[j] * Dens[,j], col = color[j], lwd = 2)
lines(xseq, mDens, col = 1, lty = 2, lwd = 2)

legend("topright", legend = c("Resting", "Low activity", "High activity"),
       col = color, lwd = 2, bty = "n")

# decoded states along a section of the time series
plotind <- 4000:8000
plot(data$logVDBA[plotind], pch = 16, cex = 0.8, col = alpha(color[mod$states[plotind]], 0.3),
     xlab = "Time", ylab = "log(VeDBA)", bty = "n")

# dev.off()


# Full REML and marginal ML -----------------------------------------------

## Get the rank and log-determinant of the penalty matrix
# needed for proper REML likelihood evaluation (happens inside qreml() as well)
g <- LaMa:::gen_inverse(S)
rank_S <- g$rank
logdet_S <- g$logdet

## Joint likelihood (has complete normal density for random effects)
# here the smoothing parameters are estimated as part of the optimisation
jnll <- function(par) {
  getAll(par, dat)
  Gamma <- tpm(eta) # transition probability matrix

  ## Compute mixture weights (multinomial logit)
  alpha <- exp(cbind(beta, rep(0, N))) # last column fixed to zero
  alpha <- alpha / rowSums(alpha)
  REPORT(alpha) # report alpha for convenience later

  ## Compute state-dependent densities
  allprobs <- matrix(1, nrow(Z), N)
  ind <- which(!is.na(x))
  allprobs[ind,] <- Z[ind,] %*% t(alpha) # sum over weighted basis functions
  nll <- -forward("stationary", Gamma, allprobs)

  ## Add negative log-density of the (improper) Gaussian random effects
  lambda <- exp(log_lambda); REPORT(lambda) # penalty strengths must be positive
  for(i in 1:N){
    # one penalty per state: 0.5 * (lambda b'Sb - rank(S) log(lambda) - log|S|_+ + rank(S) log(2 pi))
    nll <- nll + 0.5 * (lambda[i] * sum(beta[i,] * (S %*% beta[i,])) -
                          rank_S * log(lambda[i]) - logdet_S + rank_S * log(2 * pi))
  }
  nll
}

## Initial parameters and data
par <- list(eta = rep(-3, N*(N-1)),
            beta = beta,
            log_lambda = rep(log(100), 3))
dat <- list(x = data$logVDBA, Z = Z, N = 3, S = S)

## Full REML fit
# everything except log_lambda is integrated out by the Laplace approximation
t <- Sys.time()
obj_reml <- MakeADFun(jnll, par, random = names(par)[names(par) != "log_lambda"],
                      silent = TRUE)
opt_reml <- nlminb(obj_reml$par, obj_reml$fn, obj_reml$gr)
(time_REML <- Sys.time() - t)

mod_reml <- lama(obj_reml)
mod_reml$lambda # estimated smoothing parameters
edf_at(obj_reml, jnll, random = "beta", S = S) # estimated effective degrees of freedom


## Marginal ML fit
# only beta is integrated out, the remaining parameters are estimated by ML
t <- Sys.time()
obj_ml <- MakeADFun(jnll, par, random = "beta", silent = TRUE)
opt_ml <- nlminb(obj_ml$par, obj_ml$fn, obj_ml$gr)
(time_ML <- Sys.time() - t)

mod_ml <- lama(obj_ml)
mod_ml$lambda # estimated smoothing parameters
edf_at(obj_ml, jnll, random = "beta", S = S) # estimated effective degrees of freedom
