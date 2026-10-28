## Spanish energy prices case study
##
## 2-state Markov-switching GAMLSS for the Spanish energy price, with
## both the mean and the standard deviation of each state depending smoothly on
## the oil price. The model is fitted via qreml() with penalised splines, and
## the fitted state-dependent regression curves and decoded states are
## visualised. Finally, we refit the model by full REML and marginal ML (via
## hmmTMB) for comparison.

## Packages
# install.packages(c("MSwM", "RTMB, "LaMa", "hmmTMB", "scales"))
library(MSwM) # data
library(LaMa) # HMM functions and qreml
library(scales) # just for transparency in plotting
source("./case_studies/utils.R")

## Colors
color <- LaMaColors(2)

## Data
data(energy, package = "MSwM")
head(energy); nrow(energy)


# Fitting MS-GAMLSS -------------------------------------------------------

## Penalised likelihood function
pnll <- function(par) {
  getAll(par, dat) # makes all elements of par and dat available by name
  Gamma <- tpm(eta) # transition probability matrix

  ## Build state-dependent parameters from intercepts and spline coefficients
  beta <- cbind(beta0, betaspline); REPORT(beta) # mean parameter matrix
  alpha <- cbind(alpha0, alphaspline); REPORT(alpha) # sd parameter matrix
  Mu <- Z %*% t(beta) # mean, one column per state
  Sigma <- exp(Z %*% t(alpha)) # sd (log link), one column per state

  ## Compute state-dependent densities
  allprobs <- statedens(price, dnorm, Mu, Sigma)

  # penalty() does the heavy lifting and reports to qreml
  # (lambda has one entry per smooth: 2 for the mean, 2 for the sd)
  -forward("stationary", Gamma, allprobs) +
    penalty(list(betaspline, alphaspline), S, sp)
}

## Prepare model matrices (uses mgcv under the hood)
nb <- 15 # number of basis functions
modmat <- modmat(~ s(Oil, k = nb, bs = "ps"), energy) # model matrices Z and S
Z <- modmat$Z # design matrix
S <- modmat$S # penalty matrix (list)

## Prediction matrix for visualising state-dependent distributions later
xseq <- seq(min(energy$Oil), max(energy$Oil), length = 200) # sequence for prediction
Z_p <- predict(modmat, newdata = data.frame(Oil = xseq))

## Initial parameters
par <- list(eta = rep(-4, 2), # tpm parameters
            beta0 = c(2, 5), # mean intercepts
            betaspline = matrix(0, nrow = 2, ncol = nb-1), # one smooth per state
            alpha0 = c(0, 0), # sd intercepts (log scale)
            alphaspline = matrix(0, nrow = 2, ncol = nb-1))

## Data, model matrices and initial penalty strength
dat <- list(price = energy$Price, Z = Z, S = S,
            sp = rep(1e5, 4))

## Fit the model via qreml
t <- Sys.time()
mod <- qreml(pnll, par, dat, random = c("betaspline", "alphaspline"), # spline coefficients are random effects
             saveall = TRUE) # save all intermediate models for plotting
(time_qreml <- Sys.time() - t)

# estimated penalty strengths and edf
summary(mod)

## Decode states and predict state-dependent mean and sd on the grid
mod$states <- viterbi(mod) # most probable state sequence
Mu_p <- Z_p %*% t(mod$beta) # predicted means on grid
Sigma_p <- exp(Z_p %*% t(mod$alpha)) # predicted sds on grid


## Visualise results
# pdf("./case_studies/figs/energy_oil.pdf", width = 8, height = 4.5)
par(mfrow = c(1,2), mar = c(5,4,3,1))

# state-dependent mean and quantile curves on top of the data
plot(energy$Oil, energy$Price, pch = 20, bty = "n", col = alpha(color[mod$states], 0.1),
     xlab = "Oil price", ylab = "Energy price")
for(j in 1:2) lines(xseq, Mu_p[,j], col = color[j], lwd = 3)
qseq <- qnorm(seq(0.5, 0.95, length = 4)) # standard normal quantiles for the dashed curves
for(i in qseq){
  for(j in 1:2){
    lines(xseq, Mu_p[,j] + i*Sigma_p[,j], col = alpha(color[j], 0.7), lwd = 1, lty = 2)
    lines(xseq, Mu_p[,j] - i*Sigma_p[,j], col = alpha(color[j], 0.7), lwd = 1, lty = 2)
  }
}
legend("topright", bty = "n", legend = paste("State", 1:2), col = color, lwd = 3)

# decoded states along the time series (line segments coloured by state)
plot(NA, xlim = c(0, nrow(energy)), ylim = c(1,10), bty = "n",
     xlab = "Time", ylab = "Energy price")
segments(x0 = 1:(nrow(energy)-1), x1 = 2:nrow(energy),
         y0 = energy$Price[-nrow(energy)], y1 = energy$Price[-1],
         col = color[mod$states[-1]], lwd = 0.5)
# dev.off()


## Plot the model sequence of the qreml iterations
# shows how the fitted curves become smoother as lambda is updated
mods <- mod$allmods
plotind <- c(1,2,3,4,6,length(mods))

# pdf("./case_studies/figs/energy_oil_modseq.pdf", width = 8, height = 5.5)
par(mfrow = c(2,3), mar = c(5,4,3,1))
for(m in plotind){
  # decode states and predict curves for the m-th intermediate model
  states <- viterbi(mods[[m]])
  Mu_p <- Z_p %*% t(mods[[m]]$beta)
  Sigma_p <- exp(Z_p %*% t(mods[[m]]$alpha))
  plot(energy$Oil, energy$Price, pch = 20, bty = "n", col = alpha(color[states], 0.1),
       xlab = "Oil price", ylab = "Energy price", main = paste("Iteration", m-1))
  for(j in 1:2) lines(xseq, Mu_p[,j], col = color[j], lwd = 2)
  for(i in qseq){
    for(j in 1:2){
      lines(xseq, Mu_p[,j] + i*Sigma_p[,j], col = alpha(color[j], 0.7), lwd = 1, lty = 2)
      lines(xseq, Mu_p[,j] - i*Sigma_p[,j], col = alpha(color[j], 0.7), lwd = 1, lty = 2)
    }
  }
}
# dev.off()


# Full REML and marginal ML -----------------------------------------------

## Get the rank and log-determinant of the penalty matrix
# needed for proper REML likelihood evaluation (happens inside qreml() as well)
g <- LaMa:::gen_inverse(modmat$S[[1]])
rank_S <- g$rank
logdet_S <- g$logdet

## Joint likelihood (has complete normal density for random effects)
# here the smoothing parameters are estimated as part of the optimisation
jnll <- function(par) {
  getAll(par, dat)
  Gamma <- tpm(eta) # transition probability matrix

  ## Build state-dependent parameters from intercepts and spline coefficients
  beta <- cbind(beta0, betaspline); REPORT(beta) # mean parameter matrix
  alpha <- cbind(alpha0, alphaspline); REPORT(alpha) # sd parameter matrix
  Mu <- Z %*% t(beta) # mean
  Sigma <- exp(Z %*% t(alpha)) # sd

  ## Compute state-dependent densities
  allprobs <- statedens(price, dnorm, Mu, Sigma)
  nll <- -forward("stationary", Gamma, allprobs)

  ## Add negative log-density of the (improper) Gaussian random effects
  lambda <- exp(log_lambda); REPORT(lambda) # penalty strengths must be positive
  # mean coefficients use lambda[1:2]
  # each smooth: 0.5 * (lambda b'Sb - rank(S) log(lambda) - log|S|_+ + rank(S) log(2 pi))
  for(i in 1:2) {
    nll <- nll + 0.5 * (lambda[i] * sum(betaspline[i,] * (S %*% betaspline[i,])) -
                          rank_S * log(lambda[i]) - logdet_S + rank_S * log(2 * pi))
  }
  # sd coefficients use lambda[3:4]
  for(i in 1:2){
    j <- i + 2
    nll <- nll + 0.5 * (lambda[j] * sum(alphaspline[i,] * (S %*% alphaspline[i,])) -
                          rank_S * log(lambda[j]) - logdet_S + rank_S * log(2 * pi))
  }
  nll
}

## Initial parameters
par <- list(eta = rep(-3, 2),
            beta0 = c(2, 5),
            betaspline = matrix(0, nrow = 2, ncol = nb-1),
            alpha0 = c(0, 0),
            alphaspline = matrix(0, nrow = 2, ncol = nb-1),
            log_lambda = rep(log(1e5), 4))

## Data and model matrices
dat <- list(price = energy$Price, Z = Z, S = S[[1]]) # all smooths share the same penalty matrix

## Full REML fit
# everything except log_lambda is integrated out by the Laplace approximation
t <- Sys.time()
obj <- MakeADFun(jnll, par, random = names(par)[names(par) != "log_lambda"],
                 silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
(time_REML <- Sys.time() - t)

modREML <- lama(obj)
modREML$lambda # estimated smoothing parameters
edf_at(obj, jnll, random = c("betaspline", "alphaspline"),
       S = list(S[[1]], S[[1]])) # estimated effective degrees of freedom


## Marginal ML fit via hmmTMB
library(hmmTMB)

# state process: homogeneous 2-state Markov chain
hid <- MarkovChain$new(
  data = energy,
  n_states = 2
)

# observation process: spline effects of oil price on mean and sd
obs <- Observation$new(
  data = energy,
  formula = list(Price = list(
    mean = ~ s(Oil, k = nb, bs = "ps"),
    sd = ~ s(Oil, k = nb, bs = "ps"))),
  dists = list(Price = "norm"),
  par = list(Price = list(mean = c(2, 5), sd = c(1, 1))), # initial values
  n_states = 2
)

hmm <- HMM$new(obs = obs, hid = hid)

t <- Sys.time()
hmm$fit(silent = TRUE)
(time_hmmTMB <- Sys.time() - t)

hmm$lambda()$obs # estimated smoothing parameters

