## Ivory Coast elephant case study
##
## 2-state HMM for elephant movement (step length and turning angle)
## with a time-of-day (tod) dependent state process. We first fit a homogeneous
## HMM and an HMM with periodic (trigonometric) transition probabilities.
## We then model the transition probabilities as smooth, cyclic functions of
## time of day with qreml(), including confidence bands from Monte Carlo
## samples. Finally, we refit the smooth model by full REML and marginal ML
## (via hmmTMB) for comparison.

## Packages
# install.packages(c("RTMB, "LaMa", "hmmTMB", "scales"))
library(LaMa)
library(scales) # for transparent colors
source("./case_studies/utils.R")

## Colors
color <- LaMaColors(2)

## Data
data <- read.csv("./data/elephant_data.csv")
data$step[data$step == 0] <- NA # 2 zero step lengths -> set to missing instead of zigamma
head(data); nrow(data)


# Fitting homogeneous and parametric HMMs ---------------------------------

## Likelihood function for the homogeneous HMM
nll_hom <- function(par){
  getAll(par, dat) # makes all elements of par and dat available by name

  ## Transform parameters
  mu <- exp(log_mu); REPORT(mu) # mean step length
  sigma <- exp(log_sigma); REPORT(sigma) # sd of step length
  kappa <- exp(log_kappa); REPORT(kappa) # concentration of turning angle
  Gamma <- tpm(eta) # transition probability matrix from unconstrained parameters

  ## Compute state-dependent densities
  # step and angle are assumed conditionally independent given the state
  allprobs <- statedens(step, dgamma2, mu, sigma) * statedens(angle, dvm, 0, kappa)

  -forward("stationary", Gamma, allprobs) # initial distribution is stationary
}

par <- list(log_mu = log(c(0.2, 2)),
            log_sigma = log(c(0.2, 2)),
            log_kappa = log(c(0.2, 1)),
            eta = rep(-2, 2))
dat <- list(step = data$step, angle = data$angle, N = 2)

obj_hom <- MakeADFun(nll_hom, par, silent = TRUE) # build the AD objective function
opt_hom <- nlminb(obj_hom$par, obj_hom$fn, obj_hom$gr) # numerical maximum likelihood
mod_hom <- lama(obj_hom) # collect estimates and reported quantities
summary(mod_hom)

## Plot the marginal and state-dependent distributions
# pdf("./case_studies/figs/elephant_marginal.pdf", width = 8, height = 4)
par(mfrow = c(1,2))
plot(mod_hom, "statedens", stream = "step", n = 2000, legend = FALSE,
     breaks = 100, xlim = c(0,5), xlab = "Step length (km)")
legend("topright", legend = c("Encamped", "Exploratory", "Marginal"),
       col = c(color[1], color[2], "black"), lty = c(1,1,2), bty = "n")
plot(mod_hom, "statedens", stream = "angle", n = 500, legend = FALSE,
     xlab = "Turning angle (radians)")
# dev.off()


## Likelihood function with periodic variation in the state process
nll_par <- function(par){
  getAll(par, dat)
  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  kappa <- exp(log_kappa); REPORT(kappa)
  Gamma <- tpm(beta, Z) # inhomogeneous transition probability matrix, one per time point
  allprobs <- statedens(step, dgamma2, mu, sigma) * statedens(angle, dvm, 0, kappa)
  -forward("stationary", Gamma, allprobs, tod = tod) # tod picks the tpm for each observation
}

# build the trigonometric basis design matrix
# tod indexes two-hour blocks, (1:12)*2-1 are the corresponding hours of the day
Z <- cosinor((1:12)*2-1, period = c(24, 12))

par <- list(log_mu = log(c(0.2, 2)),
            log_sigma = log(c(0.2, 2)),
            log_kappa = log(c(0.2, 1)),
            beta = matrix(c(rep(-2,2), rep(0, 8)), nrow = 2)) # intercepts, then 4 trig. coefficients
dat <- list(step = data$step, angle = data$angle, N = 2, tod = data$tod)

obj_par <- MakeADFun(nll_par, par, silent = TRUE)
opt_par <- nlminb(obj_par$par, obj_par$fn, obj_par$gr)
mod_par <- lama(obj_par)
mod_par$Delta <- stationary_p(mod_par$Gamma) # periodically stationary distribution
summary(mod_par)


## Confidence intervals by sampling from the MLE distribution
# MCreport() draws parameters from the asymptotic normal distribution of the MLE
# and evaluates the REPORT()ed quantities for each draw
nSamples <- 10000
samples <- MCreport(obj_par, nSamples = nSamples)

n <- 200 # plotting resolution
tod_seq <- seq(0, 24, length = n) # sequence of time of day values
Z_p <- cosinor(tod_seq, period = c(24, 12)) # prediction design matrix

Gamma_samples <- array(dim = c(2, 2, n, nSamples))
for(b in 1:nSamples) Gamma_samples[,,,b] <- tpm(samples$beta[[b]], Z_p)
GammaCI_par <- apply(Gamma_samples, c(1,2,3), quantile, probs = c(0.025, 0.975)) # pointwise 95% CI
Gamma_plot_par <- tpm(mod_par$beta, Z_p) # point estimate


# Transition probabilities as smooth functions ----------------------------

## Penalised likelihood function
pnll <- function(par) {
  getAll(par, dat)
  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  kappa <- exp(log_kappa); REPORT(kappa)
  Gamma <- tpm(cbind(beta0, betaspline), Z) # intercepts plus spline coefficients
  allprobs <- statedens(step, dgamma2, mu, sigma) * statedens(angle, dvm, 0, kappa)
  -forward("stationary", Gamma, allprobs, tod = tod) +
    penalty(betaspline, S, lambda) # 0.5 * lambda t(b) %*% S %*% b for each b
}

## Prepare model matrices (cyclic P-spline basis)
nb <- 10 # number of basis functions
modmat <- make_matrices(~ s(tod, bs = "cp", k = nb),
                        data = data.frame(tod = (1:12)*2-1),
                        knots = list(tod = c(0, 24))) # where to wrap the basis
Z <- modmat$Z # spline design matrix
S <- modmat$S # spline penalty matrix (list)

## Initial parameters and data
N <- 2
par <- list(log_mu = log(c(0.35, 1.1)),
            log_sigma = log(c(0.25, 0.75)),
            log_kappa = log(c(0.2, 0.7)),
            beta0 = c(-2,-2), # intercepts of the two off-diagonal tpm entries
            betaspline = matrix(0, nrow = N*(N-1), ncol = nb-1)) # one smooth per row

dat <- list(step = data$step, angle = data$angle,
            tod = data$tod,
            N = 2, Z = Z, S = S,
            lambda = rep(1e5, 2)) # initial smoothing parameters, one per smooth

## Fit the model via qreml
t <- Sys.time()
mod <- qreml(pnll, par, dat, random = "betaspline", # spline coefficients are random effects
             saveall = TRUE) # save all intermediate models for plotting
(time_qreml <- Sys.time() - t)

summary(mod)

## Sample reported quantities from the MLE distribution
nSamples <- 5000
samples <- MCreport(mod$obj, nSamples = nSamples)
beta_samples <- lapply(1:nSamples, function(i){
  cbind(samples$beta0[[i]], samples$betaspline[[i]]) # same layout as in pnll()
})


## Calculate transition probabilities and stationary distribution for plotting
Z_p <- predict(modmat, newdata = data.frame(tod = tod_seq)) # prediction design matrix
Gamma_p <- tpm(mod$beta, Z_p)

# periodically stationary distribution on the prediction grid
Delta_p <- matrix(NA, length(tod_seq), N)
for(t in 1:length(tod_seq)){
  # 24 hours of tpms, starting at the current grid point
  t_seq <- (tod_seq[t] + (1:12)*2 - 1) %% 24
  Z_cont <- predict(modmat, newdata = data.frame(tod = t_seq))
  Delta_p[t,] <- stationary_p(tpm(mod$beta, Z_cont), t = 1) # distribution at the start of the cycle
}


## Compute confidence bands from Monte Carlo samples
Gamma_samples <- array(dim = c(2, 2, length(tod_seq), nSamples))
Delta_samples <- array(dim = c(length(tod_seq), 2, nSamples))
for(b in 1:nSamples) {
  Gamma_samples[,,,b] <- tpm(beta_samples[[b]], Z_p)
}

## Periodically stationary distribution uncertainty on the prediction grid
# this is slow: it computes the periodically stationary distribution for each
# grid point and each sample of the parameters
for(t in 1:length(tod_seq)) {
  # if(t %% 10 == 0) print(paste0(round(100*t/length(tod_seq)), "%"))
  t_seq <- (tod_seq[t] + (1:12)*2 - 1) %% 24
  Z_t <- predict(modmat, newdata = data.frame(tod = t_seq))
  for(b in 1:nSamples){
    Delta_samples[t,,b] <- stationary_p(tpm(beta_samples[[b]], Z_t), t = 1)
  }
}

# pointwise 95% confidence intervals
GammaCI <- apply(Gamma_samples, c(1,2,3), quantile, probs = c(0.025, 0.975), na.rm = TRUE)
DeltaCI <- apply(Delta_samples, c(1,2), quantile, probs = c(0.025, 0.975))


## Visualise results: parametric vs. nonparametric transition probabilities
# pdf("./case_studies/figs/elephant_transprobs.pdf", width = 8, height = 4)
par(mfrow = c(1,2))

# parametric fit (gamma_12 solid, gamma_21 dashed, each with confidence band)
plot(tod_seq, Gamma_plot_par[1,2,], type = "l", lwd = 1, bty = "n", main = "Parametric",
     xlab = "Time of day", ylab = "Transition probability", xaxt = "n", ylim = c(0,1))
polygon(c(tod_seq, rev(tod_seq)), c(GammaCI_par[1,1,2,], rev(GammaCI_par[2,1,2,])),
        col = alpha("black", 0.1), border = FALSE)
lines(tod_seq, Gamma_plot_par[2,1,], lwd = 1, lty = 5)
polygon(c(tod_seq, rev(tod_seq)), c(GammaCI_par[1,2,1,], rev(GammaCI_par[2,2,1,])),
        col = alpha("black", 0.1), border = FALSE)
legend(x = -0.5, y = 1.02, lty = c(1,5), y.intersp = 1.3,
       legend = c(expression(gamma[12]^(t)), expression(gamma[21]^(t))), bty = "n")
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

# nonparametric fit
plot(tod_seq, Gamma_p[1,2,], type = "l", lwd = 1, bty = "n", main = "Nonparametric",
     xlab = "Time of day", ylab = "Transition probability", xaxt = "n", ylim = c(0,1))
polygon(c(tod_seq, rev(tod_seq)), c(GammaCI[1,1,2,], rev(GammaCI[2,1,2,])),
        col = alpha("black", 0.1), border = FALSE)
lines(tod_seq, Gamma_p[2,1,], lwd = 1, lty = 5)
polygon(c(tod_seq, rev(tod_seq)), c(GammaCI[1,2,1,], rev(GammaCI[2,2,1,])),
        col = alpha("black", 0.1), border = FALSE)
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

# dev.off()


## Visualise the periodically stationary distribution
# pdf("./case_studies/figs/elephant_stationary.pdf", width = 7, height = 4.5)
par(mfrow = c(1,1), mar = c(5,4,0,2)+0.1)
plot(NA, bty = "n", ylim = c(0,1), xlim = c(0,24), ylab = "Pr(exploratory)",
     xlab = "Time of day", xaxt = "n", cex = 1.5)

# shade the night
polygon(x = c(0, 6.5, 6.5, 0), y = c(0, 0, 1, 1), col = "gray95", border = "white")
polygon(x = c(18.5, 24, 24, 18.5), y = c(0, 0, 1, 1), col = "gray95", border = "white")

# colour bar showing the sun cycle (in half-hour steps)
for(t in 0:47){
  polygon(x = c(t/2, (t+1)/2, (t+1)/2, t/2), y = c(-0.04, -0.04, -0.01, -0.01),
          col = sun_cycle_colors[t+1], border = sun_cycle_colors[t+1])
}

# estimate with confidence band
lines(tod_seq, Delta_p[,2], lwd = 1.5, col = "black")
polygon(c(tod_seq, rev(tod_seq)), c(DeltaCI[1,,2], rev(DeltaCI[2,,2])),
        col = alpha("black", 0.1), border = FALSE)
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))
# dev.off()


## Plot the model sequence of the qreml iterations
# shows how the estimated curves become smoother as lambda is updated
mods <- mod$allmods
plotind <- c(1,2,3,5,7,length(mods))

# pdf("./case_studies/figs/elephant_modseq.pdf", width = 8.5, height = 5.5)
par(mfrow = c(2,3), mar = c(5,4,3,1))
for(m in plotind){
  beta <- mods[[m]]$beta
  Gamma_p <- tpm(beta, Z_p)
  plot(tod_seq, Gamma_p[1,2,], type = "l", lwd = 1, bty = "n", main = paste("Iteration", m-1),
       xlab = "Time of day", ylab = "Transition probability", xaxt = "n", ylim = c(0,1))
  lines(tod_seq, Gamma_p[2,1,], lwd = 1, lty = 5)
  legend(x = -0.5, y = 1.02, lty = c(1,2), y.intersp = 1.3,
         legend = c(expression(gamma[12]^(t)), expression(gamma[21]^(t))), bty = "n")
  axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))
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
  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  kappa <- exp(log_kappa); REPORT(kappa)
  Gamma <- tpm(cbind(beta0, betaspline), Z)
  allprobs <- statedens(step, dgamma2, mu, sigma) * statedens(angle, dvm, 0, kappa)
  nll <- -forward("stationary", Gamma, allprobs, tod = tod)

  ## Add negative log-density of the (improper) Gaussian random effects
  lambda <- exp(log_lambda); REPORT(lambda) # penalty strengths must be positive
  for(i in 1:2){
    # one penalty per smooth: 0.5 * (lambda b'Sb - rank(S) log(lambda) - log|S|_+ + rank(S) log(2 pi))
    nll <- nll + 0.5 * (lambda[i] * sum(betaspline[i,] * (S %*% betaspline[i,])) -
                          rank_S * log(lambda[i]) - logdet_S + rank_S * log(2 * pi))
  }
  nll
}

## Initial parameters and data
par <- list(log_mu = log(c(0.35, 1.1)),
            log_sigma = log(c(0.25, 0.75)),
            log_kappa = log(c(0.2, 0.7)),
            beta0 = c(-2,-2),
            betaspline = matrix(0, nrow = N*(N-1), ncol = ncol(Z)-1),
            log_lambda = rep(log(1e5), 2))
dat <- list(step = data$step, angle = data$angle,
            tod = data$tod,
            N = 2, Z = Z, S = S[[1]]) # both smooths share the same penalty matrix

## Full REML fit
# everything except log_lambda is integrated out by the Laplace approximation
t <- Sys.time()
obj1 <- MakeADFun(jnll, par, random = names(par)[names(par) != "log_lambda"],
                  silent = TRUE)
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
(time_reml <- Sys.time() - t)

mod_reml <- lama(obj1)
mod_reml$lambda # estimated smoothing parameters
edf_at(obj1, jnll, random = "betaspline", S = dat$S) # estimated effective degrees of freedom


## Marginal ML fit via hmmTMB
library(hmmTMB)

data$hour <- data$tod*2 - 1 # hmmTMB needs the covariate on the hour scale of the knots

# state process: cyclic spline in time of day
hid <- MarkovChain$new(
  data = data,
  formula = ~ s(hour, bs = "cp", k = nb),
  n_states = 2,
  gam_args = list(knots = list(hour = c(0, 24)))
)

# observation process: same distributions and initial values as above
obs <- Observation$new(
  data = data,
  dists = list(step = "gamma2", angle = "vm"),
  par = list(step = list(mean = c(0.35, 1.1), sd = c(0.25, 0.75)),
             angle = list(mu = c(0,0), kappa = c(0.2, 0.7))),
  n_states = 2,
  fixpar = list(obs = list(angle = list(mu = rep(NA, 2)))) # fix the mean angle at 0
)

hmm <- HMM$new(obs = obs, hid = hid)

t <- Sys.time()
hmm$fit(silent = TRUE)
(time_hmmTMB <- Sys.time() - t)
