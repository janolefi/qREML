## Spanish energy prices case study

# install.packages("RTMB")
# install.packages("MSwM")
library(MSwM) # data
# devtools::install_github("janoleko/LaMa") # development version
library(LaMa) # HMM functions and qreml
# install.packages("scales")
library(scales) # just for transparency in plotting
source("./case_studies/utils.R")


## colors
color = LaMaColors(2)

## data
data(energy, package = "MSwM")
head(energy); nrow(energy) 


# Fitting MS-GAMLSS -------------------------------------------------------

## penalised likelihood
pnll = function(par) {
  getAll(par, dat)
  Gamma = tpm(eta) # computing the tpm
  beta = cbind(beta0, betaspline); REPORT(beta) # mean parameter matrix
  alpha = cbind(alpha0, alphaspline);  REPORT(alpha) # sd parameter matrix
  Mu = Z %*% t(beta) # mean
  Sigma = exp(Z %*% t(alpha)) # sd
  allprobs = statedens(price, dnorm, Mu, Sigma) # state-dependent densities
  - forward("stationary", Gamma, allprobs) +
    penalty(list(betaspline, alphaspline), S, lambda) # penalty does the heavy lifting and reports to qreml
}

## Prepwork: model matrices (uses mgcv under the hood)
nb = 15 # number of basis functions
modmat = make_matrices(~ s(Oil, k = nb, bs = "ps"), energy) # model matrices Z and S
Z = modmat$Z # design matrix
S = modmat$S # penalty matrix (list)

## prediction matrix for visualizing state-dependent distributions later
xseq = seq(min(energy$Oil), max(energy$Oil), length = 200) # sequence for prediction
Z_p = predict(modmat, newdata = data.frame(Oil = xseq))

## initial parameter list
par = list(eta = rep(-4, 2),
           beta0 = c(2, 5),
           betaspline = matrix(0, nrow = 2, ncol = nb-1),
           alpha0 = c(0, 0),
           alphaspline = matrix(0, nrow = 2, ncol = nb-1))

## data, model matrices and initial penalty strength
dat = list(price = energy$Price, Z = Z, S = S, 
           lambda = rep(1e5, 4))

## model fit
t <- Sys.time()
mod <- qreml(pnll, par, dat, random = c("betaspline", "alphaspline"), 
             saveall = TRUE)
(time_qreml <- Sys.time() - t)

# Estimated penalty strengths and edf
summary(mod)

mod$states = viterbi(mod) # decoding most probable state sequence
Mu_p = Z_p %*% t(mod$beta) # predicted means on grid 
Sigma_p = exp(Z_p %*% t(mod$alpha)) # predicted sds on grid


## Visualising results
# pdf("./case_studies/figs/energy_oil.pdf", width = 8, height = 4.5)
par(mfrow = c(1,2), mar = c(5,4,3,1))
plot(energy$Oil, energy$Price, pch = 20, bty = "n", col = alpha(color[mod$states], 0.1),
     xlab = "Oil price", ylab = "Energy price")
for(j in 1:2) lines(xseq, Mu_p[,j], col = color[j], lwd = 3)
qseq = qnorm(seq(0.5, 0.95, length = 4))
for(i in qseq){
  for(j in 1:2){
    lines(xseq, Mu_p[,j] + i*Sigma_p[,j], col = alpha(color[j], 0.7), lwd = 1, lty = 2)
    lines(xseq, Mu_p[,j] - i*Sigma_p[,j], col = alpha(color[j], 0.7), lwd = 1, lty = 2)
  }
}
legend("topright", bty = "n", legend = paste("State", 1:2), col = color, lwd = 3)

plot(NA, xlim = c(0, nrow(energy)), ylim = c(1,10), bty = "n",
     xlab = "Time", ylab = "Energy price")
segments(x0 = 1:(nrow(energy)-1), x1 = 2:nrow(energy),
         y0 = energy$Price[-nrow(energy)], y1 = energy$Price[-1], col = color[mod$states[-1]], lwd = 0.5)
# dev.off()


## plotting the model sequence
mods = mod$allmods
plotind = c(1,2,3,4,6,length(mods))

# pdf("./case_studies/figs/energy_oil_modseq.pdf", width = 8, height = 5.5)
par(mfrow = c(2,3), mar = c(5,4,3,1))
for(m in plotind){
  states = viterbi(mods[[m]])
  Mu_p = Z_p %*% t( mods[[m]]$beta)
  Sigma_p = exp(Z_p %*% t(mods[[m]]$alpha))
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

# get the rank and log-determinant of the penalty matrix
# for proper REML likelihood evaluation (happens inside qreml() as well)
g <- LaMa:::gen_inverse(modmat$S[[1]])
rank_S <- g$rank
logdet_S <- g$logdet

## joint likelihood (has complete normal density for random effects)
jnll = function(par) {
  getAll(par, dat)
  Gamma = tpm(eta) # computing the tpm
  beta = cbind(beta0, betaspline); REPORT(beta) # mean parameter matrix
  alpha = cbind(alpha0, alphaspline); REPORT(alpha) # sd parameter matrix
  Mu = Z %*% t(beta) # mean
  Sigma = exp(Z %*% t(alpha)) # sd
  allprobs = statedens(price, dnorm, Mu, Sigma)
  nll = - forward("stationary", Gamma, allprobs) 

  # negative log-density of the (improper) Gaussian random effects
  lambda = exp(log_lambda); REPORT(lambda)
  for(i in 1:2) {
    nll <- nll + 0.5 * (lambda[i] * sum(betaspline[i,] * (S %*% betaspline[i,])) -
                          rank_S * log(lambda[i]) - logdet_S + rank_S * log(2 * pi))
  }
  for(i in 1:2){
    j = i + 2
    nll <- nll + 0.5 * (lambda[j] * sum(alphaspline[i,] * (S %*% alphaspline[i,])) -
                          rank_S * log(lambda[j]) - logdet_S + rank_S * log(2 * pi))
  }
  nll
}


## initial parameter list
par = list(eta = rep(-3, 2),
           beta0 = c(2, 5),
           betaspline = matrix(0, nrow = 2, ncol = nb-1),
           alpha0 = c(0, 0),
           alphaspline = matrix(0, nrow = 2, ncol = nb-1),
           log_lambda = rep(log(1e5), 4))

## data, model matrices and initial penalty strength
dat = list(price = energy$Price, Z = Z, S = S[[1]])

### REML
t <- Sys.time()
obj <- MakeADFun(jnll, par, random = names(par)[names(par)!="log_lambda"],
                 silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr) 
(time_REML <- Sys.time() - t)

modREML <- lama(obj)
modREML$lambda  # estimated smoothing parameters
# estimated effective degrees of freemdom
edf_at(obj, jnll, random = c("betaspline", "alphaspline"), S = list(S[[1]], S[[1]]))



### Marginal ML via hmmTMB

library(hmmTMB)

hid <- MarkovChain$new(
  data = energy,
  n_states = 2
)
obs <- Observation$new(
  data = energy,
  formula = list(Price = list(
    mean = ~ s(Oil, k = nb, bs = "ps"),
    sd = ~ s(Oil, k = nb, bs = "ps"))),
  dists = list(Price = "norm"),
  par = list(Price = list(mean = c(2, 5), sd = c(1, 1))),
  n_states = 2,
)

hmm <- HMM$new(obs = obs, hid = hid)

t <- Sys.time()
hmm$fit(silent = TRUE)
(time_hmmTMB <- Sys.time() - t)

hmm$lambda()$obs
