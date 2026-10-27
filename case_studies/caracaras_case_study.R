## Caracaras case study

## packages
# install.packages("RTMB")
# install.packages("LaMa")
library(LaMa)
# install.packages("scales")
library(scales) # for transparent colors
source("./case_studies/utils.R")


## color palette
color = LaMaColors(3)

## data
data = read.csv("./data/caracara_downsampled.csv")
head(data); nrow(data)


# Simple parametric HMMs --------------------------------------------------

## likelihood function for simple homogeneous HMM
nll = function(par){
  getAll(par, dat)
  sigma = exp(log_sigma); REPORT(mu); REPORT(sigma)
  Gamma = tpm(eta)
  allprobs = statedens(x, dnorm, mu, sigma)
  -forward("stationary", Gamma, allprobs)
}

## First fit the 3 state parametric model
par = list(mu = c(-5, -4, -2), # state-dep mean
           log_sigma = rep(log(0.3), 3), # state-dep sd
           eta = rep(-2, 6)) # tpm params
dat = list(x = data$logVDBA, N = 3)

## fitting model
obj3 = MakeADFun(nll, par, silent = TRUE)
opt3 = nlminb(obj3$par, obj3$fn, obj3$gr)

## extracting parameters
mod3 = lama(obj3)
summary(mod3)

## plotting estimated state-dependent distributions
# pdf("./case_studies/figs/caracaras_simple.pdf", width = 8.5, height = 4)
par(mfrow = c(1,2))
plot(mod3, "statedens", legend = FALSE, n = 300, xlab = "log(VeDBA)")
legend("topright", legend = c("Resting", "Low activity", "High activity"), 
       col = color, lwd = 2, bty = "n")



## Then fit the 4 state parametric model
par = list(mu = c(-5, -4.5, -3, -1.5), # state-dep mean
           log_sigma = rep(log(0.3), 4), # state-dep sd
           eta = rep(-2, 12)) # tpm params
dat = list(x = data$logVDBA, N = 4)

## fitting model
obj4 = MakeADFun(nll, par, silent = TRUE)
opt4 = nlminb(obj4$par, obj4$fn, obj4$gr)

## extracting parameters
mod4 = lama(obj4)
summary(mod4)

## plotting estimated state-dependent distributions
color2 = c("#F0E442", color)
plot(mod4, "statedens", legend = FALSE, n = 300, col = color2, xlab = "log(VeDBA)")
legend("topright", legend = c("Resting 1", "Resting 2", "Low activity", "High activity"), 
       col = color2, lwd = 2, bty = "n")
# dev.off()


## Information criteria
AIC(mod3, mod4)
BIC(mod3, mod4)



# Model with nonparametric emission distributions ------------------------------

## penalised likelihood function
pnll = function(par) {
  getAll(par, dat)
  Gamma = tpm(eta) # transition probability matrix
  ## computing mixture weights
  alpha = exp(cbind(beta, rep(0, N))) # last column fixed to zero
  alpha = alpha / rowSums(alpha); REPORT(alpha) # multinomial logit
  ## computing state-dependent distributions
  allprobs = matrix(1, nrow(Z), N)
  ind = which(!is.na(x))
  allprobs[ind,] = Z[ind,] %*% t(alpha) # sum over weighted basis functions
  -forward("stationary", Gamma, allprobs) + penalty(beta, S, lambda)
}


## prepwork: model matrices and initial coefficients for smooth densities
N = 3
smoothDens = smooth_dens_construct(
  data["logVDBA"], 
  par = list(logVDBA = list(mean = c(-5, -4, -2.5), sd = rep(1.5, 3))),
  k = 30
)

Z = smoothDens$Z$logVDBA # design matrix
S = smoothDens$S$logVDBA # penalty matrix
Z_p = smoothDens$Z_predict$logVDBA # prediction design matrix
xseq = smoothDens$xseq$logVDBA # prediction grid
beta = smoothDens$coef$logVDBA # initial coefficients for basis functions

par = list(eta = rep(-3, N*(N-1)),beta = beta)
dat = list(x = data$logVDBA, 
           Z = Z, N = 3, S = S, 
           lambda = rep(100, 3))

## model fitting via qreml
t <- Sys.time()
mod <- qreml(pnll, par, dat, random = "beta")
(time_qreml <- Sys.time() - t)

summary(mod)
round(mod$lambda, 2) # estimated smoothing parameters

## extracting parameters
alpha = mod$alpha # weight matrix
delta = mod$delta # stationary distribution
mod$states = viterbi(mod) # state decoding

## AIC and BIC
AIC(mod, mod4, mod3)
BIC(mod, mod4, mod3)


## visualizing results
# pdf("./case_studies/figs/caracara_log.pdf", width = 8.5, height = 4)

Dens = Z_p %*% t(alpha) # state-dependent densities
mDens = rowSums(Dens * matrix(delta, nrow(Dens), N, byrow = TRUE)) # marginal density
par(mfrow = c(1,2), mar = c(5,4,3.5,1) + 0.1)
hist(data$logVDBA, prob = T, breaks = 50, bor = "white",
     main = "", xlab = "log(VeDBA)", ylab = "Density")
for(j in 1:N) lines(xseq, delta[j] * Dens[,j], col = color[j], lwd = 2)
lines(xseq, mDens, col = 1, lty = 2, lwd = 2)

legend("topright", legend = c("Resting", "Low activity", "High activity"), 
       col = color, lwd = 2, bty = "n")

plotind = 4000:8000
plot(data$logVDBA[plotind], pch = 16, cex = 0.8, col = alpha(color[mod$states[plotind]], 0.3), 
     xlab = "Time", ylab = "log(VeDBA)", bty = "n")

# dev.off()





# Full REML and marginal ML -----------------------------------------------

# get the rank and log-determinant of the penalty matrix
# for proper REML likelihood evaluation (happens inside qreml() as well)
g <- LaMa:::gen_inverse(S)
rank_S <- g$rank
logdet_S <- g$logdet

jnll = function(par) {
  getAll(par, dat)
  Gamma = tpm(eta) # transition probability matrix

  ## computing mixture weights
  alpha = exp(cbind(beta, rep(0, N))) # last column fixed to zero
  alpha = alpha / rowSums(alpha) # multinomial logit link
  REPORT(alpha) # report alpha for convenience later
  
  ## computing state-dependent distributions
  allprobs = matrix(1, nrow(Z), N)
  ind = which(!is.na(x))
  allprobs[ind,] = Z[ind,] %*% t(alpha) # linear: sum over weighted basis functions
  nll = -forward("stationary", Gamma, allprobs) 
  
  # negative log-density of the (improper) Gaussian random effects
  lambda = exp(log_lambda); REPORT(lambda)
  for(i in 1:N){
    nll <- nll + 0.5 * (lambda[i] * sum(beta[i,] * (S %*% beta[i,])) -
                          rank_S * log(lambda[i]) - logdet_S + rank_S * log(2 * pi))
  }
  nll
}

par = list(eta = rep(-3, N*(N-1)), 
           beta = beta, 
           log_lambda = rep(log(100), 3))
dat = list(x = data$logVDBA, Z = Z, N = 3, S = S)

## REML
t <- Sys.time()
obj_reml <- MakeADFun(jnll, par, random = names(par)[names(par)!="log_lambda"],
                      silent = TRUE)
opt_reml <- nlminb(obj_reml$par, obj_reml$fn, obj_reml$gr)
(time_REML <- Sys.time() - t)

mod_reml <- lama(obj_reml)
mod_reml$lambda
edf_at(obj_reml, jnll, random = "beta", S = S) 


## marginal ML
t <- Sys.time()
obj_ml <- MakeADFun(jnll, par, random = "beta", silent = TRUE)
opt_ml <- nlminb(obj_ml$par, obj_ml$fn, obj_ml$gr)
(time_ML <- Sys.time() - t)

mod_ml <- lama(obj_ml)
mod_ml$lambda
edf_at(obj_ml, jnll, random = "beta", S = S)