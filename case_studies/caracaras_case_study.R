## Caracaras case study

## packages
# install.packages("RTMB")
library(RTMB)
# install.packages("LaMa")
library(LaMa)
# install.packages("scales")
library(scales) # for transparent colors

## data
data = read.csv("./data/caracara_downsampled.csv")
head(data)
nrow(data)

## color palette
color = LaMaColors(3)


# Simple parametric HMM ---------------------------------------------------

## likelihood function for simple homogeneous HMM
nll = function(par){
  getAll(par, dat)
  sigma = exp(log_sigma); REPORT(mu); REPORT(sigma)
  Gamma = tpm(eta)
  delta = stationary(Gamma)
  allprobs = matrix(1, length(x), N)
  ind = which(!is.na(x))
  for(j in 1:N) allprobs[ind, j] = dnorm(x[ind], mu[j], sigma[j])
  -forward(delta, Gamma, allprobs)
}

## First fit the 3 state parametric model
N = 3
par = list(mu = c(-5, -4, -2), # state-dep mean
           log_sigma = rep(log(0.3), 3), # state-dep sd
           eta = rep(-2, 6)) # tpm params
dat = list(x = data$logVDBA, N = 3)

## fitting model
obj3 = MakeADFun(nll, par, silent = TRUE)
opt3 = nlminb(obj3$par, obj3$fn, obj3$gr)

## extracting parameters
mod3 = report(obj3)
delta3 = mod3$delta
mu3 = mod3$mu
sigma3 = mod3$sigma

## plotting estimated state-dependent distributions

# pdf("./case_studies/figs/caracaras_simple.pdf", width = 8.5, height = 4)

par(mfrow = c(1,2))

hist(data$logVDBA, breaks = 50, prob = T, bor = "white", main = "", xlab = "log(VeDBA)", ylab = "Density")
for(j in 1:N) curve(delta3[j] * dnorm(x, mu3[j], sigma3[j]), add = T, col = color[j], lwd = 2, n = 200)
curve(delta3[1] * dnorm(x, mu3[1], sigma3[1]) + 
        delta3[2] * dnorm(x, mu3[2], sigma3[2]) + 
        delta3[3] * dnorm(x, mu3[3], sigma3[3]), add = T, col = "black", lty = 2, lwd = 2, n = 200)
legend("topright", legend = c("Resting", "Low activity", "High activity"), 
       col = color, lwd = 2, bty = "n")


## Then fit the 3 state parametric model
N = 4
par = list(mu = c(-5, -4.5, -3, -1.5), # state-dep mean
           log_sigma = rep(log(0.3), 4), # state-dep sd
           eta = rep(-2, 12)) # tpm params
dat = list(x = data$logVDBA, N = 4)

## fitting model
obj4 = MakeADFun(nll, par, silent = TRUE)
opt4 = nlminb(obj4$par, obj4$fn, obj4$gr)

## extracting parameters
mod4 = report(obj4)
delta4 = mod4$delta
mu4 = mod4$mu
sigma4 = mod4$sigma

## plotting estimated state-dependent distributions
color2 = c("#F0E442", color)
hist(data$logVDBA, breaks = 50, prob = T, bor = "white", main = "", xlab = "log(VeDBA)", ylab = "Density")
for(j in 2:N) curve(delta4[j] * dnorm(x, mu4[j], sigma4[j]), add = T, col = color2[j], lwd = 2, n = 200)
curve(delta4[1] * dnorm(x, mu4[1], sigma4[1]), add = T, col = color2[1], lwd = 2, n = 200)
curve(delta4[1] * dnorm(x, mu4[1], sigma4[1]) + 
        delta4[2] * dnorm(x, mu4[2], sigma4[2]) + 
        delta4[3] * dnorm(x, mu4[3], sigma4[3])+
        delta4[4] * dnorm(x, mu4[4], sigma4[4]), add = T, col = "black", lty = 2, lwd = 2, n = 200)
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
  delta = stationary(Gamma) # stationary distribution
  ## computing mixture weights
  alpha = exp(cbind(beta, rep(0, N))) # last column fixed to zero
  alpha = alpha / rowSums(alpha) # multinomial logit link
  REPORT(alpha) # report alpha for convenience later
  ## computing state-dependent distributions
  allprobs = matrix(1, nrow(Z), N)
  ind = which(!is.na(x))
  allprobs[ind,] = Z[ind,] %*% t(alpha) # linear: sum over weighted basis functions
  -forward(delta, Gamma, allprobs) + penalty(beta, S, lambda)
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
           Z = Z, 
           N = 3, 
           S = S, 
           lambda = rep(100, 3))

## model fitting via qreml
# system.time(
#   mod <- qreml(pnll, par, dat, random = "beta", conv_crit = "relchange")
# )
system.time(
  mod <- areml(pnll, par, dat, random = "beta")
)

summary(mod)

## extracting parameters
alpha = mod$alpha # weight matrix
delta = mod$delta # stationary distribution
mod$states = viterbi(mod = mod) # state decoding

## AIC and BIC
AIC(mod, mod4, mod3)
BIC(mod, mod4, mod3)


## visualizing results
# pdf("./case_studies/figs/caracara_log.pdf", width = 8.5, height = 4)

Dens = Z_p %*% t(alpha) # state-dependent densities

par(mfrow = c(1,2), mar = c(5,4,3.5,1) + 0.1)
hist(data$logVDBA, prob = T, breaks = 50, bor = "white",
     main = "", xlab = "log(VeDBA)", ylab = "Density")
for(j in 1:N) lines(xseq, delta[j] * Dens[,j], col = color[j], lwd = 2)
lines(xseq, colSums(t(Dens)*delta), col = 1, lty = 2, lwd = 2)

legend("topright", legend = c("Resting", "Low activity", "High activity"), 
       col = color, lwd = 2, bty = "n")

plotind = 4000:8000
plot(data$logVDBA[plotind], pch = 16, cex = 0.8, col = alpha(color[mod$states[plotind]], 0.3), 
     xlab = "Time", ylab = "log(VeDBA)", bty = "n")

# dev.off()





# Full REML and marginal ML -----------------------------------------------

jnll = function(par) {
  getAll(par, dat)
  
  Gamma = tpm(eta) # transition probability matrix
  delta = stationary(Gamma) # stationary distribution
  
  ## computing mixture weights
  alpha = exp(cbind(beta, rep(0, N))) # last column fixed to zero
  alpha = alpha / rowSums(alpha) # multinomial logit link
  REPORT(alpha) # report alpha for convenience later
  
  ## computing state-dependent distributions
  allprobs = matrix(1, nrow(Z), N)
  ind = which(!is.na(x))
  allprobs[ind,] = Z[ind,] %*% t(alpha) # linear: sum over weighted basis functions
  
  lambda = exp(loglambda)
  
  -forward(delta, Gamma, allprobs, ad = T) - 
    sum(dgmrf2(beta, 0, S, lambda, log = TRUE))
}

par = list(eta = rep(-3, N*(N-1)), 
           beta = beta, 
           loglambda = rep(log(50), 3))
dat = list(x = data$logVDBA, Z = Z, N = 3, S = S)

## REML
obj1 = MakeADFun(jnll, par, random = names(par)[names(par)!="loglambda"])
system.time(
  opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
)
## marginal ML
obj2 = MakeADFun(jnll, par, random = "beta")
system.time(
  opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
)

