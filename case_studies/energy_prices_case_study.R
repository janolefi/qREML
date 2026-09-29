## Spanish energy prices case study

# install.packages("RTMB")
# install.packages("MSwM")
library(MSwM) # data
# devtools::install_github("janoleko/LaMa") # development version
library(LaMa) # HMM functions and qreml
# install.packages("scales")
library(scales) # just for transparency in plotting


## data
data(energy, package = "MSwM")
nrow(energy)
head(energy)

## colors
color = LaMaColors(2) # c("orange", "deepskyblue")


# Fitting MS-GAMLSS -------------------------------------------------------

## penalized likelihood
pnll = function(par) {
  getAll(par, dat)

  Gamma = tpm(eta) # computing the tpm
  delta = stationary(Gamma) # stationary distribution

  beta = cbind(beta0, betaspline); REPORT(beta) # mean parameter matrix
  alpha = cbind(alpha0, alphaspline);  REPORT(alpha) # sd parameter matrix
  Mu = Z %*% t(beta) # mean
  Sigma = exp(Z %*% t(alpha)) # sd
  allprobs = cbind(dnorm(price, Mu[,1], Sigma[,1]),
                   dnorm(price, Mu[,2], Sigma[,2])) # state-dependent densities
  
  - forward(delta, Gamma, allprobs) +
    penalty(list(betaspline, alphaspline), S, lambda) # penalty does the heavy lifting and reports to qreml
}

## prepwork: model matrices (uses mgcv under the hood)
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

# map = list(eta = rep(1, 2))

## model fit
system.time(
  mod <- qreml(pnll, par, dat, random = c("betaspline", "alphaspline"), 
               saveall = TRUE)
)

# Estimated penalty strengths and edf
summary(mod)


## extracting parameters
beta = mod$beta # mean parameter matrix
alpha = mod$alpha # sd parameter matrix
Gamma = mod$Gamma # t.p.m.
round(Gamma, 3)
round(c(Gamma[1,2], Gamma[2,1])^-1) # mean dwell time
(delta = mod$delta) # stationary distribution
mod$states = viterbi(mod = mod) # decoding most probable state sequence

Mu_p = Z_p %*% t(beta)
Sigma_p = exp(Z_p %*% t(alpha))


## visualizing results
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
length(mods)
plotind = c(1,2,3,4,6,length(mods))

# pdf("./case_studies/figs/energy_oil_modseq.pdf", width = 8, height = 5.5)

par(mfrow = c(2,3), mar = c(5,4,3,1))
for(m in plotind){
  beta = mods[[m]]$beta
  alpha = mods[[m]]$alpha
  
  states = viterbi(mod = mods[[m]])
  Mu_p = Z_p %*% t(beta)
  Sigma_p = exp(Z_p %*% t(alpha))
  
  plot(energy$Oil, energy$Price, pch = 20, bty = "n", col = alpha(color[states], 0.1),
       xlab = "Oil price", ylab = "Energy price", main = paste("Iteration", m))
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

## joint likelihood (has complete normal density for random effects)
jnll = function(par) {
  getAll(par, dat)
  
  Gamma = tpm(eta) # computing the tpm
  delta = stationary(Gamma) # stationary distribution
  
  beta = cbind(beta0, betaspline); REPORT(beta) # mean parameter matrix
  alpha = cbind(alpha0, alphaspline); REPORT(alpha) # sd parameter matrix
  
  Mu = Z %*% t(beta) # mean
  Sigma = exp(Z %*% t(alpha)) # sd
  
  allprobs = cbind(dnorm(price, Mu[,1], Sigma[,1]), 
                   dnorm(price, Mu[,2], Sigma[,2])) # state-dependent densities
  
  lambda = exp(loglambda)
  
  - forward(delta, Gamma, allprobs) - 
    sum(dgmrf2(rbind(betaspline, alphaspline), 0, S, lambda, log = TRUE))
}


## initial parameter list
par = list(eta = rep(-3, 2),
           beta0 = c(2, 5),
           betaspline = matrix(0, nrow = 2, ncol = nb-1),
           alpha0 = c(0, 0),
           alphaspline = matrix(0, nrow = 2, ncol = nb-1),
           loglambda = rep(log(1e5), 4))

## data, model matrices and initial penalty strength
dat = list(price = energy$Price, Z = Z, S = S[[1]])

## model fitting
## REML
t1 <- Sys.time()
obj2 <- MakeADFun(jnll, par, random = names(par)[names(par)!="loglambda"])
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr) 
Sys.time()-t1

## marginal ML
t1 <- Sys.time()
obj3 <- MakeADFun(jnll, par, random = c("betaspline", "alphaspline"))
opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr) 
Sys.time()-t1


sdr = sdreport(obj2)
round(exp(as.list(sdr, "Estimate")$loglambda), 2) # estimated lambdas 

mod2 = report(obj2)

beta = mod2$beta
alpha = mod2$alpha
mod2$states = viterbi(mod = mod2) # decoding most probable state sequence

Mu_p = Z_p %*% t(beta)
Sigma_p = exp(Z_p %*% t(alpha))


## visualizing results
# pdf("./case_studies/figs/energy_oil.pdf", width = 8, height = 4.5)

par(mfrow = c(1,2), mar = c(5,4,3,1))
plot(energy$Oil, energy$Price, pch = 20, bty = "n", col = alpha(color[mod2$states], 0.1),
     xlab = "oil price", ylab = "energy price")
for(j in 1:2) lines(xseq, Mu_p[,j], col = color[j], lwd = 3)

qseq = qnorm(seq(0.5, 0.95, length = 4))
for(i in qseq){
  for(j in 1:2){
    lines(xseq, Mu_p[,j] + i*Sigma_p[,j], col = alpha(color[j], 0.7), lwd = 1, lty = 2)
    lines(xseq, Mu_p[,j] - i*Sigma_p[,j], col = alpha(color[j], 0.7), lwd = 1, lty = 2)
  }
}
legend("topright", bty = "n", legend = paste("state", 1:2), col = color, lwd = 3)

plot(NA, xlim = c(0, nrow(energy)), ylim = c(1,10), bty = "n",
     xlab = "time", ylab = "energy price")
segments(x0 = 1:(nrow(energy)-1), x1 = 2:nrow(energy),
         y0 = energy$Price[-nrow(energy)], y1 = energy$Price[-1], col = color[mod2$states[-1]], lwd = 0.5)




# Comparison to hmmTMB ----------------------------------------------------

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

system.time(hmm$fit())

# decode the states
energy$viterbi <- factor(paste0("State ", hmm$viterbi()))

# plot mean relationship
hmm$plot(what = "obspar", var = "Oil", i = "Price.mean") +
  geom_point(aes(x = Oil, y = Price, fill = viterbi, col = viterbi),
               data = energy, alpha = 0.3) +
  theme(legend.position = "none")

# plot standard deviation
hmm$plot(what = "obspar", var = "Oil", i = "Price.sd") +
  theme(legend.position.inside = c(0.3, 0.7))
