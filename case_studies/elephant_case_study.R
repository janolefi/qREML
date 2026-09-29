## Ivory coast elephant case study

## packages
library(mvtnorm)
# install.packages("RTMB")
# devtools::install_github("janoleko/LaMa") # development version
library(LaMa)
# install.packages("scales") # for muted colors
library(scales)


## data
data = read.csv("./data/elephant_data.csv")
head(data)
nrow(data)


## colors
color = LaMaColors(2) # c("orange", "deepskyblue")


# Fitting homogeneous and parametric HMMs ---------------------------------

nll_hom = function(par){
  getAll(par, dat)
  ## transforming parameters
  mu = exp(logmu)
  sigma = exp(logsigma)
  kappa = exp(logkappa)
  Gamma = tpm(eta)
  delta = stationary(Gamma)
  ## reporting quantities of interest later
  REPORT(mu)
  REPORT(sigma)
  REPORT(kappa)
  ## computing state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) *
      dvm(angle[ind], 0, kappa[j])
  }
  ## forward algorithm
  -forward(delta, Gamma, allprobs)
}

N = 2
par = list(logmu = log(c(0.2, 2)),
           logsigma = log(c(0.2, 2)),
           logkappa = log(c(0.2, 1)),
           eta = rep(-2, 2))

dat = list(step = data$step, angle = data$angle, N = 2)

obj_hom = MakeADFun(nll_hom, par, silent = TRUE)
opt_hom = nlminb(obj_hom$par, obj_hom$fn, obj_hom$gr)

mod_hom = report(obj_hom)

mu = mod_hom$mu
sigma = mod_hom$sigma
kappa = mod_hom$kappa
delta = mod_hom$delta

mod_hom$states = viterbi(mod = mod_hom)

AIC(mod_hom)
BIC(mod_hom)


# pdf("./case_studies/figs/elephant_marginal.pdf", width = 8, height = 4)

par(mfrow = c(1,2))
hist(data$step, breaks = 100, prob = T, bor = "white", xlim = c(0,5), main = "", xlab = "Step length (km)", ylab = "Density")
curve(delta[1]*dgamma2(x, mu[1], sigma[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dgamma2(x, mu[2], sigma[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[1]*dgamma2(x, mu[1], sigma[1])+delta[2]*dgamma2(x, mu[2], sigma[2]),
        add = T, lwd = 2, lty = 2, n = 500)
legend("topright", legend = c("Encamped", "Exploratory", "Marginal"), col = c(color[1], color[2], "black"), 
       lty = c(1,1,2), bty = "n")

br <- seq(-pi, pi, length = 20)
hist(data$angle, breaks = br, prob = T, bor = "white", main = "", 
     xlab = "Turning angle (radians)", ylab = "Density")
curve(delta[1]*dvm(x, 0, kappa[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dvm(x, 0, kappa[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[1]*dvm(x, 0, kappa[1])+delta[2]*dvm(x, 0, kappa[2]), 
      add = T, lwd = 2, lty = 2, n = 500)

# dev.off()


## periodic variation in state process
nll_par = function(par){
  getAll(par, dat)
  ## transforming parameters
  mu = exp(logmu)
  sigma = exp(logsigma)
  kappa = exp(logkappa)
  Gamma = tpm(beta, Z)
  delta = stationary_p(Gamma, t = tod[1])
  ## reporting quantities of interest later
  REPORT(mu)
  REPORT(sigma)
  REPORT(kappa)
  ## computing state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) *
      dvm(angle[ind], 0, kappa[j])
  }
  ## forward algorithm
  -forward(delta, Gamma[,,tod], allprobs)
}

## prepwork: building trigonometric basis design matrix
Z = cosinor((1:12)*2-1, period = c(24, 12))

par = list(logmu = log(c(0.2, 2)),
           logsigma = log(c(0.2, 2)),
           logkappa = log(c(0.2, 1)),
           beta = matrix(c(rep(-2,2), rep(0, 8)), nrow = 2))

dat = list(step = data$step, angle = data$angle, N = 2, tod = data$tod)

obj_par = MakeADFun(nll_par, par, silent = TRUE)
opt_par = nlminb(obj_par$par, obj_par$fn, obj_par$gr)

mod_par = report(obj_par)
mod_par$Gamma = tpm(mod_par$beta, Z) # store only 12 unique tpms
mod_par$Delta = stationary_p(mod_par$Gamma)

AIC(mod_par)
BIC(mod_par)


## confidence intervals by sampling from MLE distribution
nSamples <- 10000
samples <- MCreport(obj_par, nSamples = nSamples)

n = 200 # plotting resolution
tod_seq = seq(0,24, length = n) # sequence of time of day values
Z_p = cosinor(tod_seq, period = c(24, 12)) # prediction design matrix

Gamma_boot_par = array(dim = c(2, 2, n, nSamples))
for(b in 1:nSamples){
  Gamma_boot_par[,,,b] = tpm(samples$beta[[b]], Z_p)
}
GammaCI_par = apply(Gamma_boot_par, c(1,2,3), quantile, probs = c(0.025, 0.975))
Gamma_plot_par = tpm(mod_par$beta, Z_p)



# Transition probabilities as smooth functions ----------------------------

## penalized likelihood function
pnll = function(par) {
  getAll(par, dat)
  
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)
  
  Gamma = tpm(cbind(beta0, betaspline), Z)
  delta = stationary_p(Gamma, t = tod[1])
  
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) * 
      dvm(angle[ind],0,kappa[j])
  }
  
  -forward(delta, Gamma[,,tod], allprobs) + 
    penalty(betaspline, S, lambda) # computes 0.5 * lambda t(b) S b
}

nb = 10
modmat = make_matrices(~ s(tod, bs = "cp", k = nb), # default number of basis functions
                       data = data.frame(tod = (1:12)*2-1),
                       knots = list(tod = c(0, 24))) # telling mgcv where to wrap the basis
Z = modmat$Z
S = modmat$S


# initial parameter values
N = 2 
par = list(logmu = log(c(0.35, 1.1)),
           logsigma = log(c(0.25, 0.75)),
           logkappa = log(c(0.2, 0.7)),
           beta0 = c(-2,-2),
           betaspline = matrix(0, nrow = N*(N-1), ncol = nb-1))

dat = list(step = data$step, angle = data$angle,
           tod = data$tod,
           N = 2, Z = Z, S = S, lambda = rep(1e5, 2))

system.time(
  mod <- qreml(pnll, par, dat, random = "betaspline", saveall = TRUE)
)

summary(mod)

## extracting parameters
beta = mod$beta

# Sampling reported quantities from the MLE distribution
nSamples <- 5000
samples <- MCreport(mod$obj, nSamples = nSamples)
betaSim <- lapply(1:nSamples, function(i){
  cbind(samples$beta0[[i]], samples$betaspline[[i]])
})


# calculating transition probabilities and stationary distribution for plotting
Z_p = predict(modmat, newdata = data.frame(tod = tod_seq)) # prediction design matrix
Gamma_p = tpm(beta, Z_p)

# calculating periodically stationary distribution
Delta_cont = matrix(NA, length(tod_seq), N)
for(t in 1:length(tod_seq)){
  t_seq = (tod_seq[t] + (1:12)*2 - 1) %% 24
  Z_cont = predict(modmat, newdata = data.frame(tod = t_seq))
  Delta_cont[t,] = stationary_p(tpm(beta, Z_cont), t = 1)
}


# computing confidence bands from Monte Carlo samples
# getting the mle out of the final model
Gamma_boot = array(dim = c(2, 2, length(tod_seq), nSamples))
Delta_boot = array(dim = c(length(tod_seq), 2, nSamples))
for(b in 1:nSamples){
  Gamma_boot[,,,b] = tpm(betaSim[[b]], Z_p)
}

# periodically stationary distribution uncertainty, on prediction grid
for(t in 1:length(tod_seq)) {
  if(t %% 10 == 0) print(paste0(round(100*t/length(tod_seq)), "%"))
  t_seq = (tod_seq[t] + (1:12)*2 - 1) %% 24
  Z_t = predict(modmat, newdata = data.frame(tod = t_seq))
  for(b in 1:nSamples){
    Delta_boot[t,,b] = stationary_p(tpm(betaSim[[b]], Z_t), t = 1)
  }
}
# pointwise confidence intervals
GammaCI = apply(Gamma_boot, c(1,2,3), quantile, probs = c(0.025, 0.975), na.rm = T)
DeltaCI = apply(Delta_boot, c(1,2), quantile, probs = c(0.025, 0.975))




## Visualising results: comparing parametric and nonparametric

# pdf("./case_studies/figs/elephant_transprobs.pdf", width = 8, height = 4)

par(mfrow = c(1,2))

# parametric fit
plot(tod_seq, Gamma_plot_par[1,2,], type = "l", lwd = 1, bty = "n", main = "Parametric",
     xlab = "Time of day", ylab ="Transition probability", xaxt = "n", ylim = c(0,1))
polygon(c(tod_seq, rev(tod_seq)), c(GammaCI_par[1,1,2,], rev(GammaCI_par[2,1,2,])), 
        col = alpha("black", 0.1), border = F)
lines(tod_seq, Gamma_plot_par[2,1,], lwd = 1, lty = 5)
polygon(c(tod_seq, rev(tod_seq)), c(GammaCI_par[1,2,1,], rev(GammaCI_par[2,2,1,])), 
        col = alpha("black", 0.1), border = F)
legend(x = -0.5, y = 1.02, lty = c(1,5), y.intersp = 1.3,
       legend = c(expression(gamma[12]^(t)), expression(gamma[21]^(t))), bty = "n")
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

# non-parametric fit
plot(tod_seq, Gamma_p[1,2,], type = "l", lwd = 1, bty = "n", main = "Nonparametric",
     xlab = "Time of day", ylab ="Transition probability", xaxt = "n", ylim = c(0,1))
polygon(c(tod_seq, rev(tod_seq)), c(GammaCI[1,1,2,], rev(GammaCI[2,1,2,])), 
        col = alpha("black", 0.1), border = F)
lines(tod_seq, Gamma_p[2,1,], lwd = 1, lty = 5)
polygon(c(tod_seq, rev(tod_seq)), c(GammaCI[1,2,1,], rev(GammaCI[2,2,1,])), 
        col = alpha("black", 0.1), border = F)
# legend(x = -0.5, y = 1.02, lty = c(1,2), y.intersp = 1.3,
#        legend = c(expression(gamma[12]^(t)), expression(gamma[21]^(t))), bty = "n")
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

# dev.off()


## nonparametric periodically stationary distribution

sun_cycle_colors <- c(
  "#0b0d3e",  # 00:00 - Midnight (Night)
  "#0d1046",  # 00:30
  "#0f124e",  # 01:00
  "#121557",  # 01:30
  "#14185f",  # 02:00
  "#161a67",  # 02:30
  "#191d6f",  # 03:00
  "#1b1f78",  # 03:30
  "#1d227f",  # 04:00
  "#1f2587",  # 04:30
  "#22288f",  # 05:00
  "#564d7b",  # 05:30
  "#877060",  # 06:00
  "#e67e22",  # 06:30 - Sunrise (Vibrant Orange)
  "#f28c32",  # 07:00
  "#f89b42",  # 07:30
  "#fba554",  # 08:00 - Morning
  "#fcaf67",  # 08:30
  "#fcc978",  # 09:00
  "#fcd68d",  # 09:30
  "#fce3a2",  # 10:00
  "#fde1b7",  # 10:30
  "#feeacb",  # 11:00
  "#d8eef5",  # 11:30 - Gradual shift from morning to day
  "#b7e9f9",  # 12:00 - Noon (Light Blue)
  "#a2e2f9",  # 12:30
  "#8cdcf9",  # 13:00
  "#75d5fa",  # 13:30
  "#5ecefa",  # 14:00
  "#47c7fa",  # 14:30
  "#30c0fb",  # 15:00
  "#2ca9e2",  # 15:30
  "#2593cb",  # 16:00 - Afternoon (Deeper Blue)
  "#1f7db4",  # 16:30
  "#19679d",  # 17:00
  "#7f6d63",  # 17:30
  "#e57328",  # 18:00 - Sunset (Vibrant Orange)
  "#cc6925",  # 18:30
  "#b36022",  # 19:00
  "#9a571f",  # 19:30
  "#7c4550",  # 20:00 - Nightfall
  "#603a59",  # 20:30
  "#4c2f55",  # 21:00
  "#38264b",  # 21:30
  "#291f3f",  # 22:00 - Nighttime
  "#1c1732",  # 22:30
  "#140f28",  # 23:00 - Late Night
  "#0e0a1e"   # 23:30
)



# pdf("./case_studies/figs/elephant_stationary.pdf", width = 7, height = 4.5)

par(mfrow = c(1,1), mar = c(5,4,0,2)+0.1)
plot(NA, bty = "n", ylim = c(0,1), xlim = c(0,24), ylab = "Pr(exploratory)", 
     xlab = "Time of day", xaxt = "n", cex = 1.5)
polygon(x = c(0, 6.5, 6.5, 0), y = c(0, 0, 1, 1), col = "gray95", border = "white")
polygon(x = c(18.5, 24, 24, 18.5), y = c(0, 0, 1, 1), col = "gray95", border = "white")
# polygon(x = c(0, 24, 24, 0), y = c(-0.05, -0.05, -0.01, -0.01), col = "black", border = "black")
for(t in 0:47){
  polygon(x = c(t/2, (t+1)/2, (t+1)/2, t/2), y = c(-0.04, -0.04, -0.01, -0.01), col = sun_cycle_colors[t+1], border = sun_cycle_colors[t+1])
}
# for(t in 0:95){
#   polygon(x = c(t/4, (t+1)/4, (t+1)/4, t/4), y = c(-0.04, -0.04, -0.01, -0.01), col = sun_cycle_colors[t+1], border = sun_cycle_colors[t+1])
# }
lines(tod_seq, Delta_cont[,2], lwd = 1.5, col = "black",
     bty = "n", ylim = c(0,1), xaxt = "n")
polygon(c(tod_seq, rev(tod_seq)), c(DeltaCI[1,,2], rev(DeltaCI[2,,2])), 
        col = alpha("black", 0.1), border = F)
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

# dev.off()



# Plotting the model sequence ---------------------------------------------

## plotting the model sequence
mods = mod$allmods
length(mods)
plotind = c(1,2,3,5,7,length(mods))

# pdf("./case_studies/figs/elephant_modseq.pdf", width = 8.5, height = 5.5)

par(mfrow = c(2,3), mar = c(5,4,3,1))
for(m in plotind){
  beta = mods[[m]]$beta
  
  Gamma_p = tpm(beta, Z_p)
  
  plot(tod_seq, Gamma_p[1,2,], type = "l", lwd = 1, bty = "n", main = paste("Iteration", m),
       xlab = "Time of day", ylab ="Transition probability", xaxt = "n", ylim = c(0,1))
  lines(tod_seq, Gamma_p[2,1,], lwd = 1, lty = 5)
  legend(x = -0.5, y = 1.02, lty = c(1,2), y.intersp = 1.3,
         legend = c(expression(gamma[12]^(t)), expression(gamma[21]^(t))), bty = "n")
  axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

}

# dev.off()




# Full REML and marginal ML -----------------------------------------------

## joint likelihood (has complete normal density for random effects)
jnll = function(par) {
  getAll(par, dat)
  
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)
  
  Gamma = tpm(cbind(beta0, betaspline), Z)
  delta = stationary_p(Gamma, t = tod[1])
  
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) * 
      dvm(angle[ind],0,kappa[j])
  }
  lambda = exp(loglambda)
  
  -forward(delta, Gamma[,,tod], allprobs) - 
    dgmrf(betaspline[1,], 0, lambda[1]*S, log = TRUE) -
    dgmrf(betaspline[2,], 0, lambda[2]*S, log = TRUE)
}

## initial parameters
par = list(logmu = log(c(0.35, 1.1)),
           logsigma = log(c(0.25, 0.75)),
           logkappa = log(c(0.2, 0.7)),
           beta0 = c(-2,-2),
           betaspline = matrix(0, nrow = N*(N-1), ncol = ncol(Z)-1),
           loglambda = rep(log(1e5),2))

dat = list(step = data$step, angle = data$angle, 
           tod = data$tod, 
           N = 2, Z = Z, S = as(S[[1]], "sparseMatrix"))

## creating objective function
## full REML
t1 <- Sys.time()
obj1 <- MakeADFun(jnll, par, random = names(par)[names(par)!="loglambda"])
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
Sys.time() - t1

## marginal ML
t1 <- Sys.time()
obj2 <- MakeADFun(jnll, par, random = "betaspline")
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr, control = list(iter.max = 500))
Sys.time() - t1

mod2 = report(obj1)
beta2 = mod2$beta

Gamma_plot2 = tpm(beta2, Z_p)

par(mfrow = c(1,1))
plot(tod_seq, Gamma_plot2[1,2,], type = "l", lwd = 2, bty = "n",
     xlab = "time of day", ylab = "transition probability", xaxt = "n")
lines(tod_seq, Gamma_plot2[2,1,], lwd = 2, lty = 3)
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))





# Comparison to hmmTMB ----------------------------------------------------

library(hmmTMB)

data$hour = data$tod*2
data$step[which(data$step == 0)] <- NA

hid <- MarkovChain$new(
  data = data,
  formula = ~ s(hour, bs = "cp", k = nb),
  n_states = 2,
  gam_args = list(knots = list(hour = c(0, 24)))
)

obs <- Observation$new(
  data = data,
  dists = list(step = "gamma2", angle = "vm"),
  par = list(step = list(mean = c(0.35, 1.1), sd = c(0.25, 0.75)), 
             angle = list(mu = c(0,0), kappa = c(0.2, 0.7))),
  n_states = 2,
  fixpar = list(obs = list(angle = list(mu = rep(NA, 2))))
)

hmm <- HMM$new(obs = obs, hid = hid)

system.time(hmm$fit())

hmm$plot(what = "tpm", i = )
