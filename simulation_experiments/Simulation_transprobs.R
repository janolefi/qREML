## Simulation experiments

## packages
# install.packages("RTMB")
# devtools::install_github("janoleko/LaMa") # development version
library(LaMa)
# install.packages("scales")
library(scales) # for transparent colors


# Simulation study functions ----------------------------------------------

# functions for transition probabilities (on logit scale)
sf1 = function(x) -2 + sin(3*pi*x) + exp(1.5*x)
sf2 = function(x) 2 + cos(4*pi*x) - 2*exp(x)

# function to simulate data from the inhomogeneous HMM
simHMM = function(Gamma, delta = rep(0.5,2), mu = c(1,5), sigma = c(1,3)){
  n = dim(Gamma)[3]
  # simulating Markov chain
  s = rep(NA, n)
  s[1] = sample(1:2, 1, prob = delta)
  for(t in 2:n){
    s[t] = sample(1:2, 1, prob = Gamma_sim[s[t-1],,t])
  }
  # simulating observations
  x = rnorm(n, mean = mu[s], sd = sigma[s])
  return(list(x = x, s = s))
}

# penalized likelihood function for non-parametric fit
pnll_sim = function(par){
  getAll(par, dat)
  REPORT(lambda)
  sigma = exp(logsigma) # transformation for unconstraint optimization
  delta = c(1, exp(logitdelta))
  delta = delta / sum(delta)  # initial distribution
  Gamma = tpm(cbind(beta0, betaSpline), Z) # building tpm with LaMa
  # state-dependent probabilties
  allprobs = statedens(x, dnorm, mu, sigma)
  -forward(delta, Gamma, allprobs) + penalty(betaSpline, S, lambda)
}

# (qreml() no longer returns the lambda history in all_lambda, so the trace is
#  extracted from the saved intermediate models and those are then dropped)
store_fit = function(mod, time) {
  list(beta = mod$beta,
       lambda = mod$lambda,
       lambda_trace = sapply(mod$allmods, function(m) m$lambda), # 2 x (iterations + 1)
       edf = unlist(mod$edf),
       iter = mod$iter,
       converged = mod$converged,
       time = as.numeric(time, units = "secs")) # difftime can silently switch to minutes
}


# Simulation --------------------------------------------------------------

# number of simulation runs
simiter = 250

# drawing covariate values
set.seed(1)
z = runif(1e4)

# computing tpm based on covariate
Gamma_sim = array(dim = c(2,2,1e4))
for(t in 1:1e4){
  G = diag(2)
  G[!G] = c(exp(sf2(z[t])), exp(sf1(z[t])))
  Gamma_sim[,,t] = G / rowSums(G)
}

# creating model matrices
nb = 15 # number of basis functions
modmat = make_matrices(~ s(z, bs = "ps", k = nb), data = data.frame(z = z))
Zfull = modmat$Z
S = modmat$S

# initial parameter list
par = list(mu = c(1,5), 
           logsigma = log(c(1,3)), 
           logitdelta = 0, 
           beta0 = c(-2,2),
           betaSpline = matrix(0, 2, ncol = nb-1))



### three simulation scenarios: T = 1000, 2000, 5000

## T = 1000
n = 1000
Z = Zfull[1:n,]
dat = list(Z = Z, S = S, lambda = rep(1000, 2))

modsSim1000 = list()
for(l in 1:simiter) {
  cat("Model", l, "\n")
  
  sim = simHMM(Gamma_sim[,,1:n])
  dat$x = sim$x
  
  t1 = Sys.time()
  thismod = tryCatch(qreml(pnll_sim, par, dat, random = "betaSpline",
                           silent = 2, maxiter = 100, saveall = TRUE),
                     error = function(e) "An error occurred")
  time = Sys.time() - t1
  
  if(!is.list(thismod)){
    cat("Error\n")
    next
  } else{
    modsSim1000[[l]] = store_fit(thismod, time)
  }
  gc() # free memory after each iteration
}
modsSim1000 = modsSim1000[!sapply(modsSim1000, is.null)]
saveRDS(modsSim1000[1:200], "./simulation_experiments/mods/T1000.rds")


## T = 2000
n = 2000
Z = Zfull[1:n,]
dat = list(Z = Z, S = S, lambda = rep(1000, 2))

modsSim2000 = list()
for(l in 1:simiter) {
  cat("Model", l, "\n")
  
  sim = simHMM(Gamma_sim[,,1:n])
  dat$x = sim$x
  
  t1 = Sys.time()
  thismod = tryCatch(qreml(pnll_sim, par, dat, random = "betaSpline",
                           silent = 2, maxiter = 100, saveall = TRUE),
                     error = function(e) "An error occurred")
  time = Sys.time() - t1
  
  if(!is.list(thismod)){
    cat("Error\n")
    next
  } else{
    modsSim2000[[l]] = store_fit(thismod, time)
  }
  gc() # free memory after each iteration
}
modsSim2000 = modsSim2000[!sapply(modsSim2000, is.null)]
saveRDS(modsSim2000[1:200], "./simulation_experiments/mods/T2000.rds")


## T = 5000
n = 5000
Z = Zfull[1:n,]
dat = list(Z = Z, S = S, lambda = rep(1000, 2))

modsSim5000 = list()
for(l in 1:simiter) {
  cat("Model", l, "\n")
  
  sim = simHMM(Gamma_sim[,,1:n])
  dat$x = sim$x
  
  t1 = Sys.time()
  thismod = tryCatch(qreml(pnll_sim, par, dat, random = "betaSpline",
                           silent = 2, maxiter = 100, saveall = TRUE),
                     error = function(e) "An error occurred")
  time = Sys.time() - t1
  
  if(!is.list(thismod)){
    cat("Error\n")
    next
  } else{
    modsSim5000[[l]] = store_fit(thismod, time)
  }
  gc() # free memory after each iteration
}
modsSim5000 = modsSim5000[!sapply(modsSim5000, is.null)]
saveRDS(modsSim5000[1:200], "./simulation_experiments/mods/T5000.rds")



# Results -----------------------------------------------------------------

## reading in saved models
modsSim1000 = readRDS("./simulation_experiments/mods/T1000.rds")
modsSim2000 = readRDS("./simulation_experiments/mods/T2000.rds")
modsSim5000 = readRDS("./simulation_experiments/mods/T5000.rds")

allmods = list(modsSim1000, modsSim2000, modsSim5000)
Ts = c(1000, 2000, 5000)

## true transition probabilities on a grid
zseq = seq(0, 1, length = 200)
Z_plot = predict(modmat, newdata = data.frame(z = zseq))

Gamma_true = array(dim = c(2,2,200))
for(t in 1:200){
  G = diag(2)
  G[!G] = c(exp(sf2(zseq[t])), exp(sf1(zseq[t])))
  Gamma_true[,,t] = G / rowSums(G)
}

## Note on indexing: smooth 1 (row 1 of beta, lambda[1]) belongs to gamma_21,
## smooth 2 (row 2 of beta, lambda[2]) to gamma_12

## Summary table: one row per scenario
summarise_sim = function(mods){
  time = sapply(mods, `[[`, "time")
  iter = sapply(mods, `[[`, "iter")
  conv = sapply(mods, `[[`, "converged")
  lambda = t(sapply(mods, `[[`, "lambda"))
  edf = t(sapply(mods, `[[`, "edf"))
  
  # mean absolute error of the estimated transition probabilities over the z grid
  mae = t(sapply(mods, function(m){
    G = tpm(m$beta, Z_plot)
    c(mean(abs(G[1,2,] - Gamma_true[1,2,])),
      mean(abs(G[2,1,] - Gamma_true[2,1,])))
  }))
  
  data.frame(
    fits            = length(mods),
    converged_pct   = 100 * mean(conv),
    time_mean_s     = mean(time),
    time_median_s   = median(time),
    iter_median     = median(iter),
    iter_q25        = unname(quantile(iter, 0.25)),
    iter_q75        = unname(quantile(iter, 0.75)),
    iter_le15_pct   = 100 * mean(iter <= 15),
    iter_gt25_pct   = 100 * mean(iter > 25),
    lambda12_median = median(lambda[, 2]),
    lambda21_median = median(lambda[, 1]),
    edf12_median    = median(edf[, 2]),
    edf21_median    = median(edf[, 1]),
    mae12_mean      = mean(mae[, 1]),
    mae21_mean      = mean(mae[, 2])
  )
}

results = do.call(rbind, lapply(allmods, summarise_sim))
rownames(results) = paste0("T = ", Ts)

## everything needed for the text, one scenario per row
print(t(round(results, 3)))


## sentence-ready summaries
for(m in 1:3){
  r = results[m, ]
  cat(sprintf(paste0(
    "T = %d: %d fits, %.1f%% converged; mean time %.2f s (median %.2f s); ",
    "median %g outer iterations (IQR %g-%g), %.0f%% with <= 15, %.0f%% with > 25; ",
    "median edf %.2f (gamma_12) and %.2f (gamma_21); ",
    "mean absolute error %.3f (gamma_12) and %.3f (gamma_21)\n"),
    Ts[m], r$fits, r$converged_pct, r$time_mean_s, r$time_median_s,
    r$iter_median, r$iter_q25, r$iter_q75, r$iter_le15_pct, r$iter_gt25_pct,
    r$edf12_median, r$edf21_median, r$mae12_mean, r$mae21_mean))
}


# Visualizing results -----------------------------------------------------

color <- LaMaColors(2)

## estimated functions (top row) and penalty strength traces (bottom row)
## for one transition probability
plot_sim = function(i, j, smooth, col, ylab){
  par(mfrow = c(2,3), mar = c(5,4.4,2,1)+0.1)
  
  # estimated functions
  for(m in 1:3){
    modsSim = allmods[[m]]
    plot(zseq, Gamma_true[i,j,], type = "l", lwd = 1, col = "black",
         main = paste0("T = ", Ts[m]), ylab = ylab, xlab = expression(z[t]),
         bty = "n", ylim = c(0, 1))
    for(k in seq_along(modsSim)){
      Gamma_k = tpm(modsSim[[k]]$beta, Z_plot)
      lines(zseq, Gamma_k[i,j,], col = alpha(col, 0.1), lwd = 1)
    }
    lines(zseq, Gamma_true[i,j,], lwd = 1.5, col = "black")
    if(m == 1) legend("topleft", legend = c("True function", "Estimates"),
                      col = c("black", col), lwd = c(1, 1), bty = "n")
  }
  
  # penalty strength traces
  max_iter = max(sapply(unlist(allmods, recursive = FALSE), `[[`, "iter"))
  for(m in 1:3){
    modsSim = allmods[[m]]
    plot(NA, bty = "n", xlim = c(0, max_iter), ylim = c(0, 1000),
         ylab = "Penalty strength", xlab = "Iteration")
    for(k in seq_along(modsSim)){
      trace = modsSim[[k]]$lambda_trace[smooth, ]
      lines(seq_along(trace) - 1, trace, col = alpha(col, 0.1))
    }
    abline(v = median(sapply(modsSim, `[[`, "iter")), col = "black", lwd = 1, lty = 3)
    abline(h = median(sapply(modsSim, function(x) x$lambda[smooth])),
           col = "black", lwd = 1, lty = 3)
  }
}

## Results for gamma_12
# pdf("./simulation_experiments/figs/transprobs_simulation_new.pdf", width = 7, height = 4.5)
plot_sim(i = 1, j = 2, smooth = 2, col = color[1], ylab = expression(gamma[12]^(t)))
# dev.off()

## Results for gamma_21
# pdf("./simulation_experiments/figs/transprobs_simulation2_new.pdf", width = 7, height = 4.5)
plot_sim(i = 2, j = 1, smooth = 1, col = color[2], ylab = expression(gamma[21]^(t)))
#dev.off()
