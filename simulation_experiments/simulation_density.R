## Simulation experiment: nonparametric emission densities (Supplement Section 1.2)

## packages
library(LaMa)
# devtools::install_github("janolefi/RTMBdist")
library(RTMBdist) # skew t distribution
# install.packages("scales")
library(scales) # for transparent colors

color <- LaMaColors(2)


# Simulation study functions ----------------------------------------------

## true model: 2-state HMM with skew t emissions (Azzalini type)
# state 1: right-skewed, heavier tail (df = 6) pointing towards state 2
# state 2: left-skewed, lighter tail (df = 8) pointing towards state 1
# -> sharp outer edges, heavy tails in the overlap region (overlap coefficient ~0.15)
Gamma_true <- matrix(c(0.95, 0.05,
                      0.1, 0.9), 2, byrow = TRUE)
delta_true <- stationary(Gamma_true) # (2/3, 1/3)
theta <- list(mu = c(0, 4), sigma = c(1, 1.2), skew = c(4, -4), df = c(6, 8))
dtrue <- function(x) sapply(1:2, function(j)
  dskewt(x, theta$mu[j], theta$sigma[j], theta$skew[j], theta$df[j]))

# plot true state-dependent distributions
par(mfrow = c(1,1))
curve(delta_true[1] * dskewt(x, theta$mu[1], theta$sigma[1], theta$skew[1],  theta$df[1]),
      xlim = c(-2, 6), n = 300, lwd = 2, col = color[1], bty = "n", ylab = "Density")
curve(delta_true[2] * dskewt(x, theta$mu[2], theta$sigma[2], theta$skew[2],  theta$df[2]),
      add = TRUE, n = 300, lwd = 2, col = color[2])
curve(delta_true[1] * dskewt(x, theta$mu[1], theta$sigma[1], theta$skew[1],  theta$df[1]) + 
        delta_true[2] * dskewt(x, theta$mu[2], theta$sigma[2], theta$skew[2],  theta$df[2]),
      n = 300, lwd = 1, lty = 2, add = TRUE)
legend("topright", legend = c("State 1", "State 2", "Marginal"), lwd = c(2,2,1),
       lty = c(1,1,2), col = c(color, 1), bty = "n")


# function to simulate data from the HMM
simHMM <- function(n){
  s <- x <- rep(NA, n)
  s[1] <- sample(1:2, 1, prob = delta_true) # intitial state
  x[1] <- rskewt(1, theta$mu[s[1]], theta$sigma[s[1]], theta$skew[s[1]], theta$df[s[1]])
  for(t in 2:n) { # sample next state
    s[t] <- sample(1:2, 1, prob = Gamma_true[s[t-1],])
    x[t] <- rskewt(1, theta$mu[s[t]], theta$sigma[s[t]], theta$skew[s[t]], theta$df[s[t]])
  }
  data.frame(x = x, s = s)
}

# penalized likelihood for nonparametric fit (same as caracara case study)
pnll <- function(par){
  getAll(par, dat)
  Gamma <- tpm(eta)
  alpha <- exp(cbind(beta, rep(0, 2))) # last column fixed to zero
  alpha <- alpha / rowSums(alpha) # multinomial logit link
  REPORT(alpha)
  allprobs <- Z %*% t(alpha) # linear: sum over weighted basis functions
  -forward("stationary", Gamma, allprobs) + penalty(beta, S, lambda)
}

# state-dependent densities of parametric models: Gaussian or skew t (oracle)
dpar <- function(x, p, oracle){
  sigma <- exp(p$logsigma)
  if(oracle) f <- function(j) dskewt(x, p$mu[j], sigma[j], p$skew[j], exp(p$logdf[j]))
  else f <- function(j) dnorm(x, p$mu[j], sigma[j])
  cbind(f(1), f(2))
}

# likelihood for parametric fits
nll <- function(par){
  getAll(par, dat)
  Gamma <- tpm(eta)
  allprobs <- dpar(x, par, oracle)
  -forward("stationary", Gamma, allprobs)
}

# ISE, Hellinger distance and decoding accuracy after matching states to truth
evaluate <- function(xseq, fhat, states, s){
  h <- xseq[2] - xseq[1]
  f <- dtrue(xseq)
  ise <- function(i, j) sum((fhat[,i] - f[,j])^2) * h
  perm <- if(ise(1,1) + ise(2,2) <= ise(2,1) + ise(1,2)) 1:2 else 2:1
  fhat <- fhat[, perm]
  data.frame(state = 1:2,
             ISE = colSums((fhat - f)^2) * h,
             Hellinger = sqrt(pmax(0, 1 - colSums(sqrt(fhat * f)) * h)),
             accuracy = mean(match(states, perm) == s),
             perm = perm)
}


# Simulation --------------------------------------------------------------

run_sim <- TRUE # set to TRUE to rerun (takes a while)
nsim <- 200
Ts <- c(1000, 5000)

# initial parameters for parametric fits
par_norm <- list(mu = c(1, 3.5), logsigma = log(c(1, 1)), eta = rep(-2, 2))
par_skewt <- c(par_norm, list(skew = c(1, -1), logdf = log(c(10, 10))))

if(run_sim){
  for(n in Ts){
    res <- list(); dens <- list()
    for(l in 1:nsim){
      cat("T =", n, "- run", l, "\n")
      set.seed(n + l)
      sim <- simHMM(n)

      ## qREML nonparametric fit: cubic B-splines with 2nd-order difference penalty.
      ## Number of basis functions k_l chosen for a constant knot spacing of ~0.4
      ## (spacing ~ range / (k - 3)), at least 25. A fixed k would make the resolution depend
      ## on the most extreme observations: heavy-tailed data widen the knot range and
      ## coarsen the basis near the sharp mode of state 1.
      k_l = max(25, ceiling(diff(range(sim$x)) / 0.4) + 3) # knot spacing ~0.4, at least 25 basis functions
      smoothDens = suppressMessages(smooth_dens_construct(
        data.frame(x = sim$x), par = list(x = list(mean = c(1, 3.5), sd = c(1, 1))), k = k_l))
      xseq = smoothDens$xseq$x # 500 point grid over knot range for metrics
      dat = list(x = sim$x, Z = smoothDens$Z$x, S = smoothDens$S$x, lambda = rep(20, 2))
      par = list(eta = rep(-2, 2), beta = smoothDens$coef$x)

      t1 = Sys.time()
      mod = tryCatch(qreml(pnll, par, dat, random = "beta", silent = 2),
                     error = function(e) NULL)
      time = as.numeric(Sys.time() - t1, units = "secs")
      if(is.null(mod)) next
      fhat = smoothDens$Z_predict$x %*% t(mod$alpha)
      states = viterbi(mod)
      ev = evaluate(xseq, fhat, states, sim$s)
      res[[length(res)+1]] = cbind(T = n, rep = l, method = "qREML", time = time,
                                   iter = mod$iter, converged = mod$converged, k = k_l,
                                   edf = unlist(mod$edf)[ev$perm],
                                   lambda = mod$lambda[ev$perm], ev)
      dens[[l]] = list(xseq = xseq, fhat = fhat)
      
      gc() # R + RTMB doesn't do good memory management
      
      ## parametric fits: misspecified Gaussian and correctly specified skew t
      for(oracle in c(FALSE, TRUE)){
        dat = list(x = sim$x, oracle = oracle)
        t1 = Sys.time()
        obj = MakeADFun(nll, if(oracle) par_skewt else par_norm, silent = TRUE)
        opt = tryCatch(nlminb(obj$par, obj$fn, obj$gr), error = function(e) NULL)
        time = as.numeric(Sys.time() - t1, units = "secs")
        if(is.null(opt)) next
        rep_l = lama(obj)
        fhat = dpar(xseq, obj$env$parList(opt$par), oracle)
        states = viterbi(rep_l)
        res[[length(res)+1]] = cbind(T = n, rep = l, method = ifelse(oracle, "skew t (oracle)", "Gaussian"),
                                     time = time, iter = NA, converged = opt$convergence == 0,
                                     edf = NA, lambda = NA, k = NA,
                                     evaluate(xseq, fhat, states, sim$s))
      }
      gc() # R + RTMB doesn't do good memory management
    }
    saveRDS(list(res = do.call(rbind, res), dens = dens),
            paste0("./simulation_experiments/mods/dens_T", n, ".rds"))
  }
}


# Reading in saved results ------------------------------------------------

sims = lapply(Ts, function(n) readRDS(paste0("./simulation_experiments/mods/dens_T", n, ".rds")))
res = do.call(rbind, lapply(sims, `[[`, "res"))
res$method = factor(res$method, levels = c("qREML", "Gaussian", "skew t (oracle)"))


# Summary table -----------------------------------------------------------

# number of successful fits
table(res$method[res$state == 1], res$T[res$state == 1])

msd = function(x) c(mean = mean(x), sd = sd(x))
tab = aggregate(cbind(ISE, Hellinger, accuracy, time, iter, edf) ~ method + T + state,
                data = res, FUN = msd, na.action = na.pass)
tab
# write.csv(do.call(data.frame, tab), "./revision_results/task6_density_summary.csv", row.names = FALSE)

aggregate(converged ~ method + T, data = res[res$state == 1, ], FUN = mean)

range(res$k, na.rm = TRUE)

# Visualizing results -----------------------------------------------------

xplot = seq(-3, 9, length = 500)

## estimated densities from 200 replicates vs truth (qREML)
# pdf("./simulation_experiments/figs/density_simulation.pdf", width = 8, height = 4)
par(mfrow = c(1,2), mar = c(5,4,3,1) + 0.1)
for(m in 1:2){
  plot(NA, xlim = range(xplot), ylim = c(0, 0.8), bty = "n", main = paste0("T = ", Ts[m]),
       xlab = "x", ylab = "Density")
  for(k in 1:nsim){
    d = sims[[m]]$dens[[k]]
    if(is.null(d)) next
    # align states with truth via ISE
    f = dtrue(d$xseq)
    if(sum((d$fhat - f)^2) > sum((d$fhat[,2:1] - f)^2)) d$fhat = d$fhat[,2:1]
    for(j in 1:2) lines(d$xseq, d$fhat[,j], col = alpha(color[j], 0.15))
  }
  for(j in 1:2) lines(xplot, dtrue(xplot)[,j], lwd = 1.5)
}
legend("topright", legend = c("True densities", "State 1", "State 2"),
       col = c("black", color), lwd = 1, bty = "n")
# dev.off()

## Hellinger distance by method and T
# pdf("./simulation_experiments/figs/density_simulation_hellinger.pdf", width = 8, height = 4)
par(mfrow = c(1,2), mar = c(5,4,3,1) + 0.1)
for(j in 1:2){
  boxplot(Hellinger ~ method + T, data = res[res$state == j,], main = paste("State", j),
          col = alpha(color[j], 0.5), border = color[j], frame = FALSE,
          xlab = "", ylab = "Hellinger distance", las = 2, cex.axis = 0.7)
}
# dev.off()

