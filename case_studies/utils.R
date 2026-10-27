## Effective degrees of freedom at the estimate of a fitted RTMB object
# obj:    fitted RTMB object (e.g. full REML) with log smoothing parameters
# nll:    the joint negative log-likelihood that obj was built from
# random: names of the spline coefficient parameters (matrix: one smooth per row)
# S:      penalty matrix, or list of penalty matrices (one per element of random)
# lsp:    name of the log smoothing parameter vector
edf_at <- function(obj, nll, random, S, lsp = "log_lambda"){
  if(!is.list(S)) S <- list(S)
  if(length(S) == 1) S <- rep(S, length(random)) # one penalty matrix for all
  
  # all parameters (fixed and random) at the estimate
  pars <- obj$env$parList(par = obj$env$last.par.best)
  lam <- exp(pars[[lsp]])
  
  # same likelihood with lambda fixed: its Hessian is the penalized Hessian J
  map <- setNames(list(factor(rep(NA, length(lam)))), lsp)
  obj_fix <- RTMB::MakeADFun(nll, pars, map = map, silent = TRUE)
  J <- obj_fix$he(obj_fix$par)
  Jinv <- tryCatch(solve(J), error = function(e){
    solve(J + diag(1e-8 * max(abs(diag(J))), nrow(J))) # ridge for (near-)flat directions
  })
  
  # edf_i = K_i - lambda_i * tr([J^-1]_ii S_i), smooth by smooth
  edf <- c()
  l <- 1 # running index into lambda
  for(r in seq_along(random)){
    b <- pars[[random[r]]]
    n_smooth <- if(is.matrix(b)) nrow(b) else 1
    idx <- matrix(which(names(obj_fix$par) == random[r]), nrow = n_smooth)
    for(j in 1:n_smooth){
      ii <- idx[j, ]
      edf <- c(edf, length(ii) - lam[l] * sum(diag(Jinv[ii, ii] %*% S[[r]])))
      l <- l + 1
    }
  }
  edf
}



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