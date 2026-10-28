## Caracaras case study: comparison with the Bayesian spline-based HMM of
## Chen & Finkenstädt (2024, JASA), "Bayesian spline-based hidden Markov models with
## applications to actimetry data and sleep analysis"
##
## The middle part is the authors' own code for their shark application (adsp_sharkN3.R
## from the supplementary material of the paper), applied to the caracara data.
## All changes to their code are marked with "## CHANGED"; the sampler itself is untouched.
## As in their shark application: N = 3, 120,000 MCMC iterations, of which 70,000 burn-in,
## thinning by 10, 21 initial interior knots (25 basis functions).
## One chain (their seed) takes about 2 hours on an Apple M2.

## packages
# install.packages(c("LaMa", "scales", "truncnorm", "fda", "extraDistr", "splines2", "DTMCPack"))
library(LaMa)
library(scales) # for transparent colors

## color palette
color = LaMaColors(3)

## data
data = read.csv("./data/caracara_downsampled.csv")


# qREML fit (as in caracaras_case_study.R) --------------------------------

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
t0 <- Sys.time()
mod <- qreml(pnll, par, dat, random = "beta")
(time_qreml <- Sys.time() - t0)

Dens = Z_p %*% t(mod$alpha) # state-dependent densities
states_qreml = viterbi(mod) # global decoding


# Chen & Finkenstädt: their code (adsp_sharkN3.R) -------------------------

## sharkN3 computed on bsu laptop
library(truncnorm)
library(fda)
library(extraDistr)
library(splines2)
library(splines)

library(DTMCPack) ## CHANGED: statdistr() is used below but DTMCPack is not loaded in the original
x <- read.csv("./data/caracara_downsampled.csv")$logVDBA ## CHANGED: caracara data instead of load("~/sharkdata.Rdata")
set.seed(1000)
n<-length(x)


Llk<-function(parvect,x,N){ # evalute -(log-likelihood) for a spline-based HMM
  as.t<-exp(parvect[1:(N*nb-N)])
  as<-vector("list")
  l=1
  for (k in 1:N){as[[k]]<-c(as.t[((l-1)*(nb-1)/2+1):(l*(nb-1)/2)], 1, as.t[(l*(nb-1)/2+1):((l+1)*(nb-1)/2)]); as[[k]]<-as[[k]]/sum(as[[k]]); l=l+2}
  delta.t<-exp(parvect[(N*nb-N+1):(N*nb-1)]) # (delta2/del1,del3/del1)
  delta<-c(1,delta.t) 
  delta<-delta/sum(delta)
  egamma<-diag(N)
  egamma[!egamma]<-exp(parvect[(N*nb):length(parvect)])
  egamma<-t(egamma)
  egamma<-egamma/apply(egamma,1,sum)
  allprobs<-matrix(NA,length(x),N)
  for (j in 1:N){allprobs[,j]<-B%*%as[[j]]}
  omega1 <- sum(delta*allprobs[1,])
  foo <- (delta*allprobs[1,])/omega1
  lscale <- log(omega1)
  for (i in 2:n){
    foo <- foo%*%egamma*allprobs[i,]
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
  }
  return(-(lscale))
}

mle2<-function(x,as0,delta0, gams0,N){ # initialize parameters based on MLE
  gams0<-t(gams0)
  parvect <- c(log(as0),log(delta0[-1]/delta0[1]),log(gams0[!diag(N)]/rep(diag(gams0),each=N-1))) ### changed
  mod     <- nlm(Llk,parvect,x,N,iterlim=1000)
  Llk    <- mod$minimum
  as.t    <- exp(mod$estimate[1:(N*nb-N)])
  as      <- vector("list")
  l=1
  for (k in 1:N){as[[k]]<-c(as.t[((l-1)*(nb-1)/2+1):(l*(nb-1)/2)], 1, as.t[(l*(nb-1)/2+1):((l+1)*(nb-1)/2)]); as[[k]]<-as[[k]]/sum(as[[k]]); l=l+2}
  delta.t<-exp(mod$estimate[(N*nb-N+1):(N*nb-1)])
  delta<-c(1,delta.t)
  delta<-delta/sum(delta)
  egamma<-diag(N)
  egamma[!egamma]<-exp(mod$estimate[(N*nb):length(mod$estimate)])
  egamma<-t(egamma)
  egamma<-egamma/apply(egamma,1,sum)
  list(Llk=Llk,as=as,gamma=egamma,delta=delta)
} 

P<- function (t) { # prior on K
  ddunif(x=t, min=3, max=50, log = FALSE)
}

b_k<- function (t) { # probability for choosing a birth move
  if (t>3) {
    0.5
  } else {
    1
  }
}

H_b <- function (z, inknot, knot, a, b) { # used in the evaluation of acceptance probability
  h_b<-rep(NA,length(inknot))
  for (i in 1:length(inknot)) {
    h_b[i]<-dtruncnorm(z,a=a, b=b, mean=inknot[i], sd=(knot[i+2]-knot[i])^alpha) # knot here is (a,inknots,b)
  }
  sum(h_b)/length(inknot)
}

cLLk<-function(Z,allprobs,gamma){ # compute complete log-likelihood of a spline-based HMM
  cllk<- log(delta[Z[1]])+sum(log(t(gamma)[N*(Z[1:(n-1)]-1)+Z[-1]]))+sum(log(t(allprobs)[N*(0:(n-1))+Z]))
  return(cllk)
}

dlogga<-function(t, shape, rate) { # prior for the reparameterized spline coefficients - log-gamma distribution 
  (rate^shape)*exp(shape*t-rate*exp(t))/gamma(shape)
}

# algorithm set-ups
N<-3 # number of states
Tb<-70000 # total number of iterations of MCMC 
T<-120000
AR1<-rep(NA,T) # record total number of acceptance for a MCMC step
AR2<-rep(NA,T)
AR4<-rep(NA,T)
accept_b<-0
accept_d<-0
M<-10 # for adaptive MCMC
A_t<-vector("list")
Inknots<-vector("list") # store sampled interior knot points
A<-vector("list") # store sampled spline coefficients
Gamma<-vector("list") # store sampled transition probabilities
Z<-vector("list")
K<-rep(NA,T) # store the number of interior knot points
logllk<-rep(NA,T) # store log-likelihood of the HMM generated during MCMC simulation
LG<-rep(NA,T) # store log joint posterior density generated during MCMC simulation
Omg<-rep(NA,T) # store the simulated shape parameter of the log-gamma prior
delta<-rep(1/N,N)
# initialization of the algorithm (specify constants, tuning parameters and obtian initial parameter values)
ord=4 # order of the B-splines
degree=ord-1 # degree of the B-splines
lambda<-21; nb<-lambda+degree+1 # initial length of inknots (interior knots); corresponds to 11+3+1=15 basis elements
tau_r<-rep(NA,T) # tuning parameters in M-H steps
tau_2<-rep(NA,T)
tau_o<-rep(NA,T)
alpha<-2 # tuning parameters in the birth move
a<- -8.25 # left boundary of knot sequence ## CHANGED: shark data: -5.5
b<- 1 # right boundary of knot sequence ## CHANGED: shark data: -1
L<- b-a  
a0<-a
b0<-b
gams0<-matrix(c(0.8,0.1,0.1,0.1,0.8,0.1,0.1,0.1,0.8),nrow=3)
as0.1<-c(rep(5,(nb-1)/2),rep(0.2,(nb-1)/2))        
as0.2<-c(rep(0.6,(nb-1)/2),rep(0.5,(nb-1)/2))
as0.3<-c(rep(0.1,(nb-1)/2),rep(2,(nb-1)/2))
as0<-c(as0.1,as0.2,as0.3) 
delta0<-rep(1/N,N)
Inknots[[1]] <- as.vector(quantile(x,probs=seq(0,1,length=lambda)))
knots<-c(rep(a,4),Inknots[[1]],rep(b,4)) 
ws<-(ibs(x=b, knots = Inknots[[1]], degree = 3, intercept = T, Boundary.knots = c(a,b))[1:(lambda+4)])^(-1) #
B <- t(t(spline.des(knots,x,degree+1,outer.ok=T)$design)*ws)

init.start<-Sys.time() ## CHANGED: also time the initialisation
mod_ini<-mle2(x,as0,delta0,gams0,N) # obtain MLEs for the model parameters based on initial knot configuration
init.time<-round(difftime(Sys.time(), init.start, units = "secs"),2) ## CHANGED

K[1]<-lambda # number of interior knot points
logllk[1]<--mod_ini$Llk
A_t[[1]] <-rep(NA,N*nb) # stacked vector of spline coefficients (reparameterized)
for (k in 1:N) {
  A_t[[1]][(1+(k-1)*nb):(k*nb)]<-log(rtruncnorm(nb,a=0,b=1,mean=mod_ini$as[[k]],sd=1e-02))
}
A[[1]]<-rep(NA,N*nb)
for (k in 1:N) {
  A[[1]][(1+(k-1)*nb):(k*nb)]<-exp(A_t[[1]][(1+(k-1)*nb):(k*nb)])/sum(exp(A_t[[1]][(1+(k-1)*nb):(k*nb)]))
}
Gamma[[1]]<-mod_ini$gamma
Omg[1]<-0.8 
tau_r[1:(M+1)]<-1
tau_2[1:(M+1)]<-1
tau_o[1:(M+1)]<-1

start.time<-Sys.time()
for (i in 2:Tb) {
  inknots<-Inknots[[i-1]]
  nb<-K[i-1]+4
  knots<-c(rep(a,4),inknots,rep(b,4))
  B <- mSpline(x, knots = inknots, degree = 3, intercept = TRUE,Boundary.knots = c(a,b))
  as0<-A[[i-1]]
  omg<-Omg[i-1]
  allprobs<-B%*%matrix(as0,nrow = nb, ncol=N)
  # FFBS 
  FF<-matrix(NA,N,n)
  FF[,1]<-delta*allprobs[1,]/sum(delta*allprobs[1,])
  for (k in 2:n) {
    for (j in 1:N) {
      FF[j,k]<-allprobs[k,j]*sum(FF[,k-1]*Gamma[[i-1]][,j])
    }
    FF[,k]<-FF[,k]/sum(FF[,k])
  }
  Z[[i]]<-rep(NA,n)
  Z[[i]][n]<-sample(1:N,size = 1, replace = T, prob=FF[,n]/sum(FF[,n]))
  for (k in (n-1):1) {Z[[i]][k]<-sample(1:N,size = 1, replace = T, prob=FF[,k]*Gamma[[i-1]][,Z[[i]][k+1]]/sum(FF[,k]*Gamma[[i-1]][,Z[[i]][k+1]]))}
  Z_counts<-matrix(table(factor(paste0(head(Z[[i]],-1),tail(Z[[i]],-1)),levels = c('11','12','13','21','22','23','31','32','33'))),N,N,byrow=T)
  # update transition matrix
  Gamma[[i]]<-matrix(NA,N,N)
  for (k in 1:N) {Gamma[[i]][k,]<-rdirichlet(1, alpha=1+Z_counts[k,])}
  # update r_k
  Llk0<-cLLk(Z=Z[[i]],allprobs,gamma=Gamma[[i]])
  kreloc_ind<-sample(x=1:length(inknots), size=1, replace = TRUE)
  krenew<-rtruncnorm(1, a=a0, b=b0, mean=inknots[kreloc_ind], sd=tau_r[i-1])
  knots<-sort(c(rep(a,4),inknots[-kreloc_ind],krenew,rep(b,4)))
  B <- mSpline(x, knots = knots[5:(K[i-1]+4)], degree = 3, intercept = TRUE,Boundary.knots = c(a,b))
  allprobs<-B%*%matrix(as0,nrow = nb, ncol=N)
  Llk1<-cLLk(Z=Z[[i]],allprobs,gamma=Gamma[[i]])
  AR1[i]<-min(1,exp(Llk1-Llk0)*(dtruncnorm(x=inknots[kreloc_ind], a=a0, b=b0, mean=krenew, sd=tau_r[i-1])/dtruncnorm(x=krenew, a=a0, b=b0, mean=inknots[kreloc_ind], sd=tau_r[i-1])))
  # acceptance-rejection
  if (runif(1)<AR1[i]) {
    Inknots[[i]]<-sort(c(inknots[-kreloc_ind],krenew))
    Llk0<-Llk1
  } else {
    Inknots[[i]]<-Inknots[[i-1]]
  }
  if (i>(M+1)) {tau_r[i]<-ifelse(mean(AR1[(i-M):i])>0.4,tau_r[i-1]+min(0.01,1/sqrt(i)),max(tau_r[i-1]-min(0.01,1/sqrt(i)),0.000001))}
  # update A_t (using RW-MH)
  knots<-c(rep(a,4),Inknots[[i]],rep(b,4))
  B <- mSpline(x, knots = knots[5:(K[i-1]+4)], degree = 3, intercept = TRUE,Boundary.knots = c(a,b))
  a_t<-A_t[[i-1]]+rnorm(N*nb,0,tau_2[i-1])
  as1<-rep(NA,N*nb)
  for (k in 1:N) {
    as1[(1+(k-1)*nb):(k*nb)]<-exp(a_t[(1+(k-1)*nb):(k*nb)])/sum(exp(a_t[(1+(k-1)*nb):(k*nb)]))
  }
  allprobs<-B%*%matrix(as1,nrow = nb, ncol=N)
  Llk1<-cLLk(Z=Z[[i]],allprobs,gamma=Gamma[[i]])
  AR2[i]<-min(1,exp(Llk1-Llk0)*(prod(dlogga(a_t, shape=omg, rate=1)/dlogga(A_t[[i-1]], shape=omg, rate=1))))
  # acceptance-rejection
  if (runif(1)<AR2[i]) {
    A_t[[i]]<-a_t
    Llk0<-Llk1
    as0<-as1
  } else {
    A_t[[i]]<-A_t[[i-1]]
  }
  if (i>(M+1)) {tau_2[i]<-ifelse(mean(AR2[(i-M):i])>0.24,tau_2[i-1]+min(0.01,1/sqrt(i)),max(tau_2[i-1]-min(0.01,1/sqrt(i)),0.000001))}
  # update omg (hyperprior ga(1,1))
  omg1<-exp(log(omg)+rnorm(1,0,tau_o[i-1]))
  AR4[i]<-min(1,(dgamma(omg1, shape=1, rate=1)/dgamma(omg, shape=1, rate=1))*(prod(dlogga(A_t[[i]], shape=omg1, rate=1)/dlogga(A_t[[i]], shape=omg, rate=1)))*(omg1/omg))
  # acceptance-rejection
  if (runif(1)<AR4[i]) {
    Omg[i]<-omg1
    omg<-omg1
  } else {
    Omg[i]<-Omg[i-1]
  }
  if (i>(M+1)) {tau_o[i]<-ifelse(mean(AR4[(i-M):i])>0.4,tau_o[i-1]+min(0.01,1/sqrt(i)),tau_o[i-1]-min(0.01,1/sqrt(i)))}
  A[[i]]<-as0
  ## RJMCMC step
  u<-runif(1)
  # birth step 
  if (u<b_k(K[i-1])) {
    kadd_ind<-sample(x=1:length(Inknots[[i]]), size=1, replace = TRUE)
    knew<-rtruncnorm(1, a=a0, b=b0, mean=Inknots[[i]][kadd_ind], sd=(knots[kadd_ind+5]-knots[kadd_ind+3])^alpha)      
    kn<-which(sort(c(knots,knew))==knew)-5 # 
    r<-rep(NA,2) # shared across states
    r[1]<-(knew-knots[kn+2])/(knots[kn+5]-knots[kn+2])
    r[2]<-(knew-knots[kn+3])/(knots[kn+6]-knots[kn+3])
    knots<-sort(c(knots,knew))
    B <- mSpline(x, knots = knots[5:(K[i-1]+5)], degree = 3, intercept = TRUE,Boundary.knots = c(a,b))
    a_t<-rep(NA,N*(nb+1))
    as0<-rep(NA,N*(nb+1))
    for (k in 1:N) {
      u<-runif(1)
      a_t[(1+(k-1)*(nb+1)):(k*(nb+1))][kn+2]<-r[1]*(A_t[[i]][(1+(k-1)*nb):(k*nb)][kn+2])+(1-r[1])*(A_t[[i]][(1+(k-1)*nb):(k*nb)][kn+1])
      a_t[(1+(k-1)*(nb+1)):(k*(nb+1))][kn+3]<-r[2]*(A_t[[i]][(1+(k-1)*nb):(k*nb)][kn+3])+(1-r[2])*(A_t[[i]][(1+(k-1)*nb):(k*nb)][kn+2])
      a_t[(1+(k-1)*(nb+1)):(k*(nb+1))][kn+4]<-u*(A_t[[i]][(1+(k-1)*nb):(k*nb)][kn+4])+(1-u)*(A_t[[i]][(1+(k-1)*nb):(k*nb)][kn+3])
      a_t[(1+(k-1)*(nb+1)):(k*(nb+1))][-c(kn+2,kn+3,kn+4)]<-A_t[[i]][(1+(k-1)*nb):(k*nb)][-c(kn+2,kn+3)]
      a_tt<-a_t[(1+(k-1)*(nb+1)):(k*(nb+1))]-max(a_t[(1+(k-1)*(nb+1)):(k*(nb+1))])
      as0[(1+(k-1)*(nb+1)):(k*(nb+1))]<-exp(a_tt)/sum(exp(a_tt))
    }
    allprobs<-B%*%matrix(as0,nrow = nb+1, ncol=N)
    Llk_b<-cLLk(Z=Z[[i]],allprobs,gamma=Gamma[[i]])
    # acceptance-rejection
    MHR<-(P(K[i-1]+1)/P(K[i-1]))*exp(Llk_b-Llk0)*(1/L)*(prod(dlogga(a_t[1:(N*nb)],shape=omg,rate=1)/dlogga(A_t[[i]],shape=omg,rate=1))*prod(dlogga(a_t[(N*nb+1):(N*nb+N)],shape=omg,rate=1)))*((1-b_k(K[i-1]+1))/(b_k(K[i-1])*H_b(z=knew, inknot=Inknots[[i]], knot=c(a,Inknots[[i]],b), a0, b0)))
    Ja<-rep(NA,N)
    for (k in 1:N) {
      Ja[k]<-A_t[[i]][(1+(k-1)*nb):(k*nb)][kn+4]-A_t[[i]][(1+(k-1)*nb):(k*nb)][kn+3]
    }
    J<-abs(prod(Ja)*(r[1]*r[2])^N) # product of state-specific Js
    if (runif(1)<min(1,MHR*J)) {
      K[i]<-K[i-1]+1
      logllk[i]<-Llk_b
      Inknots[[i]]<-sort(c(Inknots[[i]],knew))  
      A_t[[i]]<-a_t 
      A[[i]]<-as0
      accept_b<-accept_b+1
    } else {
      K[i]<-K[i-1]
      logllk[i]<-Llk0
    } 
  } else {
    kd<-sample(x=1:length(Inknots[[i]]), size=1, replace = TRUE)
    knots<-c(rep(a,4),Inknots[[i]][-kd],rep(b,4))
    kn<-kd-1
    r<-rep(NA,2)
    r[1]<-(Inknots[[i]][kd]-knots[kn+2])/(knots[kn+5]-knots[kn+2])
    r[2]<-(Inknots[[i]][kd]-knots[kn+3])/(knots[kn+6]-knots[kn+3])
    B <- mSpline(x, knots = knots[5:(K[i-1]+3)], degree = 3, intercept = TRUE,Boundary.knots = c(a,b))
    a_t<-rep(NA,N*(nb-1))
    as0<-rep(NA,N*(nb-1))
    for (k in 1:N) {
      a_t[(1+(k-1)*(nb-1)):(k*(nb-1))][1:kd]<-A_t[[i]][(1+(k-1)*nb):(k*nb)][1:kd]
      a_t[(1+(k-1)*(nb-1)):(k*(nb-1))][kd+1]<-(A_t[[i]][(1+(k-1)*nb):(k*nb)][kd+1]-(1-r[1])*a_t[(1+(k-1)*(nb-1)):(k*(nb-1))][kd])/r[1]
      a_t[(1+(k-1)*(nb-1)):(k*(nb-1))][kd+2]<-(A_t[[i]][(1+(k-1)*nb):(k*nb)][kd+2]-(1-r[2])*a_t[(1+(k-1)*(nb-1)):(k*(nb-1))][kd+1])/r[2]
      a_t[(1+(k-1)*(nb-1)):(k*(nb-1))][-(1:(kd+2))]<-A_t[[i]][(1+(k-1)*nb):(k*nb)][(kd+4):nb]
      a_tt<- a_t[(1+(k-1)*(nb-1)):(k*(nb-1))]-max(a_t[(1+(k-1)*(nb-1)):(k*(nb-1))]) #
      as0[(1+(k-1)*(nb-1)):(k*(nb-1))]<-exp(a_tt)/sum(exp(a_tt))
    }
    allprobs<-B%*%matrix(as0,nrow = nb-1, ncol=N)
    Llk_d<-cLLk(Z=Z[[i]],allprobs,gamma=Gamma[[i]])
    if (is.nan(Llk_d)) {Llk_d<--Inf}
    # acceptance-rejection
    MHR<-(P(K[i-1]-1)/P(K[i-1]))*exp(Llk_d-Llk0)*L*(prod(dlogga(a_t,shape=omg,rate=1)/dlogga(A_t[[i]][1:(N*nb-N)],shape=omg,rate=1))/prod(dlogga(A_t[[i]][(N*nb-N+1):(N*nb)],shape=omg,rate=1)))*(b_k(K[i-1]-1))*H_b(z=Inknots[[i]][kd], inknot=Inknots[[i]][-kd], knot=c(a,Inknots[[i]][-kd],b), a0, b0)/(1-b_k(K[i-1]))
    Ja<-rep(NA,N)
    for (k in 1:N) {
      Ja[k]<-a_t[(1+(k-1)*(nb-1)):(k*(nb-1))][kd+3]-a_t[(1+(k-1)*(nb-1)):(k*(nb-1))][kd+2]
    }
    J<-1/abs(prod(Ja)*(r[1]*r[2])^N) #
    if (runif(1)<min(1,MHR*J)) {
      K[i]<-K[i-1]-1
      logllk[i]<-Llk_d
      Inknots[[i]]<-Inknots[[i]][-kd]
      A_t[[i]]<-a_t 
      A[[i]]<-as0
      accept_d<-accept_d+1
    } else {
      K[i]<-K[i-1]
      logllk[i]<-Llk0
    } 
  }
  Z[i]<-list(NULL) ## CHANGED: free burn-in state sequences (storing all 120k needs ~4 GB)
}

# stop adaptive tuning from now on:
Tau_r<-mean(tau_r[(Tb-5000):Tb])
Tau_2<-mean(tau_2[(Tb-5000):Tb])
Tau_o<-mean(tau_o[(Tb-5000):Tb])

for (i in (Tb+1):(T)) {
  inknots<-Inknots[[i-1]]
  nb<-K[i-1]+4
  knots<-c(rep(a,4),inknots,rep(b,4))
  B <- mSpline(x, knots = inknots, degree = 3, intercept = TRUE,Boundary.knots = c(a,b))
  as0<-A[[i-1]]
  omg<-Omg[i-1]
  allprobs<-B%*%matrix(as0,nrow = nb, ncol=N)
  # FFBS 
  FF<-matrix(NA,N,n)
  FF[,1]<-delta*allprobs[1,]/sum(delta*allprobs[1,])
  for (k in 2:n) {
    for (j in 1:N) {
      FF[j,k]<-allprobs[k,j]*sum(FF[,k-1]*Gamma[[i-1]][,j])
    }
    FF[,k]<-FF[,k]/sum(FF[,k])
  }
  Z[[i]]<-rep(NA,n)
  Z[[i]][n]<-sample(1:N,size = 1, replace = T, prob=FF[,n]/sum(FF[,n]))
  for (k in (n-1):1) {Z[[i]][k]<-sample(1:N,size = 1, replace = T, prob=FF[,k]*Gamma[[i-1]][,Z[[i]][k+1]]/sum(FF[,k]*Gamma[[i-1]][,Z[[i]][k+1]]))}
  Z_counts<-matrix(table(factor(paste0(head(Z[[i]],-1),tail(Z[[i]],-1)),levels = c('11','12','13','21','22','23','31','32','33'))),N,N,byrow=T)
  # update transition matrix
  Gamma[[i]]<-matrix(NA,N,N)
  for (k in 1:N) {Gamma[[i]][k,]<-rdirichlet(1, alpha=1+Z_counts[k,])}
  # update r_k
  Llk0<-cLLk(Z=Z[[i]],allprobs,gamma=Gamma[[i]])
  kreloc_ind<-sample(x=1:length(inknots), size=1, replace = TRUE)
  krenew<-rtruncnorm(1, a=a0, b=b0, mean=inknots[kreloc_ind], sd=Tau_r)
  knots<-sort(c(rep(a,4),inknots[-kreloc_ind],krenew,rep(b,4)))
  B <- mSpline(x, knots = knots[5:(K[i-1]+4)], degree = 3, intercept = TRUE,Boundary.knots = c(a,b))
  allprobs<-B%*%matrix(as0,nrow = nb, ncol=N)
  Llk1<-cLLk(Z=Z[[i]],allprobs,gamma=Gamma[[i]])
  AR1[i]<-min(1,exp(Llk1-Llk0)*(dtruncnorm(x=inknots[kreloc_ind], a=a0, b=b0, mean=krenew, sd=Tau_r)/dtruncnorm(x=krenew, a=a0, b=b0, mean=inknots[kreloc_ind], sd=Tau_r)))
  # acceptance-rejection
  if (runif(1)<AR1[i]) {
    Inknots[[i]]<-sort(c(inknots[-kreloc_ind],krenew))
    Llk0<-Llk1
  } else {
    Inknots[[i]]<-Inknots[[i-1]]
  }
  # update A_t (using RW-MH)
  knots<-c(rep(a,4),Inknots[[i]],rep(b,4))
  B <- mSpline(x, knots = knots[5:(K[i-1]+4)], degree = 3, intercept = TRUE,Boundary.knots = c(a,b))
  a_t<-A_t[[i-1]]+rnorm(N*nb,0,Tau_2)
  as1<-rep(NA,N*nb)
  for (k in 1:N) {
    as1[(1+(k-1)*nb):(k*nb)]<-exp(a_t[(1+(k-1)*nb):(k*nb)])/sum(exp(a_t[(1+(k-1)*nb):(k*nb)]))
  }
  allprobs<-B%*%matrix(as1,nrow = nb, ncol=N)
  Llk1<-cLLk(Z=Z[[i]],allprobs,gamma=Gamma[[i]])
  AR2[i]<-min(1,exp(Llk1-Llk0)*(prod(dlogga(a_t, shape=omg, rate=1)/dlogga(A_t[[i-1]], shape=omg, rate=1))))
  # acceptance-rejection
  if (runif(1)<AR2[i]) {
    A_t[[i]]<-a_t
    Llk0<-Llk1
    as0<-as1
  } else {
    A_t[[i]]<-A_t[[i-1]]
  }
  # update omg (hyperprior ga(1,1))
  omg1<-exp(log(omg)+rnorm(1,0,Tau_o))
  AR4[i]<-min(1,(dgamma(omg1, shape=1, rate=1)/dgamma(omg, shape=1, rate=1))*(prod(dlogga(A_t[[i]], shape=omg1, rate=1)/dlogga(A_t[[i]], shape=omg, rate=1)))*(omg1/omg))
  # acceptance-rejection
  if (runif(1)<AR4[i]) {
    Omg[i]<-omg1
    omg<-omg1
  } else {
    Omg[i]<-Omg[i-1]
  }
  A[[i]]<-as0
  ## RJMCMC step
  u<-runif(1)
  # birth step 
  if (u<b_k(K[i-1])) {
    kadd_ind<-sample(x=1:length(Inknots[[i]]), size=1, replace = TRUE)
    knew<-rtruncnorm(1, a=a0, b=b0, mean=Inknots[[i]][kadd_ind], sd=(knots[kadd_ind+5]-knots[kadd_ind+3])^alpha)      
    kn<-which(sort(c(knots,knew))==knew)-5 # 
    r<-rep(NA,2) # shared across states
    r[1]<-(knew-knots[kn+2])/(knots[kn+5]-knots[kn+2])
    r[2]<-(knew-knots[kn+3])/(knots[kn+6]-knots[kn+3])
    knots<-sort(c(knots,knew))
    B <- mSpline(x, knots = knots[5:(K[i-1]+5)], degree = 3, intercept = TRUE,Boundary.knots = c(a,b))
    a_t<-rep(NA,N*(nb+1))
    as0<-rep(NA,N*(nb+1))
    for (k in 1:N) {
      u<-runif(1)
      a_t[(1+(k-1)*(nb+1)):(k*(nb+1))][kn+2]<-r[1]*(A_t[[i]][(1+(k-1)*nb):(k*nb)][kn+2])+(1-r[1])*(A_t[[i]][(1+(k-1)*nb):(k*nb)][kn+1])
      a_t[(1+(k-1)*(nb+1)):(k*(nb+1))][kn+3]<-r[2]*(A_t[[i]][(1+(k-1)*nb):(k*nb)][kn+3])+(1-r[2])*(A_t[[i]][(1+(k-1)*nb):(k*nb)][kn+2])
      a_t[(1+(k-1)*(nb+1)):(k*(nb+1))][kn+4]<-u*(A_t[[i]][(1+(k-1)*nb):(k*nb)][kn+4])+(1-u)*(A_t[[i]][(1+(k-1)*nb):(k*nb)][kn+3])
      a_t[(1+(k-1)*(nb+1)):(k*(nb+1))][-c(kn+2,kn+3,kn+4)]<-A_t[[i]][(1+(k-1)*nb):(k*nb)][-c(kn+2,kn+3)]
      a_tt<-a_t[(1+(k-1)*(nb+1)):(k*(nb+1))]-max(a_t[(1+(k-1)*(nb+1)):(k*(nb+1))])
      as0[(1+(k-1)*(nb+1)):(k*(nb+1))]<-exp(a_tt)/sum(exp(a_tt))
    }
    allprobs<-B%*%matrix(as0,nrow = nb+1, ncol=N)
    Llk_b<-cLLk(Z=Z[[i]],allprobs,gamma=Gamma[[i]])
    # acceptance-rejection
    MHR<-(P(K[i-1]+1)/P(K[i-1]))*exp(Llk_b-Llk0)*(1/L)*(prod(dlogga(a_t[1:(N*nb)],shape=omg,rate=1)/dlogga(A_t[[i]],shape=omg,rate=1))*prod(dlogga(a_t[(N*nb+1):(N*nb+N)],shape=omg,rate=1)))*((1-b_k(K[i-1]+1))/(b_k(K[i-1])*H_b(z=knew, inknot=Inknots[[i]], knot=c(a,Inknots[[i]],b), a0, b0)))
    Ja<-rep(NA,N)
    for (k in 1:N) {
      Ja[k]<-A_t[[i]][(1+(k-1)*nb):(k*nb)][kn+4]-A_t[[i]][(1+(k-1)*nb):(k*nb)][kn+3]
    }
    J<-abs(prod(Ja)*(r[1]*r[2])^N) # product of state-specific Js
    if (runif(1)<min(1,MHR*J)) {
      K[i]<-K[i-1]+1
      logllk[i]<-Llk_b
      Inknots[[i]]<-sort(c(Inknots[[i]],knew))  
      A_t[[i]]<-a_t 
      A[[i]]<-as0
      accept_b<-accept_b+1
    } else {
      K[i]<-K[i-1]
      logllk[i]<-Llk0
    } 
  } else {
    kd<-sample(x=1:length(Inknots[[i]]), size=1, replace = TRUE)
    knots<-c(rep(a,4),Inknots[[i]][-kd],rep(b,4))
    kn<-kd-1
    r<-rep(NA,2)
    r[1]<-(Inknots[[i]][kd]-knots[kn+2])/(knots[kn+5]-knots[kn+2])
    r[2]<-(Inknots[[i]][kd]-knots[kn+3])/(knots[kn+6]-knots[kn+3])
    B <- mSpline(x, knots = knots[5:(K[i-1]+3)], degree = 3, intercept = TRUE,Boundary.knots = c(a,b))
    a_t<-rep(NA,N*(nb-1))
    as0<-rep(NA,N*(nb-1))
    for (k in 1:N) {
      a_t[(1+(k-1)*(nb-1)):(k*(nb-1))][1:kd]<-A_t[[i]][(1+(k-1)*nb):(k*nb)][1:kd]
      a_t[(1+(k-1)*(nb-1)):(k*(nb-1))][kd+1]<-(A_t[[i]][(1+(k-1)*nb):(k*nb)][kd+1]-(1-r[1])*a_t[(1+(k-1)*(nb-1)):(k*(nb-1))][kd])/r[1]
      a_t[(1+(k-1)*(nb-1)):(k*(nb-1))][kd+2]<-(A_t[[i]][(1+(k-1)*nb):(k*nb)][kd+2]-(1-r[2])*a_t[(1+(k-1)*(nb-1)):(k*(nb-1))][kd+1])/r[2]
      a_t[(1+(k-1)*(nb-1)):(k*(nb-1))][-(1:(kd+2))]<-A_t[[i]][(1+(k-1)*nb):(k*nb)][(kd+4):nb]
      a_tt<- a_t[(1+(k-1)*(nb-1)):(k*(nb-1))]-max(a_t[(1+(k-1)*(nb-1)):(k*(nb-1))]) #
      as0[(1+(k-1)*(nb-1)):(k*(nb-1))]<-exp(a_tt)/sum(exp(a_tt))
    }
    allprobs<-B%*%matrix(as0,nrow = nb-1, ncol=N)
    Llk_d<-cLLk(Z=Z[[i]],allprobs,gamma=Gamma[[i]])
    if (is.nan(Llk_d)) {Llk_d<--Inf}
    # acceptance-rejection
    MHR<-(P(K[i-1]-1)/P(K[i-1]))*exp(Llk_d-Llk0)*L*(prod(dlogga(a_t,shape=omg,rate=1)/dlogga(A_t[[i]][1:(N*nb-N)],shape=omg,rate=1))/prod(dlogga(A_t[[i]][(N*nb-N+1):(N*nb)],shape=omg,rate=1)))*(b_k(K[i-1]-1))*H_b(z=Inknots[[i]][kd], inknot=Inknots[[i]][-kd], knot=c(a,Inknots[[i]][-kd],b), a0, b0)/(1-b_k(K[i-1]))
    Ja<-rep(NA,N)
    for (k in 1:N) {
      Ja[k]<-a_t[(1+(k-1)*(nb-1)):(k*(nb-1))][kd+3]-a_t[(1+(k-1)*(nb-1)):(k*(nb-1))][kd+2]
    }
    J<-1/abs(prod(Ja)*(r[1]*r[2])^N) #
    if (runif(1)<min(1,MHR*J)) {
      K[i]<-K[i-1]-1
      logllk[i]<-Llk_d
      Inknots[[i]]<-Inknots[[i]][-kd]
      A_t[[i]]<-a_t 
      A[[i]]<-as0
      accept_d<-accept_d+1
    } else {
      K[i]<-K[i-1]
      logllk[i]<-Llk0
    } 
  }
  if ((i-Tb-1) %% 10 != 0) Z[i]<-list(NULL) ## CHANGED: free state sequences dropped by the thinning below (memory)
}

end.time<-Sys.time()
runtime<-round(difftime(end.time, start.time, units = "secs"),2)



################################################################




#
K_trace <- K ## CHANGED: keep the full trace of K (incl. burn-in)
nthin <- 10
Gamma <- Gamma[-(1:Tb)]
Gamma <- Gamma[rep(c(TRUE, FALSE), c(1, nthin-1))] # thin = 10
Z <- Z[-(1:Tb)]
Z <- Z[rep(c(TRUE, FALSE), c(1, nthin-1))]
A <- A[-(1:Tb)]
A <- A[rep(c(TRUE, FALSE), c(1, nthin-1))]
Inknots <- Inknots[-(1:Tb)]
Inknots <- Inknots[rep(c(TRUE, FALSE), c(1, nthin-1))]
K <- K[-(1:Tb)]
K <- K[rep(c(TRUE, FALSE), c(1, nthin-1))]
Tp<-length(K) # number of samples for inference
## CHANGED: removed preliminary density plot (repeated below with credible bands)

ipm<-matrix(c(1:N^2),nrow=N, byrow = T)
pm<-c(1,2,3) 
spm<-order(pm)
for (i in 1:Tp) {
  Gamma[[i]]<-matrix(as.vector(t(Gamma[[i]]))[as.vector(t(ipm[spm,spm]))],N,N,byrow=T)
}

pZ<-vector("list")
for (i in 1:Tp) {
  pZ[[i]]<-rep(NA,n)
  for (k in 1:N) {pZ[[i]][which(Z[[i]]==k)]<-pm[k]}
}

for (i in 1:Tp) {
  A[[i]]<-as.vector(matrix(A[[i]],nrow = K[i]+4, ncol=N)[,spm])
}

table(K)

# decoding  
stateprob<-matrix(NA,N,n)
for (i in 1:N) {
  for (j in 1:n){
    sp<-rep(NA,Tp)
    for (k in 1:(Tp)) {sp[k]<-1*(pZ[[k]][j]==i)}
    stateprob[i,j]<-sum(sp)/(Tp)
  }
}
ld<-rep(NA,n)
for (i in 1:n){ld[i]<-which.max(stateprob[,i])}

sp2<-seq(a,b,length=5000)
dsp<-vector("list")
edsp<-vector("list")
edsp_L95<-vector("list")
edsp_U95<-vector("list")
for (j in 1:N) {
  for (i in 1:(Tp)){
    B <- mSpline(sp2, knots = Inknots[[i]], degree = 3, intercept = TRUE, Boundary.knots = c(a,b))
    dsp[[i]]<-B%*%A[[i]][(1+(j-1)*(K[i]+4)):(j*(K[i]+4))]
  }
  edsp[[j]] <- rep(NA,length(sp2))
  edsp_L95[[j]] <- rep(NA,length(sp2))
  edsp_U95[[j]] <- rep(NA,length(sp2))
  for (i in 1:length(sp2)) {edsp[[j]][i]<-mean(sapply(dsp, `[[`, i))}
  for (i in 1:length(sp2)) {edsp_L95[[j]][i]<-quantile(sapply(dsp, `[[`, i),c(0.025))}
  for (i in 1:length(sp2)) {edsp_U95[[j]][i]<-quantile(sapply(dsp, `[[`, i),c(0.975))}
}

## CHANGED: removed plot with shark-specific axis limits
Gamma_hat <- Reduce("+", Gamma)/length(Gamma)
pi<-statdistr(Gamma_hat)
sdgamma<-rep(NA,N^2)
for (j in 1:(N^2)){
  gamma0<-rep(NA,Tp)
  for (i in 1:Tp) {
    gamma0[i]<-Gamma[[i]][j]
  }
  sdgamma[j]<-sd(gamma0)
}
sdgamma1<-matrix(nrow=N,ncol=N,sdgamma)
## CHANGED: removed plots with shark-specific axis limits and the per-chain output for
## their multi-chain PSRF script; the comparison with qREML follows below


# Comparison with qREML ---------------------------------------------------

## computation time (their runtime covers the MCMC only)
c(qREML_sec = as.numeric(time_qreml, units = "secs"),
  CF_init_sec = as.numeric(init.time), CF_MCMC_sec = as.numeric(runtime),
  CF_total_hours = (as.numeric(init.time) + as.numeric(runtime)) / 3600)

## number of interior knots: trace incl. burn-in and posterior after burn-in
plot(K_trace, type = "s", xlab = "Iteration", ylab = "K", bty = "n")
abline(v = Tb, lty = 2)
table(K) / Tp

## state labels are fixed in their code (pm = 1:3): check that states are ordered as in qREML
sapply(edsp, function(f) sum(f * sp2) / sum(f)) # posterior mean density: state means
colSums(Dens * xseq) / colSums(Dens) # qREML: state means

## transition probabilities and stationary distribution
round(Gamma_hat, 3); round(mod$Gamma, 3)
round(pi, 3); round(mod$delta, 3)

## state-dependent densities: Hellinger distance between posterior mean and qREML estimate
Dens_CF = sapply(1:N, function(j) approx(sp2, edsp[[j]], xseq)$y)
sqrt(0.5 * colSums((sqrt(Dens_CF) - sqrt(Dens))^2) * diff(xseq[1:2]))

## decoded states: proportion of agreement
states_local_qreml = apply(stateprobs(mod = mod), 1, which.max)
allprobs_CF = sapply(1:N, function(j) approx(sp2, edsp[[j]], x)$y) # posterior mean densities at data
states_CF = viterbi(as.vector(pi), Gamma_hat, allprobs_CF) # plug-in global decoding
c(their_local_vs_qreml_viterbi = mean(ld == states_qreml),
  local_vs_local = mean(ld == states_local_qreml),
  viterbi_vs_viterbi = mean(states_CF == states_qreml))

## visualizing results
pdf("./case_studies/figs/caracara_chen.pdf", width = 8.5, height = 4)
par(mfrow = c(1,2), mar = c(5,4,3.5,1) + 0.1)
hist(data$logVDBA, prob = TRUE, breaks = 50, bor = "white",
     main = "", xlab = "log(VeDBA)", ylab = "Density")
for(j in 1:N){
  lines(xseq, mod$delta[j] * Dens[,j], col = color[j], lwd = 2)
  lines(sp2, pi[j] * edsp[[j]], col = color[j], lwd = 2, lty = 2)
}
legend("topright", legend = c("Resting", "Low activity", "High activity"), 
       col = color, lwd = 2, bty = "n")

plot(NA, xlim = range(xseq), ylim = c(0, 1.3 * max(unlist(edsp_U95))), bty = "n",
     xlab = "log(VeDBA)", ylab = "State-dependent density")
for(j in 1:N){
  polygon(c(sp2, rev(sp2)), c(edsp_L95[[j]], rev(edsp_U95[[j]])), 
          col = alpha(color[j], 0.3), border = NA)
  lines(xseq, Dens[,j], col = color[j], lwd = 2)
  lines(sp2, edsp[[j]], col = color[j], lwd = 2, lty = 2)
}
legend("topright", legend = c("qREML", "Bayesian (posterior mean)", "95% credible band"),
       lwd = c(2, 2, 8), lty = c(1, 2, 1), col = c(1, 1, alpha(1, 0.3)), bty = "n")
dev.off()
