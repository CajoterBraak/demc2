#  Adaptive MCMC via demc
rm(list=ls(all=TRUE))

#library(LaplacesDemon) # or geoR for dinvchisq (and rinvchisq, rmvn)
library(geoR)
library(demc)
# generate some data

set.seed(13451)
n= 20
rho = -0.8
# generate two correlated predictors
Sigma = matrix(c(1,rho,rho,1),nrow = 2,ncol=2)
X = rmvn(n, mu = rep(0, nrow(Sigma)), Sigma = Sigma)
Beta = c(1,2);sigma = 0.5
y = X%*%Beta + sigma*rnorm(n)
Data = data.frame(y,X); names(Data)
# 
logposterior_sigma2 <- function(theta, data = Data){
  # make data variates local 
  y = data$y; x1= data$X1 ; x2 = data$X2
  # make the parameters available
  b1 = theta[1]
  b2 = theta[2]
  sigma2 = theta[3] # sigma squared
  if (sigma2 >0) {
    # prior
    logprior = dnorm(b1, 0, sqrt(v1), log= TRUE)+dnorm(b2,0,sqrt(v2), log= TRUE)+geoR::dinvchisq(sigma2, nu, s2, log= TRUE) 
    # likelihood 
    loglikelihood = sum(dnorm(y, mean = b1*x1+b2*x2, sd = sqrt(sigma2),log= TRUE))
    logpost = logprior + loglikelihood  
  } else logpost = -Inf # because logprior = -Inf
  logpost
} # end of function logposterior

logposterior_logsigma2 <- function(theta, data = Data){ 
  theta_original = theta
  theta_original[3] =  exp(theta[3]) #    sigma2 
  logposterior_sigma2(theta_original, data = data) + theta[3]   # theta[3] = log(sigma2)
  # the second term stems from the Jacobian, slide LL122, section 2.9 of the LL book
}


logPosterior <- logposterior_logsigma2
ls()
v1= 10; v2= 5; # fixed values for the parametes in the prior for beta1 and beta2
nu = 2;  s2 = 3 # fixed values for the parameters in the prior for sigma2

tau2 = s2
#Ex DEMC.1
c(v1,v2,nu,tau2) # parameters of priors
d = 3 # number of parameter
(Nchain = 2*d) # 
# DEMC needs an initial population of Nchain samples,
# for example, sample an initial population from the prior
b1 = b1_init = rnorm(Nchain,0,sqrt(v1))
b2 = b2_init = rnorm(Nchain,0,sqrt(v2))
logsigma2 = logsigma2_init = log(rinvchisq(Nchain,nu,tau2))
Z = rbind(b1,b2, logsigma2)

# check whether the logPosterior function works
logPosterior(Z[,1],data = Data)
#Ex DEMC.2
N.iter = 6000

Draws <- demc(Z, logPosterior, data=Data, n.generation = floor(N.iter/Nchain), n.thin = 1, n.burnin =0 )
summary(Draws)
summary(Draws, keep.all = FALSE)


# DEMCzs needs a large initial population
# start from fresh 
nval = 10*d # 10* number of parameters
b0 = b0_init = rnorm(nval,0,sqrt(v1))
b1 = b1_init = rnorm(nval,0,sqrt(v2))
logsigma2 = logsigma2_init = log(rinvchisq(nval,nu,tau2))
Z = rbind(b0,b1, logsigma2)
# or use result from previos
nZ = ncol(Draws$Draws)
Z_index = floor(seq(from =nZ/2, to = nZ, length.out = nval ))
Z = Draws$Draws[,Z_index]


Nchain = 3
Draws2 <- demc_zs(Nchain, Z, logPosterior, data=Data, n.generation = floor(N.iter/Nchain), n.thin = 1, n.burnin = 0 )
#monitor.DE.MC(Draws.DEMCZS3$Draws, keep.all = TRUE , m= Nchain)
summary(Draws2, keep.all = TRUE)


# Ex DEMC.5
Nchain = 3
blocks = list(); 
# three blocks
blocks[[1]]=  1 
blocks[[2]] = 2
blocks[[3]] = 3
Draws3 <- demc_zs(Nchain, Z, logPosterior, blocks = blocks, data=Data, n.generation = floor(N.iter/Nchain), n.thin = 1, n.burnin = 0 )
#monitor.DE.MC(Draws.DEMCZS3$Draws, keep.all = TRUE , m= Nchain)
summary(Draws3, keep.all = TRUE)


# Ex DEMC.5
Nchain = 3
blocks = list(); 
# two blocks
# beta1 and beta2 need to be in one block as they are highly correlated
blocks[[1]]=  c(1,2) 
blocks[[2]] = 3
Draws4 <- demc_zs(Nchain, Z, logPosterior, blocks = blocks, data=Data, n.generation = floor(N.iter/Nchain), n.thin = 1, n.burnin = 0 )
#monitor.DE.MC(Draws.DEMCZS2$Draws, keep.all = TRUE , m= Nchain)
summary(Draws4, keep.all = TRUE)



blocks[[1]]=  c(1,3) 
blocks[[2]] = 2
Draws4B <- demc_zs(Nchain, Z, logPosterior, blocks = blocks, data=Data, n.generation = floor(N.iter/Nchain), n.thin = 1, n.burnin = 0 )
#monitor.DE.MC(Draws.DEMCZS2B$Draws, keep.all = TRUE , m= Nchain)
summary(Draws4B, keep.all = TRUE)

# 
# How can check you this idea numerically?
# look at accept rates per block, convergence, autocorrelation or Rhat (Brooks-Gelman)
# Perform the check.
# The accept rates for the beta1&2 block is higher than when separated
# the rates are close to optimal for two parameters (d = 2)
# I did not calculate auto correlatios yet
# The Rhat statistics for the three blocks are lower (even <1) than with two blocks
# This suggests that three blocks may be doing fine! (This is unexpected, and probably not the right conclusion)

C3= as.coda.demc(Draws)
C2= demc2coda(Draws.DEMCZS2)
C2B= demc2coda(Draws.DEMCZS2B)

gelman.plot(C3,  max.bins = 20,confidence = 0.95, transform = FALSE, autoburnin=FALSE, auto.layout = TRUE)
gelman.plot(C2,  max.bins = 20,confidence = 0.95, transform = FALSE, autoburnin=FALSE, auto.layout = TRUE)
# note the different vertical scale in next plot
gelman.plot(C2B,  max.bins = 20,confidence = 0.95, transform = FALSE, autoburnin=FALSE, auto.layout = TRUE)
gelman.diag(C3, mult =TRUE)
gelman.diag(C2, mult =TRUE)
gelman.diag(C2B, mult =TRUE)

C3 = window(C3, start= 1000)
C2 = window(C2, start= 1000)
C2B = window(C2B, start= 1000)

#plot(C3)
acfplot(C3)
acfplot(C2)
acfplot(C2B)
