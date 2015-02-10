#  Adaptive MCMC via demc
rm(list=ls(all=TRUE))

library(LaplacesDemon) # or geoR for dinvchisq (and rinvchisq, rmvn)
library(geoR)
library(demc)
# example of Bayesian linear regresssion without intercept
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
# parameters of priors (used with these names in the function defining the logposterior)
v1= 10; v2= 5; # fixed values for the parametes in the prior for beta1 and beta2
nu = 2;  tau2 = 3 # fixed values for the parameters in the prior for sigma2
# end of set-up of Bayesian model and data for linear regression

# define a demc run to find the posterior of the parameters vector (b1,b2,sigma2)
# step 1: write a function that gives 
# the value of the logposterior for each value of the parmeter vector
# here it is predefined; see what it looks like
logposterior_sigma2
# check whether the logposterior function works
logposterior_sigma2(c(1,1,1),data = Data)
logposterior_sigma2(c(1,1,-1),data = Data) # -Inf as it should

# step 2: define the number of parallel chains
d = 3 # number of parameters (=length of parameter vector)
(Nchain = 2*d) # 
#  DEMC needs an initial population Z of Nchain samples
# step 3: generate initial population as start of each chain
# for example, sample an initial population from the prior
b1 = b1_init = rnorm(Nchain,0,sqrt(v1))
b2 = b2_init = rnorm(Nchain,0,sqrt(v2))
sigma2 = sigma2_init = rinvchisq(Nchain,nu,tau2)
Z = rbind(b1,b2, sigma2)
# note that 

N.iter = 6000
Draws0 <- demc(Z, logposterior_sigma2, data=Data, n.generation = floor(N.iter/Nchain), n.thin = 1, n.burnin =0 )
summary(Draws0)
summary(Draws0, keep.all = FALSE)



logsigma2 = log(sigma2_init)
Z = rbind(b1,b2, logsigma2)

# check whether the logPosterior function works
logPosterior <- logposterior_logsigma2
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


C3= as.coda(Draws)
C2= as.coda(Draws4)
C2B= as.coda(Draws4B)

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

acfplot(C3)
acfplot(C2)
acfplot(C2B)
