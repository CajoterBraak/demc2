#  Adaptive MCMC via demc
rm(list=ls(all=TRUE))

library(LaplacesDemon) # or geoR for dinvchisq (and rinvchisq, rmvn)
library(geoR)
library(demc2)
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
v1 = 10; v2 = 5; # fixed values for the parametes in the prior for beta1 and beta2
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

N.iter = 6000
Draws0 <- demc(Z, logposterior_sigma2, data=Data, n.generation = floor(N.iter/Nchain), n.thin = 1, n.burnin =0 )
summary(Draws0)
summary(Draws0, keep.all = FALSE)


# Recommendation: use parameters with unbounded support
logsigma2 = log(sigma2_init)
Z = rbind(b1,b2, logsigma2)
logposterior_logsigma2
logPosterior <- logposterior_logsigma2
logPosterior(Z[,1],data = Data)
#Ex DEMC.2
N.iter = 60000
Draws1 <- demc(Z, logPosterior, data=Data, n.generation = floor(N.iter/Nchain), n.thin = 1, n.burnin =0 )
summary(Draws1)
dim(Draws1$Draws)
#summary(Draws1, keep.all = FALSE)
C1= as.coda(Draws1)
library(coda)
gelman.plot(C1,  max.bins = 20,confidence = 0.95, transform = FALSE, autoburnin=FALSE, auto.layout = TRUE)
C1.2 = window(C1, start= N.iter/Nchain/2) # start = 500 from the 1000 generations
acfplot(C1.2)

## DEMCzs needs a large initial population
## start from fresh
## NOT run
# nval = 10*d ## 10* number of parameters
# b0 = b0_init = rnorm(nval,0,sqrt(v1))
# b1 = b1_init = rnorm(nval,0,sqrt(v2))
# logsigma2 = logsigma2_init = log(rinvchisq(nval,nu,tau2))
# Z = rbind(b0,b1, logsigma2)
##
##
## or use result from previous run, either explicitly
nval = 10*d ## 10* number of parameters
nZ = ncol(Draws1$Draws)
Z_index = floor(seq(from =nZ/2, to = nZ, length.out = nval ))
Z = Draws1$Draws[,Z_index]
Nchain = 3
set.seed(101)
system.time(Draws2.0 <- demc_zs(Nchain, Z, logPosterior, data=Data, n.generation = floor(N.iter/Nchain), n.thin = 1, n.burnin = 0 ))
summary(Draws2.0, keep.all = TRUE)




#################
### Parallel ####
#################
require(doParallel)

# By default, it runs in serial (result not exactly demc_zs because all random numbers are generated beforehand)
nval = 10*d ## 10* number of parameters
nZ = ncol(Draws1$Draws)
Z_index = floor(seq(from =nZ/2, to = nZ, length.out = nval ))
Z = Draws1$Draws[,Z_index]
Nchain = 3
set.seed(101)
system.time(Draws2.0 <- demc_zs_p(Nchain, Z, logPosterior, data=Data, n.generation = floor(N.iter/Nchain), n.thin = 1, n.burnin = 0 ))
summary(Draws2.0, keep.all = TRUE)



# If we want parallelism, we need to create our own cluster
nval = 10*d ## 10* number of parameters
nZ = ncol(Draws1$Draws)
Z_index = floor(seq(from =nZ/2, to = nZ, length.out = nval ))
Z = Draws1$Draws[,Z_index]
Nchain = 3
cl <- makeCluster(3)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("logPosterior","v1","v2","nu","tau2"))
set.seed(101)
Draws2.0 <- demc_zs_p(Nchain, Z, logPosterior, data=Data, n.generation = floor(N.iter/Nchain),
                      n.thin = 1, n.burnin = 0, pClust = cl)
summary(Draws2.0, keep.all = TRUE)
stopCluster(cl)
