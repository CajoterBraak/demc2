demc <- function(X, FUN, blocks = list(seq_len(nrow(X))), f = -2.38, n.generation = 1000, n.thin = 1, n.burnin = 0, eps = 0, ...){
# Differential Evolution Markov Chain applied to X with logposterior specified by FUN
# X is the initial population: a matrix of number of parameters by number of individuals (k x Npop)
# value
# FUN(theta ,...) with theta a k vector of parameter values and ... are other arguments to FUN (e.g. NULL, data, ..)
# blocks: list of sets of indices for parameters, where each set contains the
#              indices of variables in theta that are updated jointly
#              If all elements of theta need to be sampled, 
#              the union of sets needs to be equal to 1:k
#              So, for  element iblock of the list,  
#              blocks[[iblock]] contains indices of parameters that are jointly updated,
#              e.g. blocks= list(); blocks[[1]] = c(1,3); blocks[[2]] = c(2,4)
# f  scale factor for d = 1; for d>1 divided by sqrt(2d); if negative, scale is set to 1 (0.98 in fact) to allow for multimodality.
# Value
# $Draws k x (Npop*n) array  with n the number of retained simulations  [post process with monitor.DE.MC]
#         matrix(Draws[p,],nrow= Npop) is an Npop x n array (for each p in 1:k)
# $ accept.prob
# $ X.final
# $ logfitness.X.final
# Reference:
# ter Braak, C. J. F. (2006). A Markov Chain Monte Carlo version of the genetic algorithm Differential Evolution:
# easy Bayesian computing for real parameter spaces. Statistics and Computing, 16, 239-249.
Npop = ncol(X)
Npar = nrow(X)
chainset = seq_len(Npop)
Nblock = length(blocks)
F2 =  matrix(abs(f),nrow=1,ncol=Nblock)
if (min(f)>0) F1 = F2 else F1 = rep(0.98,Nblock)
Draws = NULL
accept = matrix(rep(NA,n.generation*Nblock),nrow = n.generation, ncol= Nblock)
logfitness_X = apply (X, 2, FUN, ...)
Naccept = matrix(0,nrow = Nblock, ncol = Npop, dim=list(blocks=paste("block", seq_len(Nblock)), chains=paste("chain", chainset) ))
for (iter in 1:n.generation) {
   for (i in chainset){
     for (iblock in 1:Nblock)  {
     # set of parameters to update
     parset = blocks[[iblock]]
     dsub = length(parset)  # dsubspace dimension of sampling
      # select to random different individuals (and different from i) in rr, a 2-vector
     rr = sample(chainset[-i], 2, replace = FALSE)
     if (iter%%10) F = F2/sqrt(2*dsub) else F = F1
     x_prop = X[,i]
     x_prop_sub = X[parset,i] + F[iblock]*(X[parset,rr[1]]-X[parset,rr[2]])  +  eps*rnorm(dsub,0,1)
     x_prop[parset] = x_prop_sub
     logfitness_x_prop = FUN(x_prop,  ...)
     if ((logfitness_x_prop - logfitness_X[i] )> log(runif(1)) ){
        Naccept[iblock,i] = Naccept[iblock,i]+1
        X[,i] = x_prop
        logfitness_X[i] = logfitness_x_prop
     }
     } # iblock
  } # i loop
  if (!(iter%%n.thin) && iter > n.burnin) {  # retain sample
     Draws = cbind(Draws,X)
  }
} # n.generation
 out = list(Draws= Draws, accept.prob.mat= Naccept/n.generation, X.final = X, logfitness.X.final = logfitness_X, Nchain=Npop)
 class(out) <- c("demc")
 invisible(out)
}
