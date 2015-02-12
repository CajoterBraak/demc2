#' @title generates MCMC samples using Differential Evolution Markov Chain (ter Braak 2008)
#'
#' @description
#' \code{demc} implements multi-chain adaptive MCMC on real parameter spaces via
#'  Differential Evolution Markov Chain. It allows restart from a previous run.
#' It is the only (?) adaptive MCMC method that is really Markovian and
#' not just ergodic. Required input: starting position for each chain (X) 
#' and a unnormalized logposterior function (FUN).
#'
#'@param X matrix of initial values or \code{demc} object resulting from a previous run. 
#'If matrix, initial parameter values for N chains (columns) for each of the parameters (rows).
#' See Details for choice of N. 
#'@param FUN function specifying the (unnormalized) logposterior or target. 
#'FUN(theta ,...) with theta a d vector of parameter values and ...  other arguments to FUN (e.g. NULL, data, ..). The value of FUN should a single numeric value.
#'@param blocks list of sets of parameters, where each set contains the
#'indices (or names if X has rownames)  of variables in theta that are updated jointly,
#' e.g.
#' blocks= list(); blocks[[1]] = c(1,3); blocks[[2]] = c(2,4). The default uses a single block (joint update of all parameters).
#'@param f value specifying the scale of the updates. 
#' The scale used is f/sqrt(2k), where k the number of parameters in a block 
#' (k=d if there only one block); 
#' f can be a vector of length(blocks) to set differential scaling factors for blocks. 
#'@param n.generation integer, number of generations, each one generating ncol(X) samples.
#'@param n.thin integer, thinning number, e.g 3 stores every 3rd generation (default 1, no thinning).
#'@param n.burnin integer, number of initial generations that is not stored (default 0).
#'@param eps.add real value . Standard deviation of the (additive) independent random normal step added to the DE update. 
#' that guarantees that all parameter values can be reached.
#' Must be postive to garantee that all positions in the space can be reached.
#' For targets (posteriors) without gaps, it can set small or even to 0 (default 0). 
#'@param p.f.is.1 probability (default 0.1) of using scale 0.98 in the parallel update instead of f/sqrt(2d). Allows exploration of multi-modal posteriors.
#'@param verbose logical (default FALSE). No used currently
#'@param logfitness_X Can be missing. If set, logposterior values for the columns of X; can save some computation if known.
#'@return an S3 object of class \code{demc} which is a list with 
#'
#' $Draws d x (Nchain*n) matrix  with n the number of retained simulations  [post process with summary or as.coda].
#' matrix(Draws[p,],nrow= Nchain) is an Nchain x n array (for each parameter p in 1:d)
#' 
#' $accept.prob.mat accept probability matrix of blocks by chains. 
#' 
#' $X.final matrix of parameter values of the last generation (same shape as initial X).
#' 
#' $logfitness.X.final vector of logposterior values for columns of X.final (useful for restarting).
#' 
#' $Nchain number of chains.
#' 
#' $demc_zs  logical set of FALSE.
#'@details \code{demc} This function implements the method of ter Braak (2006) using the differential evolution update (also known as the paralled direction sampling update). 
#' See also \code{demc_zs} for an egodic variant with learning from the past. 
#' The number of initial positions (N, columns of X) should be larger than the number of paramaters (d) in each block. The advice is N=2d. So fewer starting positions are required by specifying blocks with few (or even single) parameters in each block. 
#'@references ter Braak, C. J. F. (2006). A Markov Chain Monte Carlo version of the genetic algorithm Differential Evolution: easy Bayesian computing for real parameter spaces. Statistics and Computing, 16, 239-249.\url{http://edepot.wur.nl/26336}
#'



demc <- function(X, FUN, blocks, f = 2.38, n.generation = 1000, n.thin = 1, n.burnin = 0, eps.add = 0, p.f.is.1= 0.1, verbose = FALSE,logfitness_X, ...){
if ("demc"%in%class(X)) {
  is.update = TRUE
  X = X$X.final
  logfitness_X = logfitness.X.final
} else {
  if(!(is.matrix(X)&& is.numeric(X))) stop("demc: X must be a numeric matrix")
  is.update = FALSE
}
if (missing(blocks)) blocks= list(seq_len(nrow(X)))
if (nrow(X)> max(sapply(blocks, length))) {
  stop("demc: ncol of initial population (", 
       ncol(X),") smaller than maximum block size (",max(sapply(blocks, length)),")\n Increase the number of starting positions (columns) or specify (smaller) blocks.")
} 
Npop = ncol(X)
Npar = nrow(X)
chainset = seq_len(Npop)
Nblock = length(blocks)
F2 =  matrix(abs(f),nrow=1,ncol=Nblock)
F1 = rep(0.98,Nblock)
Draws = NULL
accept = matrix(rep(NA,n.generation*Nblock),nrow = n.generation, ncol= Nblock)
if (missing(logfitness_X) && !is.update) logfitness_X = apply (X, 2, FUN, ...)
Naccept = matrix(0,nrow = Nblock, ncol = Npop, dim=list(blocks=paste("block", seq_len(Nblock)), chains=paste("chain", chainset) ))
for (iter in 1:n.generation) {
   for (i in chainset){
     for (iblock in 1:Nblock)  {
     # set of parameters to update
     parset = blocks[[iblock]]
     dsub = length(parset)  # dsubspace dimension of sampling
      # select to random different individuals (and different from i) in rr, a 2-vector
     rr = sample(chainset[-i], 2, replace = FALSE)
     if ( runif(1) < p.f.is.1 )  F = F1 else  F2/sqrt(2*dsub)
     x_prop = X[,i]
     x_prop_sub = X[parset,i] + F[iblock]*(X[parset,rr[1]]-X[parset,rr[2]])  +  eps.add*rnorm(dsub,0,1)
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
 out = list(Draws= Draws, accept.prob.mat= Naccept/n.generation, X.final = X, logfitness.X.final = logfitness_X, Nchain=Npop, demc_zs = FALSE)
 class(out) <- c("demc")
 invisible(out)
}
