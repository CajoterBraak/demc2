#' @title generates MCMC samples using Differential Evolution Markov Chain (ter Braak & Vrugt 2008)
#'
#' @description
#' \code{demc_zs} implements multi-chain adaptive MCMC on real parameter spaces via
#'  Differential Evolution Markov Chain with learning from the past and the snooker update.
#'  It allows restart from a previous run.
#'  Required input: starting position for each chain (X)
#'  and a unnormalized logposterior function (FUN).
#'
#'@param Nchain number of (parallel) chains.
#'@param Z matrix of initial values or \code{demc} object resulting from a previous run.
#'If matrix, initial parameter values for N chains (columns) for each of the parameters (rows).
#' See Details for choice of N.
#'@param X optional matrix of initial values (d x Nchain); default: X.final of a previous run or the first Nchain columns of Z if Z is a matrix.
#'If matrix, initial parameter values for N chains (columns) for each of the parameters (rows).
#' See Details for choice of N.
#'@param f value specifying the scale of the updates (default 2.38).
#' The scale used is f/sqrt(2k), where k the number of parameters in a block
#' (k=d if there only one block) in the parallel update and k = 1 in the snooker update.
#' f can be a vector of length(blocks) to set differential scaling factors for blocks.
#'@param pSnooker probability of snooker update (default 0.1); 1-pSnooker: prob of the differential evolution (parallel) update
#'@param eps.mult real value (default 0.2). It adds scaled uncorrelated noise
#' to the proposal, by multiplying the scale of the update with a d-dimension uniform
#' variate [1-eps,1+eps].
#' @inheritParams demc
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
#'@details \code{demc_zs} This function implements the method of ter Braak & Vrugt (2008) using learning of the past, the differential evolution (parelle direction) update
#' and the snooker update.
#' The number of initial positions (N, columns of Z) should be much larger than the number of paramaters (d)
#' in each block. The advice is N=10d. So fewer initial values are required by specifying blocks with few (or even single) parameters in each block.
#'
#'@references ter  Braak C. J. F., and Vrugt J. A. (2008). Differential Evolution Markov Chain with snooker updater and fewer chains. Statistics and Computing, 18 (4), 435-446. \url{http://dx.doi.org/10.1007/s11222-008-9104-9}

demc_zs_p <- function(Nchain = 3, Z, FUN, X,
                      blocks, f = 2.38, pSnooker= 0.1, p.f.is.1 = 0.1, n.generation = 1000,
 n.thin, n.burnin = 0, eps.mult =0.2,eps.add = 0, verbose = FALSE, logfitness_X, pClust = NULL,...){
# Differential Evolution Markov Chain applied to X with logposterior specified by FUN
# Z is the initial population: a matrix of number of parameters by number of individuals (d x m0)
# X is the initial active population that will be evolved: a matrix of number of parameters by number of individuals (d x N)
# value
# FUN(theta ,...) with theta a d vector of parameter values and ... are other arguments to FUN (e.g. NULL, data, ..)
# blocks: list of sets of indices for parameters, where each set contains the
#              indices of variables in theta that are updated jointly
#              If all elements of theta need to be sampled,
#              the union of sets needs to be equal to 1:k
#              So, for  element iblock of the list,
#              blocks[[iblock]] contains indices of parameters that are jointly updated,
#              e.g. blocks= list(); blocks[[1]] = c(1,3); blocks[[2]] = c(2,4)
##
# f = scale factor for d =1 (default f=2.38);
#    f can be a tuning vector with one value for each block, if blocks is specified
# eps.mult >0 gives d-dimensional variation around gamma. It adds scaled uncorrelated noise to the proposal.
#             Its advantage over eps.add is that its effect scales with the differences of vectors in the population
#             whereas eps.add does not. if the variance of a dimenstion is close to 0, eps.mult gives smaller changes.
# eps.add >0  is needed to garantee that all positions in the space can be reached.
#             For targets without gaps, it can set small or even to 0.
#
# Value
# $Draws d x (Nchain*n) array  with n the number of retained simulations  [post process with monitor.DE.MC]
#         matrix(Draws[j,],nrow= Nchain) is an Nchain x n array (for each j in 1:d)
# $ accept.prob
# $ X.final
# $ logfitness.X.final
#
#  ter  Braak C. J. F., and Vrugt J. A. (2008). Differential Evolution Markov Chain
#  with snooker updater and fewer chains. Statistics and Computing,
#  18, 435-446.  http://dx.doi.org/10.1007/s11222-008-9104-9 .
#
#    identifier   symbol in paper
#     Nchain          N
#     Npar          d
#     gamma_par     gamma
#     gamma_snooker gamma
#     M0            M0
#     mZ            M
#     n.thin        K
#
# see also
# ter Braak, C. J. F. (2006). A Markov Chain Monte Carlo version of
#    the genetic algorithm Differential Evolution: easy Bayesian computing
#    for real parameter spaces. Statistics and Computing, 16, 239-249.
#
#
#  Cajo ter Braak, Biometris, Wageningen UR, cajo.terbraak@wur.nl   23 October 2008
#
#
if ("demc"%in%class(Z)) {
      is.update = TRUE
      was.demc_zs = Z$demc_zs
      Nchain0 = Z$Nchain
      if (missing(X)) {
        X = Z$X.final
        if (Nchain > Nchain0) {
           Nchain.add = Nchain - Nchain0
           warning("demc_zs: Nchain of previous run was smaller;",Nchain.add, " additional chain-heads sampled from second halve of previous run")

          indx = ncol(Z$Draws)/2  + sample.int(ncol(Z$Draws)/2,Nchain.add)
          X.add = Z$Draws[,indx]
          } else X.add = NULL
        index = 1: min(ncol(X),Nchain)
        X = cbind(X[,index, drop = FALSE],X.add)
        if (is.null(X.add))logfitness_X = Z$logfitness.X.final[index]
        else    logfitness_X = apply (X, 2, FUN, ...)
      }
      Z = Z$Draws
      if (n.burnin>0 && was.demc_zs) warning("n.burnin>0 invalids ergodicity on restarts/updates;/n use n.thin instead")
} else { # no update
      if(!(is.matrix(Z)&& is.numeric(Z))) stop("demc: Z and X must be a numeric matrix")
      is.update = FALSE
      if (nrow(Z)>ncol(Z)) {
        message("demc_zs: nrow of initial population (",
                nrow(Z),") larger than ncol (",ncol(Z),")\n Initial matrix transposed on the assumption that there are ", ncol(Z)," parameters")
        Z = t(Z)
      }
      if (missing(X)){
        X = Z[,1:Nchain,drop=FALSE]
      } else
        if (ncol(X) !=  Nchain){
          if (nrow(X) ==  Nchain && nrow(Z)==ncol(X)){
            message("demc_zs: nrow of initial population X (=",
                    nrow(X),") equal to Nchain (=",Nchain,")\n  Matrix X transposed on the assumption that there are ", ncol(X)," parameters")
            X = t(X)
          } else stop("demc_zs: ncol of initial population X(",
                     nrow(X),") not equal to Nchain (=",Nchain,")\n  and there appear ", nrow(Z)," parameters")
        }
}

if(!(is.matrix(X)&& is.numeric(X))) stop("demc: X must be a numeric matrix")

if (missing(blocks)) blocks= list(seq_len(nrow(X)))
if (missing(logfitness_X)) logfitness_X = apply (X, 2, FUN, ...)

M0 = mZ = ncol(Z)
#Nchain = ncol(X)
Npar = nrow(X)
#Npar12  =(Npar - 1)/2  # factor for Metropolis ratio DE Snooker update
F1 = 0.98
accept = rep(NA,n.generation)
chainset = seq_len(Nchain)
rr = NULL
r_extra = 0
Nblock = length(blocks)
f = matrix(f,nrow=1,ncol=Nblock)
Fs = f/sqrt(2)
Naccept = matrix(0,nrow = Nblock, ncol = Nchain, dim=list(blocks=paste("block", seq_len(Nblock)), chains=paste("chain", chainset) ))

# If a parallel cluster is provided, export necessary variables that will not change
# with iteration
if(!is.null(pClust)) {
  clusterExport(cl = pClust, varlist = c("Nblock","blocks","Npar","pSnooker","mZ",
                "Z","Fs","eps.mult","p.f.is.1","F1","f","eps.add","FUN"), envir = environment())
  clusterExport(cl = pClust, varlist = names(list(...)), envir = parent.env(environment()))
}
for (iter in 1:n.generation) {
  accepti = 0
  # Generate all necessary random numbers before running the chains
  # There is one sequence for each usage of runif or rnorm inside the calc_block function
  prob_snooker = runif(length(chainset)*Nblock)
  prob_acceptance = runif(length(chainset)*Nblock)
  prob_p.f.is.1 = runif(length(chainset)*Nblock)
  err_gamma_snooker = runif(length(chainset)*Nblock, min=1-eps.mult,max=1+eps.mult)
  err_gamma_par = runif(length(chainset)*Npar, min=1-eps.mult,max=1+eps.mult)
  err_eps.add = rnorm(length(chainset)*Npar, 0,1)
  rr_vec = vector("list",Nchain)
  for(i in chainset) {
    rr_vec[[i]] = vector("list",Nblock)
    for(j in 1:Nblock) {
      rr_vec[[i]][[j]] = sample(1:mZ, 3, replace = FALSE)
    }
  }
  # Use nested foreach if a parallel cluster is provided
  if(!is.null(pClust)) {
    # Export variables that will change in each iteration
    clusterExport(cl = pClust, varlist = c("X","logfitness_X","Naccept",
                  "prob_snooker","prob_acceptance","prob_p.f.is.1","err_gamma_snooker",
                  "err_gamma_par","err_eps.add","rr_vec"),
                 envir = environment())
    # Parallel for-loop. User is responsible for setting up the parallel cluster
    # and its management. This returns a list with the values of X, logfitness and Naccept for each chain
    out = foreach(i = chainset,.errorhandling = "stop") %dopar% {
        calc_block(X, logfitness_X, Naccept, Nblock, Npar, blocks,i,pSnooker,mZ,Z,Fs,eps.mult,p.f.is.1,F1,f,
                   eps.add,FUN, prob_snooker, prob_acceptance, prob_p.f.is.1,err_gamma_snooker,err_gamma_par,
                   err_eps.add, rr_vec,...)
    }
    for(i in 1:length(out)) {
      X[,i] = out[[i]][["X"]]
      logfitness_X[i] = out[[i]][["logfitness_X"]]
      Naccept[,i] = out[[i]][["Naccept"]]
    }
    # If no parallel cluster is provided, use nested for loops
  } else {
    for (i in chainset){
      out = calc_block(X, logfitness_X, Naccept, Nblock, Npar, blocks,i,pSnooker,mZ,Z,Fs,eps.mult,p.f.is.1,F1,f,
                      eps.add,FUN, prob_snooker, prob_acceptance, prob_p.f.is.1,err_gamma_snooker,err_gamma_par,
                      err_eps.add, rr_vec,...)
      X[,i] = out[["X"]]
      logfitness_X[i] = out[["logfitness_X"]]
      Naccept[,i] = out[["Naccept"]]
    }
  }
  if (!(iter%%n.thin) ){
    Z = cbind(Z,X)
    mZ = ncol(Z)
  }
} # n.generation
 if (!is.update) Z = Z[,-(1:(M0 + Nchain* floor(n.burnin/n.thin))),drop = FALSE]
 out = list(Draws=  Z , accept.prob.mat= Naccept/n.generation, X.final = X, logfitness.X.final = logfitness_X, Nchain=Nchain, n.generation= ncol(Z)/Nchain, demc_zs = TRUE)
 class(out) <- c("demc")
 invisible(out)
}


# Function that implements the loop over blocks
calc_block = function(X, logfitness_X, Naccept, Nblock, Npar, blocks,i,pSnooker,mZ,Z,Fs,eps.mult,p.f.is.1,F1,f,
                      eps.add,FUN, prob_snooker, prob_acceptance, prob_p.f.is.1,err_gamma_snooker,
                      err_gamma_par, err_eps.add, rr_vec,...) {
  # cumulative variable required for multiplicative error of gamma_par
  cdsub = 1
  for (iblock in 1:Nblock)  {
    # set of parameters to update
    parset = blocks[[iblock]]
    dsub = length(parset)  # dsubspace dimension of sampling
    # select to random different individuals from Z in rr, a 2-vector
    if ( prob_snooker[(i-1)*Nblock + iblock] < pSnooker ) { # prob_snooker[(i-1)*Nblock + iblock] = runif(1)
      # DE-Snooker update
      # if (Nchain >1) { z =  X[,sample(chainset[-i],1)]} else { # no real advantage and precludes parallel computing
      rr = rr_vec[[i]][[iblock]][1:3] #sample(1:mZ, 3, replace = FALSE)
      z = Z[,rr[3]]
      # }                                                  # no real advantage and precludes parallel computing
      x_z = X[ parset,i] - z[parset]
      D2 = max(sum(x_z*x_z),1.0e-300)
      #gamma_snooker = runif(1, min=Fs[iblock]-0.5,max=Fs[iblock]+0.5)
      gamma_snooker = Fs[iblock] * err_gamma_snooker[(i-1)*Nblock + iblock] # runif(1, min=1-eps.mult,max=1+eps.mult)
      #gamma_snooker =1.7
      projdiff = sum((Z[ parset,rr[1]] -Z[ parset,rr[2]]) *x_z)/D2  # inner_product of difference with x_z / squared norm x_z
      x_prop = X[,i]
      x_prop_sub = X[parset,i] + (gamma_snooker * projdiff) * x_z
      x_prop[parset] = x_prop_sub
      x_z = x_prop - z
      D2prop = max(sum(x_z*x_z), 1.0e-30)
      Npar12  =(dsub - 1)/2
      r_extra = Npar12 * (log(D2prop) - log(D2))   # Npar12 = extra term in logr for accept - reject ratio
    } else {
      # DE-parallel direction update
      if ( prob_p.f.is.1[(i-1)*Nblock + iblock] < p.f.is.1 ) { # prob_p.f.is.1[(i-1)*Nblock + iblock] = runif(1)
        gamma_par = F1 # to be able to jump between modes
      } else {
        gamma_par = (f[iblock]/sqrt(2*dsub)) * err_gamma_par[(i-1)*Npar + cdsub:dsub] # runif(dsub, min=1-eps.mult, max=1+eps.mult)   # multiplicative error to be applied to the difference
        # gamma_par = F2
      }
      rr = rr_vec[[i]][[iblock]][1:2] #sample(1:mZ, 2, replace = FALSE)

      if (eps.add ==0) {  # avoid generating normal random variates if possible
        x_prop_sub = X[parset,i] + gamma_par * (Z[parset,rr[1]]-Z[parset,rr[2]]) } else {
          x_prop_sub = X[parset,i] + gamma_par * (Z[parset,rr[1]]-Z[parset,rr[2]])  +  eps.add*err_eps.add[(i-1)*Npar + cdsub:dsub] # rnorm(dsub,0,1)
        }
      x_prop = X[,i]
      x_prop[parset] = x_prop_sub
      r_extra = 0
    }
    logfitness_x_prop = FUN(x_prop,  ...)
    #logfitness_x_prop = FUN(x_prop,  data = Data)
    logr =  logfitness_x_prop - logfitness_X[i]
    # print(c(logfitness_X[i], logfitness_x_prop ,logr,r_extra))
    if (!is.na(logr) & (logr + r_extra)> log(prob_acceptance[(i-1)*Nblock + iblock]) ){ #  prob_acceptance[(i-1)*Nblock + iblock]) = runif(1)
      Naccept[iblock,i] = Naccept[iblock,i]+1
      X[,i] = x_prop
      logfitness_X[i] = logfitness_x_prop
    }
    # update cumulative dsub
    cdsub = cdsub + dsub
  }
  return(list(X = X[,i], logfitness_X = logfitness_X[i], Naccept = Naccept[,i]))
}