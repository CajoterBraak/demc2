demc_zs <- function(Nchain = 3, Z, FUN, X= matrix(Z[,1:Nchain,drop=FALSE], ncol = Nchain),
                     blocks = list(seq_len(nrow(X))), f = 2.38, pSnooker= 0.1, pGamma1 = 0.1, n.generation = 10, 
 n.thin = 1, n.burnin = 0, eps.mult =0.2,eps.add = 0, ...){
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
M0 = mZ = ncol(Z)
#Nchain = ncol(X)
#Npar = nrow(X)
#Npar12  =(Npar - 1)/2  # factor for Metropolis ratio DE Snooker update
F1 = 1.0
accept = rep(NA,n.generation)
chainset = seq_len(Nchain)
rr = NULL
r_extra = 0
logfitness_X = apply (X, 2, FUN, ...)
Nblock = length(blocks)
f = matrix(f,nrow=1,ncol=Nblock)
Fs = f/sqrt(2)
Naccept = matrix(0,nrow = Nblock, ncol = Nchain, dim=list(blocks=paste("block", seq_len(Nblock)), chains=paste("chain", chainset) ))
for (iter in 1:n.generation) {
   accepti = 0
     for (i in chainset){
       for (iblock in 1:Nblock)  {
        # set of parameters to update
        parset = blocks[[iblock]]
        dsub = length(parset)  # dsubspace dimension of sampling 
       # select to random different individuals from Z in rr, a 2-vector
       if ( runif(1)< pSnooker ) {  
       # DE-Snooker update
         # if (Nchain >1) { z =  X[,sample(chainset[-i],1)]} else { # no real advantage and precludes parallel computing
          rr = sample(1:mZ, 3, replace = FALSE)
          z = Z[,rr[3]]
         # }                                                  # no real advantage and precludes parallel computing
          x_z = X[ parset,i] - z[parset]
          D2 = max(sum(x_z*x_z),1.0e-300)
          gamma_snooker = runif(1, min=Fs[iblock]-0.5,max=Fs[iblock]+0.5)
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
         if ( runif(1)< pGamma1 ) { gamma_par = F1 # to be able to jump between modes
          } else {
          gamma_par = (f[iblock]/sqrt(2*dsub)) * runif(dsub, min=1-eps.mult, max=1+eps.mult)    # multiplicative error to be applied to the difference
           # gamma_par = F2 
         }
         rr = sample(1:mZ, 2, replace = FALSE)
     
         if (eps.add ==0) {  # avoid generating normal random variates if possible
           x_prop_sub = X[parset,i] + gamma_par * (Z[parset,rr[1]]-Z[parset,rr[2]]) } else {
           x_prop_sub = X[parset,i] + gamma_par * (Z[parset,rr[1]]-Z[parset,rr[2]])  +  eps.add*rnorm(dsub,0,1)
         }
         x_prop = X[,i]
         x_prop[parset] = x_prop_sub
         r_extra = 0
        }
       logfitness_x_prop = FUN(x_prop,  ...)
       logr =  logfitness_x_prop - logfitness_X[i]
      # print(c(logfitness_X[i], logfitness_x_prop ,logr,r_extra))
       if (!is.na(logr) & (logr + r_extra)> log(runif(1)) ){
          Naccept[iblock,i] = Naccept[iblock,i]+1
          X[,i] = x_prop
          logfitness_X[i] = logfitness_x_prop
       }
    } # iblock         
  } # i loop
  if (!(iter%%n.thin) ){   
   Z = cbind(Z,X)
   mZ = ncol(Z)
   }

} # n.generation
 out = list(Draws= Z[,-(1:(M0 + Nchain* floor(n.burnin/n.thin))),drop = FALSE] , accept.prob.mat= Naccept/n.generation, X.final = X, logfitness.X.final = logfitness_X, Nchain=Nchain)
 class(out) <- c("demc")
 invisible(out)
}
