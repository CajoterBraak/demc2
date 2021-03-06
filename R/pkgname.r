#'Adaptive MCMC on real parameter spaces via Differential Evolution Markov Chain 
#'
#'demc2 implements both the origal Markovian \code{demc} version of ter Braak (2006) and
#' the egodic \code{demc_zc} version of ter Braak & Vrugt (2008) using learning from the past and the snooker update.
#' Convergence can be checked with coda.
#' The demc package allows restarts.
#' 
#' The functions that you're likely need from \pkg{demc} are
#' \code{\link{demc}} and  \code{\link{demc_zs}} for generating samples and
#' \code{\link{as.coda}} to be able to use package \code{coda} for convergence checks

#'@references ter Braak, C. J. F. (2006). A Markov Chain Monte Carlo version of the genetic algorithm Differential Evolution: easy Bayesian computing for real parameter spaces. Statistics and Computing, 16, 239-249.\url{http://edepot.wur.nl/26336}
#' ter  Braak C. J. F., and Vrugt J. A. (2008). Differential Evolution Markov Chain with snooker updater and fewer chains. Statistics and Computing, 18 (4), 435-446. \url{http://dx.doi.org/10.1007/s11222-008-9104-9}
#' 
#'@docType package
#'@name demc2
NULL