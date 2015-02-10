summary.demc <- function (a, trans=NULL, keep.all=TRUE, Rupper.keep=FALSE) {
  # adapted from monitor from Andrew Gelman's monitor
  # a in this interface is k x (m x k) array [ from DE.MC] instead of n x m x k
  # with k the number of parameters, and m the number of chains
  #
  # If keep.all=TRUE:  a is a n x m x k array:
  #   m sequences of length n, k variables measured
  # If keep.all=FALSE:  a is a 2n x m x k array (first half will be discarded)
  #
  # trans is a vector of length k:  "" if no transformation, or "log" or "logit"
  # (If trans is not defined, it will be set to "log" for parameters that
  # are all-positive and 0 otherwise.)
  #
  # If Rupper.keep=TRUE:  keep Rupper.  (Otherwise don't display it.)
  if ("demc"%in%class(a)) {
    
    m = a$Nchain
    accept = a$accept.prob.mat
    a = a$Draws
    } else 
    stop("argument to summary.demc must be of class demc")
  output <- NULL
  #  nparams <- ifelse (length(dim(a))<3, 1, dim(a)[length(dim(a))])
  nparams <- nrow(a)
  #  if (length(dim(a))==2) a <- array (a, c(dim(a),1))
  if (!keep.all){
    half <- floor(ncol(a)/2)
    a <- a[, (half+1): (2*half),drop=FALSE]
  }
if (m>1){
  if (is.null(trans))
    trans <- ifelse ((apply (a<=0, 1, sum))==0, "log", "")
  for (i in 1:nparams){
    #    if (trans[i]=="log") conv.p <- conv.par.log(a[,,i])
    if (trans[i]=="log") conv.p <- conv.par.log(matrix(a[i,],nrow= m))
    else if (trans[i]=="logit") conv.p <- conv.par.logit(matrix(a[i,],nrow= m))
    else conv.p <- conv.par(matrix(a[i,],nrow= m))
    output <- rbind (output, c(mean(a[i,]), sqrt(var(as.vector(a[i,]))),
                               conv.p$quantiles, conv.p$confshrink, round.sci(conv.p$n.eff,2)))
  }

} else {
  for (i in 1:nparams){
    quant = quantile(a[i,], probs=c(.025,.25,.5,.75,.975))
    output <- rbind (output, c(mean(a[i,]), sqrt(var(as.vector(a[i,,drop=FALSE]))),
                               quant,  c(rep(NA,3))))
  }
}
dimnames(output) <- list(dimnames(a)[[1]],c("mean","sd",
                                            "2.5%","25%","50%","75%","97.5%","Rhat","Rupper","n.eff"))
if (!Rupper.keep)                     # discard Rupper (nobody ever uses it)
  output = output[,-(ncol(output)-1)]
list(parameters = output, accept.prob = accept, mean.accept = mean(accept))
}
