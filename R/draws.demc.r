#' @title Extracts the draws from a demc object
#'
#'@param object \code{demc} object
#'@param keep.all logical or fraction of draws to keep. If fraction, the one-complement is discarded as burnin.
#'If logical the first half of the draws are discarded. Default 0.9.
#'@param thin integer (default 1), keep every thin-th draw
#'@return a dataframe with parameters as variables and as rows the kept samples
#'@export draws.demc

draws.demc <- function(object, keep.all=0.9, thin = 1){
  if ("demc"%in%class(object)) {
    a = object$Draws
    m = object$Nchain
    } else 
    stop("argument to draws.demc must be of class demc")
  if (is.logical(keep.all)) keep.all = ifelse(keep.all,1,0.5)
  if (keep.all < 1){
    n.gen = ncol(a)/m
    start = m * (floor((1-keep.all)*n.gen)+1)
    a <- a[, (start+1) : ncol(a),drop=FALSE]
    # thin: but in batches of m = Nchain  
    indices = c(t(sapply(seq_len(m), function(i){seq(from=i, to = ncol(a), by = thin*m)})))
  } else indices = seq_len(ncol(a))
  as.data.frame(t(a[,indices,drop=FALSE]))
}