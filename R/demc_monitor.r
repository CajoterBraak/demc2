# help functions for monitor.demc
conv.par <- function (x) {
  alpha <- .05                     # 95% intervals
  m <- ncol(x)
  n <- nrow(x)

# We compute the following statistics:
#
#  xdot:  vector of sequence means
#  s2:  vector of sequence sample variances (dividing by n-1)
#  W = mean(s2):  within MS
#  B = n*var(xdot):  between MS.
#  muhat = mean(xdot):  grand mean; unbiased under strong stationarity
#  varW = var(s2)/m:  estimated sampling var of W
#  varB = B^2 * 2/(m+1):  estimated sampling var of B
#  covWB = (n/m)*(cov_monitor(s2,xdot^2) - 2*muhat*cov_monitor(s^2,xdot)):
#                                               estimated sampling cov(W,B)
#  sig2hat = ((n-1)/n))*W + (1/n)*B:  estimate of sig2; unbiased under
#                                               strong stationarity
#  quantiles:  emipirical quantiles from last half of simulated sequences

  xdot <- apply(x,2,mean)
  s2 <- apply(x,2,var)
  W <- mean(s2)
  B <- n*var(xdot)
  muhat <- mean(xdot)
  varW <- var(s2)/m
  varB <- B^2 * 2/(m-1)
  covWB <- (n/m)*(cov_monitor(s2,xdot^2) - 2*muhat*cov_monitor(s2,xdot))
  sig2hat <- ((n-1)*W + B)/n
  quantiles <- quantile (as.vector(x), probs=c(.025,.25,.5,.75,.975))

  if (W > 1.e-8) {            # non-degenerate case

# Posterior interval post.range combines all uncertainties
# in a t interval with center muhat, scale sqrt(postvar),
# and postvar.df degrees of freedom.
#
#       postvar = sig2hat + B/(mn):  variance for the posterior interval
#                               The B/(mn) term is there because of the
#                               sampling variance of muhat.
#       varpostvar:  estimated sampling variance of postvar

    postvar <- sig2hat + B/(m*n)
    varpostvar <- max (0,
      (((n-1)^2)*varW + (1+1/m)^2*varB + 2*(n-1)*(1+1/m)*covWB)/n^2)
    post.df <- min (chisqdf (postvar, varpostvar), 1000)
    post.range <- muhat + sqrt(postvar) * qt(1-alpha/2, post.df)*c(-1,0,1)

# Estimated potential scale reduction (that would be achieved by
# continuing simulations forever) has two components:  an estimate and
# an approx. 97.5% upper bound.
#
# confshrink = sqrt(postvar/W),
#     multiplied by sqrt(df/(df-2)) as an adjustment for the
###      CHANGED TO sqrt((df+3)/(df+1))
#     width of the t-interval with df degrees of freedom.
#
# postvar/W = (n-1)/n + (1+1/m)(1/n)(B/W); we approximate the sampling dist.
# of (B/W) by an F distribution, with degrees of freedom estimated
# from the approximate chi-squared sampling dists for B and W.  (The
# F approximation assumes that the sampling dists of B and W are independent;
# if they are positively correlated, the approximation is conservative.)

    varlo.df <- chisqdf (W, varW)
    confshrink.range <- sqrt (c(postvar/W,
      (n-1)/n + (1+1/m)*(1/n)*(B/W) * qf(.975, m-1, varlo.df)) *
      (post.df+3)/(post.df+1))

# Calculate effective sample size:  m*n*min(sigma.hat^2/B,1)
# This is a crude measure of sample size because it relies on the between
# variance, B, which can only be estimated with m degrees of freedom.
    
    n.eff <- m*n*min(sig2hat/B,1)
    list(post=post.range, quantiles=quantiles, confshrink=confshrink.range,
         n.eff=n.eff)
  }
  else {      # degenerate case:  all entries in "data matrix" are identical
    list (post=muhat*c(1,1,1), quantiles=quantiles, confshrink=c(1,1),
          n.eff=1)
  }
}


cov_monitor <- function (a,b) {
  m <- length(a)
  ((mean ((a-mean(a))*(b-mean(b)))) * m)/(m-1)}
logit <- function (x) {log(x/(1-x))}
invlogit <- function (x) {1/(1+exp(-x))}

chisqdf <- function (A, varA) {
# Chi-squared degrees of freedom estimated by method of moments
# (Assume A has a gamma sampling distribution and varA is an unbiased
# estimate of the sampling variance of A.)
  2*(A^2/varA)}

round.sci <- function (a, digits){
# Round to the specified number of significant digits
  round (a, min(0,digits-1-floor(log10(a))))
}
conv.par.log <- function (r) {
  conv.p <- conv.par(log(r))
  list (post=exp(conv.p$post), quantiles=exp(conv.p$quantiles),
    confshrink=conv.p$confshrink, n.eff=conv.p$n.eff)}

conv.par.logit <- function (r) {
  conv.p <- conv.par(logit(r))
  list (post=invlogit(conv.p$post), quantiles=invlogit(conv.p$quantiles),
    confshrink=conv.p$confshrink, n.eff=conv.p$n.eff)}

