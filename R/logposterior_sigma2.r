logposterior_sigma2 <- function(theta, data){
    # make data variates local 
  y = data$y; x1= data$X1 ; x2 = data$X2
  # make the parameters available
  b1 = theta[1]
  b2 = theta[2]
  sigma2 = theta[3] # sigma squared
  if (sigma2 >0) {
    # prior
    logprior = dnorm(b1, 0, sqrt(v1), log= TRUE)+dnorm(b2,0,sqrt(v2), log= TRUE)+geoR::dinvchisq(sigma2, nu, tau2, log= TRUE) 
    # likelihood 
    loglikelihood = sum(dnorm(y, mean = b1*x1+b2*x2, sd = sqrt(sigma2),log= TRUE))
    logpost = logprior + loglikelihood  
  } else logpost = -Inf # because logprior = -Inf
  logpost
} # end of function logposterior
