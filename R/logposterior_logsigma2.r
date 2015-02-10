logposterior_logsigma2 <- function(theta, data){ 
  theta_original = theta
  theta_original[3] =  exp(theta[3]) #    sigma2 
  logposterior_sigma2(theta_original, data = data) + theta[3]   # theta[3] = log(sigma2)
  # the second term stems from the Jacobian, slide LL122, section 2.9 of the LL book
}