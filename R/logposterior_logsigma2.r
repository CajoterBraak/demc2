logposterior_logsigma2 <- function(theta, data, v1 = 10, v2 = 5, nu = 2, tau2 = 3){ 
  theta_original = theta
  theta_original[3] =  exp(theta[3]) #    sigma2 
  logposterior_sigma2(theta_original, data = data, v1,v2,nu,tau2) + theta[3]   # theta[3] = log(sigma2)
  # the second term stems from the Jacobian, slide LL122, section 2.9 of the LL book
}