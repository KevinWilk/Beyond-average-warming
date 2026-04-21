######################################################################
####################### functions ####################################
######################################################################

mu = function(x){
  x = 2*x -1
  return(3*sin(1.5*pi*(x))*exp(-2*abs(x)))
}

delta = function(x,constant = T){
  if(constant == T) return(rep(2, length(x)))
  else{ x = 2*x -1
  return(-sin(pi*(x))*exp(-2*abs(x))+2)}
}

mu_d = function(x,...){
  mu(x)-delta(x,...)
}


##########################################################
################ Dependent Ornstein-Uhlenbeck  ###########
##########################################################

sim.d.OU = function(n,t = seq(0, 1, len = 201),theta = 1, sigma = 1, rho = 0 ,tau = 1, rho_B = 0) {
  
  stopifnot(abs(rho) < 1, abs(rho_B) < 1)
  
  dt = diff(t)[1]
  m  = length(t)
  
  mu = numeric(n)
  mu[1] = rnorm(1, mean = 0, sd = tau / sqrt(1 - rho^2))  
  for (k in 2:n) {
    mu[k] = rho * mu[k - 1] + rnorm(1, 0, tau)
  }
  
  d = 0:(n-1)
  Sigma = toeplitz(rho_B^d)
  L = t(chol(Sigma))  
  
  X = matrix(NA_real_, nrow = n, ncol = m)
  X[, 1] = rnorm(n, mean = mu, sd = sigma / sqrt(2 * theta))
  
  for (k in 2:m) {
    phi = exp(-theta * dt)
    sd_incr = sqrt((1 - phi^2) * sigma^2 / (2 * theta))
    Z = matrix(rnorm(n), ncol = 1)
    Z_cor = as.vector(L %*% Z)
    X[, k] = mu + (X[, k - 1] - mu) * phi + sd_incr * Z_cor
  }
  
  colnames(X) = sprintf("%.3f", t)
  return(X)
}


cov.d.OU = function(p.eval,theta = 1, sigma = 1, rho = 0, tau = 1, rho_B  = 0, m_lag = 0) {
  
  stopifnot(abs(rho) < 1, abs(rho_B) < 1)
  
  
  t = seq(0, 1, length.out = p.eval)
  
  Tdiff = abs(outer(t, t, `-`))
  Tsum  = outer(t, t, `+`)
  K1 = exp(-theta * Tdiff)          
  K2 = exp(-theta * Tsum)          
  
  Jp   = matrix(1, p.eval, p.eval)
  
  AR_part = (tau^2 * rho^m_lag/ (1 - rho^2)) * Jp
  OU_part = (sigma^2 / (2 * theta)) * (rho_B^m_lag * (K1 - K2) + if (m_lag == 0) K2 else 0 )
  
  Sigma = AR_part + OU_part
  
  return(Sigma)
}

LR.cov.d.OU = function(m_lag,...){
  results = future_lapply(0:m_lag,function(i){(if (i == 0) 1 else 2 )*cov.d.OU(m_lag = i,...)},future.seed = TRUE)
  return(Reduce(`+`, results))
}

