
# Optimal bandwidth selection with known function for univariate local polynomial estimator 

h.optim = function(func,n, p, p.eval, N, h.seq, m = 2){
  
  x      = (1:p - 0.5)/p
  x.eval = (1:p.eval - 0.5)/p.eval
  
  sup.error = future_sapply(1:length(h.seq), function(i){
    weights = locPolWeights(x = x, bw = h.seq[i], deg = m, xeval = x.eval, kernel=EpaK)$locWeig
    
    max.errors = sapply(1:N, function(k){
      Z_mean = apply(sim.d.OU(n, t = x,rho_B = 0.5,tau = 0,sigma = 4),2,mean)
      e_mean = rnorm(p, mean=0, sd = 0.1/sqrt(n))
      sample = func(x) + e_mean + Z_mean
      max( abs( as.vector( weights %*% (sample)) - func(x.eval)))
    })
    
    mean(max.errors)
  }, future.seed = T)
  
  return(sup.error)
}


# Optimal bandwidth selection with known function for bivariate local polynomial estimator 

h.optim.cov = function(N, n, p, p.eval, h.seq, max.lag = 0, m = 1, w.parallel = T, func){
  
  x.design = (1:p - 0.5)/p
  
  if(max.lag == 0){help = function(h){
    w_h = local.polynomial.weights(p, h, p.eval = p.eval, m = m, parallel = w.parallel, parallel.environment = F, grid.type = "less")
    max_diff = numeric(1)
    
    true.kernel = cov.d.OU(p.eval = p.eval, rho_B = 0.5, tau = 0, sigma = 4, m_lag = 0)
    
    future_replicate(N, {
      Y        = t(func(x.design) + t(sim.d.OU(n, t = x.design,rho_B = 0.5,tau = 0,sigma = 4)) + matrix(rnorm(length(x.design)*n, 0, 0.1),length(x.design), n))
      pred     = eval_weights(w_h, observation.transformation(Y, grid.type = "less"))
      max_diff = max(abs((true.kernel - pred)[!as.logical(diag(p.eval))]))
      
      max_diff
    }, future.seed = T)
  }}else{
    help = function(h){
      w_h = local.polynomial.weights(p, h, p.eval = p.eval, m = m, parallel = w.parallel, parallel.environment = F, grid.type = "lesseq")
      max_diff = matrix(numeric(max.lag), nrow = max.lag)
      
      future_replicate(N, {
        Y             = t(func(x.design) + t(sim.d.OU(n, t = x.design,rho_B = 0.5,tau = 0,sigma = 4)) + matrix(rnorm(length(x.design)*n, 0, 0.1),length(x.design), n))
        
        for(j in 1:max.lag){
          pred         = eval.weights(w_h, observation.transformation(Y, grid.type = "lesseq", lag = j), lag = j)
          true.kernel = cov.d.OU(p.eval = p.eval, rho_B = 0.5, tau = 0, sigma = 4, m_lag = j)
          max_diff[j,1] = max(abs((true.kernel + t(true.kernel) - pred - t(pred))))
        }
        
        max_diff
      }, future.seed = T)
    }
    
    mean_sup = sapply(h.seq, help)
    
    M        = matrix(mean_sup,nrow = max.lag)
    M.row    = nrow(M)
    M.col    = ncol(M)
    hblock   = M.col / N
    arr      = array(M, dim = c(M.row, N, hblock))
    res.mean = apply(arr, c(1, 3), mean)
    
    bw.optim = h.seq[apply(res.mean,1,which.min)]
    names(bw.optim) = paste0("lag ", 1:max.lag)
    return(bw.optim) 
  }
  
  mean_sup = sapply(h.seq, help)
  
  h.seq[which.min(apply(mean_sup, 2, mean))]
}



###########################################################
##### hv-block cross validation in a K-fold framework #####
###########################################################

# K-Fold hv-block cross validation for univariate local polynomial estimator 

k.fold.hc.cv = function(N, n, p, h.seq, K = 5,gap = 0, m = 2, func, diff = F,...){
  
  grp = rep(1:K, each = ceiling(n/K))
  x.design = (1:p - 0.5)/p
  
  if(diff == F){help = function(h){
    w_h = locPolWeights(x = x.design, bw = h, deg = m, xeval= x.design, kernel=EpaK)$locWeig
    max_diff = numeric(K)
    
    future_replicate(N, {
      Y = t(func(x.design) + t(sim.d.OU(n, t = x.design,rho_B = 0.5,tau = 0,sigma = 4)) + matrix(rnorm(length(x.design)*n, 0, 0.1),length(x.design), n))
      for(kk in 1:K){
        test     = setdiff(1:n,setdiff(1:n,which(grp == kk)))
        l.gap    = max(min(test) - gap, 1) : (min(test) - 1)
        l.gap    = l.gap[l.gap >= 1 & l.gap <= n]
        r.gap    = (max(test) + 1) : min(max(test) + gap, n)
        r.gap    = r.gap[r.gap >= 1 & r.gap <= n]
        full.gap = c(l.gap, r.gap)
        train    = setdiff(1:n, c(test, full.gap))
        test_grp     = apply(Y[test,],2,mean)
        pred         = w_h %*% apply(Y[train,],2,mean)
        max_diff[kk] = max(abs(test_grp - pred))
      }
      mean(max_diff)
    }, future.seed = T)
  }
  
  mean_sup = sapply(h.seq, help)
  
  return(h.seq[apply(mean_sup, 1, which.min)])}
  else{help = function(h,...){
    
    arg <- list(...)
    nd = arg[[1]]
    pd = arg[[2]]
    h.dense = arg[[3]]
    func2   = arg[[4]]
    
    w.dense = locPolWeights(x = (1:pd - 0.5)/pd, bw = h.dense, deg = m, xeval= x.design, kernel=EpaK)$locWeig
    w_h = locPolWeights(x = x.design, bw = h, deg = m, xeval= x.design, kernel=EpaK)$locWeig
    max_diff = numeric(K)
    
    future_replicate(N, {
      Y  = t(func(x.design)  + t(sim.d.OU(n, t = x.design,rho_B = 0.5,tau = 0,sigma = 4)) + matrix(rnorm(length(x.design)*n, 0, 0.1),length(x.design), n))
      Yd = t(func2((1:pd - 0.5)/pd) + t(sim.d.OU(nd, t = (1:pd - 0.5)/pd,rho_B = 0.5,tau = 0,sigma = 4)) + matrix(rnorm(length((1:pd - 0.5)/pd)*nd, 0, 0.1),length((1:pd - 0.5)/pd), nd))
      res = t(t(Y) - c(w.dense %*% apply(Yd,2,mean)))
      
      for(kk in 1:K){
        test     = setdiff(1:n,setdiff(1:n,which(grp == kk)))
        l.gap    = max(min(test) - gap, 1) : (min(test) - 1)
        l.gap    = l.gap[l.gap >= 1 & l.gap <= n]
        r.gap    = (max(test) + 1) : min(max(test) + gap, n)
        r.gap    = r.gap[r.gap >= 1 & r.gap <= n]
        full.gap = c(l.gap, r.gap)
        train    = setdiff(1:n, c(test, full.gap))
        
        test_grp     = apply(res[test,],2,mean)
        pred         = w_h %*% apply(res[train,],2,mean)
        max_diff[kk] = max(abs(test_grp - pred))
      }
      mean(max_diff)
    }, future.seed = T)
  }
  
  mean_sup = sapply(h.seq, function(x){help(x,...)})
  
  return(h.seq[apply(mean_sup, 1, which.min)])}
}



#K-Fold hv-block cross validation for bivariate local polynomial estimator 

k.fold.hc.cv.cov = function(N, n, p, h.seq, K = 5,gap = 0, max.lag = 0, m = 1, w.parallel = T, func){
  
  grp = rep(1:K, each = ceiling(n/K))
  x.design = (1:p - 0.5)/p
  
  if(max.lag == 0){help = function(h){
    w_h = local.polynomial.weights(p, h, p.eval = p, m = m, parallel = w.parallel, parallel.environment = F, grid.type = "less")
    max_diff = numeric(K)
    
    future_replicate(N, {
      Y = t(func(x.design) + t(sim.d.OU(n, t = x.design,rho_B = 0.5,tau = 0,sigma = 4)) + matrix(rnorm(length(x.design)*n, 0, 0.1),length(x.design), n))
      for(kk in 1:K){
        test     = setdiff(1:n,setdiff(1:n,which(grp == kk)))
        l.gap    = max(min(test) - gap, 1) : (min(test) - 1)
        l.gap    = l.gap[l.gap >= 1 & l.gap <= n]
        r.gap    = (max(test) + 1) : min(max(test) + gap, n)
        r.gap    = r.gap[r.gap >= 1 & r.gap <= n]
        full.gap = c(l.gap, r.gap)
        train    = setdiff(1:n, c(test, full.gap))
        
        test_grp     = matrix(observation.transformation(Y[test,], grid.type = "full"), p, p)
        train_grp    = observation.transformation(Y[train,], grid.type = "less")
        pred         = eval_weights(w_h, train_grp)
        max_diff[kk] = max(abs((test_grp - pred)[!as.logical(diag(p))]))
      }
      mean(max_diff)
    }, future.seed = T)
  }}
  else{help = function(h){
    w_h = local.polynomial.weights(p, h, p.eval = p, m = m, parallel = w.parallel, parallel.environment = F, grid.type = "lesseq")
    max_diff = matrix(numeric(K*max.lag), nrow = max.lag)
    
    future_replicate(N, {
      Y = t(func(x.design) + t(sim.d.OU(n, t = x.design,rho_B = 0.5,tau = 0,sigma = 4)) + matrix(rnorm(length(x.design)*n, 0, 0.1),length(x.design), n))
      for(kk in 1:K){
        test     = setdiff(1:n,setdiff(1:n,which(grp == kk)))
        l.gap    = max(min(test) - gap, 1) : (min(test) - 1)
        l.gap    = l.gap[l.gap >= 1 & l.gap <= n]
        r.gap    = (max(test) + 1) : min(max(test) + gap, n)
        r.gap    = r.gap[r.gap >= 1 & r.gap <= n]
        full.gap = c(l.gap, r.gap)
        train    = setdiff(1:n, c(test, full.gap))
        
        for(j in 1:max.lag){
          test_grp       = matrix(observation.transformation(Y[test,], grid.type = "full", lag = j), p, p)
          train_grp      = observation.transformation(Y[train,], grid.type = "lesseq", lag = j)
          pred           = eval.weights(w_h, train_grp, lag = j)
          max_diff[j,kk] = max(abs((test_grp + t(test_grp) - pred - t(pred))))
        }
      }
      apply(max_diff,1,mean)
    }, future.seed = T)
  }
  
  mean_sup = sapply(h.seq, help)
  print("done")
  bw.list = matrix(h.seq[apply(mean_sup, 1, which.min)],nrow = max.lag)
  row.names(bw.list) = paste0("lag ", 1:max.lag)
  return(bw.list)
  
  }
  
  mean_sup = sapply(h.seq, help)
  print("done")
  h.seq[apply(mean_sup, 1, which.min)]
}
