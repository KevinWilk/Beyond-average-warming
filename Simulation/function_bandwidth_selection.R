
# Optimal bandwidth selection with known function for univariate local polynomial estimator 

h.optim = function(func,n, p, p.eval, N, h.seq, m = 2, norm = "sup-norm"){
  
  x      = (1:p - 0.5)/p
  x.eval = (1:p.eval - 0.5)/p.eval
  
  w.h.list = lapply(h.seq, function(h){locPolWeights(x = x, bw = h, deg = m, xeval = x.eval, kernel = EpaK)$locWeig})  
  true.val = func(x.eval) 
  
  err = future_sapply(1:N, function(k){
    Z_mean = apply(sim.d.OU(n, t = x,rho_B = 0.5,tau = 0,sigma = 4),2,mean)
    e_mean = rnorm(p, mean=0, sd = 0.1/sqrt(n))
    sample = func(x) + e_mean + Z_mean
    
    sapply(w.h.list, function(w_h){ if(norm == "sup-norm"){ max(abs(drop(w_h %*% sample) - true.val)  )      }
                                    if(norm == "L2-norm" ){ sum(   (drop(w_h %*% sample) - true.val)^2)^(1/2)}   } ) 
  }, future.seed = TRUE)
  
  sup.error = rowMeans(err)
  
  return( list(h.min = h.seq[which.min(sup.error)], h.list = data.frame(h = h.seq, mean.sup.error = sup.error)) )
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

cv.sim = function(N, n, p, h.seq, K = 5, gap = 0, m = 2, R = 5, func, diff = F, norm = "sup-norm",...){
  
  grp = rep(1:K, each = ceiling(n/K))
  x.design = (1:p - 0.5)/p
  
  w.h.list = lapply(h.seq, function(h){locPolWeights(x = x.design, bw = h, deg = m, xeval = x.design, kernel = EpaK)$locWeig})
  
  
  if(diff == T){
    
    arg     = list(...)
    nd      = arg[[1]]
    pd      = arg[[2]]
    h.dense = arg[[3]]
    func2   = arg[[4]]
    
    w.dense = locPolWeights(x = (1:pd - 0.5)/pd, bw = h.dense, deg = m, xeval= x.design, kernel=EpaK)$locWeig
    
  }
  
  ## R-shifted K-fold partitions
  L = ceiling(n/K)
  make.folds = function(shift){
    grp = ((0:(n-1) + shift) %/% L) %% K + 1                       
    lapply(1:K, function(kk){
      test = which(grp == kk)
      excl = unique(unlist(lapply(test, function(i) max(1, i - gap):min(n, i + gap))))
      list(test = test, train = setdiff(1:n, excl))
    })
  }
  shifts = round(seq(0, L, length.out = R + 1))[1:R]                
  folds  = do.call(c, lapply(shifts, make.folds))                   
  
  mean_sup = future_sapply(1:N, function(r){
    
    Y = t(func(x.design)  + t(sim.d.OU(n, t = x.design,rho_B = 0.5,tau = 0,sigma = 4)) + matrix(rnorm(length(x.design)*n, 0, 0.1),length(x.design), n))
    
    if(diff == T){
      
      Yd = t(func2((1:pd - 0.5)/pd) + t(sim.d.OU(nd, t = (1:pd - 0.5)/pd,rho_B = 0.5,tau = 0,sigma = 4)) + matrix(rnorm(length((1:pd - 0.5)/pd)*nd, 0, 0.1),length((1:pd - 0.5)/pd), nd)) 
      Y = t(t(Y) - c(w.dense %*% colMeans(Yd)))
    }
    
    test.Kfold   = sapply(folds, function(f) colMeans(Y[f$test,  , drop = FALSE]))  
    train.Kfold  = sapply(folds, function(f) colMeans(Y[f$train, , drop = FALSE]))   
    
    if(norm == "sup-norm"){ sapply(w.h.list, function(w_h){ mean(apply(abs(test.Kfold - as.matrix(w_h %*% train.Kfold))  , 2, max))} ) }
    if(norm == "L2-norm"){  sapply(w.h.list, function(w_h){ mean(apply(   (test.Kfold - as.matrix(w_h %*% train.Kfold))^2, 2, sum)^(1/2))} ) }
    
  }, future.seed = TRUE)                                                         
  
  mean_sup = t(mean_sup)                                                             
  h.seq[apply(mean_sup, 1, which.min)]
  
}



# for one sample

kfold.cv = function(Y, h.seq, K = 5, gap = 0, m = 2, R = 5, parallel = F, norm = "sup-norm"){
  
  n = dim(Y)[1]
  p = dim(Y)[2]
  
  x.design = (1:p - 0.5)/p
  
  L = ceiling(n/K)
  
  make.folds = function(shift){
    grp = ((0:(n-1) + shift) %/% L) %% K + 1
    lapply(1:K, function(kk){
      test = which(grp == kk)
      excl = unique(unlist(lapply(test, function(i) max(1, i - gap):min(n, i + gap))))
      list(test = test, train = setdiff(1:n, excl))
    })
  }
  
  shifts = round(seq(0, L, length.out = R + 1))[1:R]
  folds  = do.call(c, lapply(shifts, make.folds))
  
  test.Kfold   = sapply(folds, function(f) colMeans(Y[f$test,  , drop = FALSE]))  
  train.Kfold  = sapply(folds, function(f) colMeans(Y[f$train, , drop = FALSE]))   
  
  help = function(h){
    w_h = locPolWeights(x = x.design, bw = h, deg = m, xeval = x.design, kernel = EpaK)$locWeig
    if(norm == "sup-norm"){ return(mean(apply(abs(test.Kfold - w_h %*% train.Kfold)  , 2, max)))}
    if(norm == "L2-norm" ){ return(mean(apply(   (test.Kfold - w_h %*% train.Kfold)^2, 2, sum)))}   
  }                           
  
  mean_sup = if(parallel){ future_sapply(h.seq, help)}else{ sapply(h.seq, help)}
  
  return(h.seq[which.min(mean_sup)])
  
}

                        
#K-Fold hv-block cross validation for bivariate local polynomial estimator
                        
cv.sim.cov = function(N, n, p, h.seq, K = 5, gap = 0, R = 5, max.lag = 0, m = 1, w.parallel = T, func){
  
  x.design = (1:p - 0.5)/p
  
  if(max.lag == 0){ grid.type = "less" }else{ grid.type = "lesseq" }
  w.h.list = lapply(h.seq, function(h){ local.polynomial.weights(p, h, p.eval = p, m = m, parallel = w.parallel, parallel.environment = F, grid.type = grid.type) })
  
  
  ## R-shifted K-fold partitions
  L = ceiling(n/K)
  make.folds = function(shift){
    grp = ((0:(n-1) + shift) %/% L) %% K + 1                       
    lapply(1:K, function(kk){
      test = which(grp == kk)
      excl = unique(unlist(lapply(test, function(i) max(1, i - gap):min(n, i + gap))))
      list(test = test, train = setdiff(1:n, excl))
    })
  }
  shifts = round(seq(0, L, length.out = R + 1))[1:R]                
  folds  = do.call(c, lapply(shifts, make.folds))                   
  
  if(max.lag == 0){
    
    mean_sup = future_sapply(1:N, function(r){
    
    Y = t(func(x.design) + t(sim.d.OU(n, t = x.design,rho_B = 0.5,tau = 0,sigma = 4)) + matrix(rnorm(length(x.design)*n, 0, 0.1),length(x.design), n))
    
    test.Kfold   = lapply(folds, function(f) matrix(observation.transformation(Y[f$test,  , drop = FALSE], grid.type = "full"), p, p) )  
    train.Kfold  = lapply(folds, function(f)        observation.transformation(Y[f$train, , drop = FALSE], grid.type = grid.type)     )   
    
    sapply(w.h.list, function(w_h){ mean(sapply(1:(K*R),function(k){ max(abs((test.Kfold[[k]] - eval_weights(w_h, train.Kfold[[k]]))[!as.logical(diag(p))]))} ) )} ) 
    
    }, future.seed = TRUE)   
  
    mean_sup = t(mean_sup)
  
    return(h.seq[apply(mean_sup, 1, which.min)])
  
  }else{ 
  
    bw.list = future_sapply(1:N, function(r){
    
      Y = t(func(x.design) + t(sim.d.OU(n, t = x.design,rho_B = 0.5,tau = 0,sigma = 4)) + matrix(rnorm(length(x.design)*n, 0, 0.1),length(x.design), n))
    
      sapply(1:max.lag, function(j){
        
        test.Kfold   = lapply(folds, function(f) matrix(observation.transformation(Y[f$test,  , drop = FALSE], grid.type = "full", lag = j), p, p) )  
        train.Kfold  = lapply(folds, function(f)        observation.transformation(Y[f$train, , drop = FALSE], grid.type = grid.type, lag = j)     )  
        
        mean_sup = sapply(w.h.list, function(w_h){ mean(sapply(1:(K*R),function(k){ max(abs((test.Kfold[[k]] - eval.weights(w_h, train.Kfold[[k]], lag = j))[!as.logical(diag(p))]))} ))} ) 
    
        h.seq[which.min(mean_sup)]
        
      })
      
    }, future.seed = TRUE)   
  
    bw.list = matrix(bw.list, nrow = max.lag, dimnames = list(paste0("lag ", 1:max.lag), NULL))
  
    return(bw.list)
  }
}









