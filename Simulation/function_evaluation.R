est.results = function(data.sparse,data.dense,bandwidth,p.eval){
  
  p  = dim(data.sparse)[2]
  pd = dim(data.dense)[2] 
  
  eval.grid =  (1:p.eval-0.5)/p.eval  
  
  eval.tibble = tibble(x = eval.grid)
  
  bw.dense  = bandwidth[1]
  bw.delta  = bandwidth[2]
  
  L.dense = tibble(TIME = colnames(Y.d), MEAN = apply(Y.d,2,mean,na.rm = T))
  
  L.dense.eval          = eval.tibble
  L.dense.eval$ESTIMATE = locPolSmootherC(x = (1:pd-0.5)/pd, y = L.dense$MEAN, xeval = eval.grid, bw = bw.dense, deg = 2, EpaK)$beta0
  
  L.sparse = tibble(TIME  = colnames(Y.s), MEAN = apply(Y.s,2,mean,na.rm = T))
  dense.hat               = locPolSmootherC(x = (1:pd-0.5)/pd,y = L.dense$MEAN, xeval = (1:p-0.5)/p, bw = bw.dense,deg = 2, EpaK)$beta0
  res                     = L.sparse$MEAN - dense.hat  
  
  L.sparse.eval           = eval.tibble
  dense.hat.eval          = locPolSmootherC(x = (1:pd-0.5)/pd, y = L.dense$MEAN, xeval= eval.grid ,bw = bw.dense,deg = 2, EpaK)$beta0
  L.sparse.eval$ESTIMATE  = dense.hat.eval + locPolSmootherC(x = (1:p-0.5)/p, y = res, xeval = eval.grid, bw = bw.delta,deg = 2, EpaK)$beta0
  
  L.delta.eval            = eval.tibble
  L.delta.eval$ESTIMATE   = locPolSmootherC(x = (1:p-0.5)/p, y = res,xeval = eval.grid ,bw = bw.delta,deg = 2, EpaK)$beta0 
  
  int.delta.eval          = eval.tibble
  int_delta_hat           = integrate(function(x)locPolSmootherC(x = (1:p-0.5)/p, y = res,xeval = x ,bw = bw.delta,deg = 2, EpaK)$beta0,lower = 0,upper =  1, stop.on.error = FALSE)$value 
  
  int.delta.eval$ESTIMATE = rep(int_delta_hat, times = length(eval.tibble))
  
  return(list(delta = L.delta.eval, dense = L.dense.eval, sparse = L.sparse.eval, delta_int = int.delta.eval)) 
}



q.MB = function(Ys, Yd, Ys.est,  Yd.est, cov, x.eval, bandwidth, alpha = 0.9, B = 1000, depend = F, int = F, H0 = 0){
  
  grid.s  = (1:dim(Ys)[2]-0.5)/dim(Ys)[2]
  grid.d  = (1:dim(Yd)[2]-0.5)/dim(Yd)[2]
  
  w.s = locPolWeights(x= grid.s, bw = bandwidth[2], deg = 2, xeval= x.eval, kernel=EpaK)$locWeig
  w.d = locPolWeights(x= grid.d, bw = bandwidth[1], deg = 2, xeval= grid.s, kernel=EpaK)$locWeig
  
  ls = list(sample.s = t(Ys), weight.s = w.s, sample.s.mean = unlist(colMeans(Ys)), est.s = Ys.est$ESTIMATE)
  ld = list(sample.d = t(Yd), weight.d = w.d, sample.d.mean = unlist(colMeans(Yd)), est.d = Yd.est$ESTIMATE) 
  
  if(int == F){
    ls = modifyList(ls, list(constant = H0))
    ld = modifyList(ld, list(constant = H0))
  }
  
  if(int == T){
    
    integral.d = integrate(function(x)locPolSmootherC(x = grid.d, y = colMeans(Yd), xeval = x ,bw = bandwidth[1],deg = 2, kernel = EpaK)$beta0,lower = 0,upper =  1, stop.on.error = FALSE)$value
    
    res        = colMeans(Ys) - locPolSmootherC(x = grid.d, y = colMeans(Yd), xeval = grid.s, bw = bandwidth[1],deg = 2, EpaK)$beta0
    integral   = integral.d   + integrate(function(x)locPolSmootherC(x = grid.s, y = res, xeval = x ,bw = bandwidth[2],deg = 2, kernel = EpaK)$beta0,lower = 0,upper =  1, stop.on.error = FALSE)$value
    
    ls = modifyList(ls, list(int.s = integral))
    ld = modifyList(ld, list(int.d = integral.d))
    
    
    ls = modifyList(ls, list(bw.s = bandwidth[2]))    
    ld = modifyList(ld, list(bw.d = bandwidth[1]))   
  }
  
  delta_Bootstrap = unlist(future_lapply(1:B, function(i) {MB(ls,ld,cov, dependent = depend,int = int) },future.seed = T))
  
  q = quantile(delta_Bootstrap, probs = alpha, Type = 2)
  
  return(list(quantile = q, sample = delta_Bootstrap))
} 



MB = function(list1, list2, cov, dependent = F, int = F){
  
  sample                   = list1[[1]] 
  sample_d                 = list2[[1]]
  weights                  = list1[[2]]
  weights_d                = list2[[2]]  
  Mean                     = list1[[3]] 
  Mean_d                   = list2[[3]] 
  T_val1                   = unlist(list1[[4]]) 
  T_val2                   = unlist(list2[[4]])
  
  int1                   = list1[[5]] 
  int2                   = list2[[5]]
  
  n  = dim(sample)[2]
  p  = dim(sample)[1]
  nd = dim(sample_d)[2]
  pd = dim(sample_d)[1]
  
  if(int == T){
    h                      = list1[[6]]
    h_d                    = list2[[6]]
  }
  
  if(dependent == T){
    
    
    ln_func = function(n){floor(1.25*n^(1/3))}
    k1 = function(h, n, func) {
      L = func(n)
      ifelse(abs(h) < L, 1 / (2*L - 1), 0)
    }
    
    
    q_n = 1/(2*ln_func(n)-1)
    q_nd = 1/(2*ln_func(nd)-1)
    
    w_n  = rnorm(3*n, mean = 0, sd = 1/sqrt(q_n))
    w_nd = rnorm(3*nd, mean = 0, sd = 1/sqrt(q_nd))
    
    g_n  = numeric(n)
    g_nd = numeric(nd)
    
    for(j in 1:n){
      g_n[j] = sum(sapply((-ln_func(n)):ln_func(n), function(h) k1(h, n, ln_func)) * w_n[j:(j+2*ln_func(n))])
    }
    for(j in 1:nd){
      g_nd[j] = sum(sapply((-ln_func(nd)):ln_func(nd), function(h) k1(h, nd, ln_func)) * w_nd[j:(j+2*ln_func(nd))])
    }
    g_n  = g_n  - mean(g_n)   # Remark 1 of Bücher: -1
    g_nd = g_nd - mean(g_nd) }
  else{
    g_n  = rnorm(n)
    g_nd = rnorm(nd)}
  
  if(int == T){ f1_int = sapply(1:n,  function(i){mean(locPolSmootherC((1:p-0.5)/p,   sample[, i],   seq(0, 1, length.out = 1000), h,   2, EpaK)$beta0)})
                f2_int = sapply(1:nd, function(i){mean(locPolSmootherC((1:pd-0.5)/pd, sample_d[, i], seq(0, 1, length.out = 1000), h_d, 2, EpaK)$beta0)})}
  
  if(int == F){return(max(abs((1/sqrt(n-1)*( (weights %*% sample - T_val1) - int1)) %*% g_n  -
                                (n/(nd*sqrt(n-1))*( (weights %*% drop(weights_d %*% sample_d) - T_val2) - int2)) %*% g_nd)/sqrt(cov))) }
  
  else{        return(max(abs(1/sqrt(n-1)*(sweep(weights %*% sample, 2, f1_int, "-") - T_val1 + int1) %*% g_n  -
                                (n/(nd*sqrt(n-1))*(weights %*% drop(sweep(weights_d %*% sample_d, 2, f2_int, "-")) - T_val2 + int2)) %*% g_nd)/sqrt(cov)))}
}
