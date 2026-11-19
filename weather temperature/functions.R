
#############################################################################################################################
################ Mean and difference function estimation ####################################################################
#############################################################################################################################


est.results = function(data.s,data.d,bandwidth,from = 1,to = 12, eval.grid = NA){
  
  p  = dim(data.s)[2] - 4
  pd = dim(data.d)[2] - 4
  
  n  = dim(data.s)[1]
  nd = dim(data.d)[1]
  
  if(any(is.na(eval.grid))) {eval.grid =  (20:(136-20))/136} # every 15 Minute from -19:00 till +05:00: 4*5 + (4*24+1) + 4*5 = 20 + 97 + 20
  
  eval.tibble = tibble(TIME = hms::as_hms(seq(from = as.POSIXct("1970-01-01 00:00:00"),to   = as.POSIXct("1970-01-01 23:45:00"),by   = "15 min")))
  eval.tibble = eval.tibble |> add_row(TIME = as_hms("23:59:59"))
  
  result = future_sapply(from:to, function(m){
    
    bw.d      = bandwidth[1,m]
    bw.delta  = bandwidth[2,m]
    bw.s      = bandwidth[3,m]
    
    S.d  = data.d |> filter(MONTH == month.name[m]) |> dplyr::select(4:dim(data.d)[2])
    S.d  = S.d[rowSums(is.na(S.d)) == 0,]
    
    L.d = tibble(TIME = colnames(S.d), MEAN = apply(S.d,2,mean,na.rm = T))
    
    L.d.eval          = eval.tibble
    L.d.eval$ESTIMATE = locPolSmootherC(x = (0:pd)/pd, y = L.d$MEAN, xeval = eval.grid, bw = bw.d, deg = 2, EpaK)$beta0
    L.d.eval$MONTH    = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    S.s  = data.s |> filter(MONTH == month.name[m]) |> dplyr::select(4:dim(data.s)[2])
    S.s  = S.s[rowSums(is.na(S.s)) == 0,]
    
    L.s = tibble(TIME  = colnames(S.s), MEAN = apply(S.s,2,mean,na.rm = T))
    dense.hat               = locPolSmootherC(x = (0:pd)/pd,y = L.d$MEAN, xeval = (0:p)/p, bw = bw.d,deg = 2, EpaK)$beta0
    res                     = dense.hat - L.s$MEAN 
    
    L.s.eval           = eval.tibble
    dense.hat.eval          = locPolSmootherC(x = (0:pd)/pd, y = L.d$MEAN, xeval= eval.grid ,bw = bw.d,deg = 2, EpaK)$beta0
    L.s.eval$ESTIMATE  = dense.hat.eval - locPolSmootherC(x = (0:p)/p, y = res, xeval = eval.grid, bw = bw.delta,deg = 2, EpaK)$beta0
    L.s.eval$MONTH     = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    L.delta.eval            = eval.tibble
    L.delta.eval$ESTIMATE   = locPolSmootherC(x = (0:p)/p, y = res,xeval = eval.grid ,bw = bw.delta,deg = 2, EpaK)$beta0 
    L.delta.eval$MONTH      = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    int.delta.eval          = eval.tibble
    int.delta.hat           = integrate(function(x)locPolSmootherC(x = (0:p)/p, y = res,xeval = x ,bw = bw.delta,deg = 2, EpaK)$beta0,lower = eval.grid[1],upper =  eval.grid[length(eval.grid)], stop.on.error = FALSE)$value/(eval.grid[length(eval.grid)]-eval.grid[1])  
    
    int.delta.eval$ESTIMATE = rep(int.delta.hat, times = length(eval.tibble))
    int.delta.eval$MONTH    = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    
    L.compare = tibble(TIME = colnames(S.s), MEAN = apply(S.s,2,mean,na.rm = T))
    
    L.compare.eval          = eval.tibble
    L.compare.eval$ESTIMATE = L.d.eval$ESTIMATE - locPolSmootherC(x = (0:p)/p, y = L.compare$MEAN, xeval = eval.grid, bw = bw.s, deg = 2, EpaK)$beta0
    L.compare.eval$MONTH    = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    
    int.sparse              = integrate(function(x)locPolSmootherC(x = (0:p)/p  , y = L.compare$MEAN, xeval = x, bw = max(bw.s,0.08), deg = 2, EpaK)$beta0,lower = eval.grid[1],upper =  eval.grid[length(eval.grid)], stop.on.error = FALSE)$value/(eval.grid[length(eval.grid)]-eval.grid[1])  
    int.dense               = integrate(function(x)locPolSmootherC(x = (0:pd)/pd, y = L.d$MEAN  , xeval = x, bw = bw.d ,           deg = 2, EpaK)$beta0,lower = eval.grid[1],upper =  eval.grid[length(eval.grid)], stop.on.error = FALSE)$value/(eval.grid[length(eval.grid)]-eval.grid[1])  
    
    L.compare.eval$Integral = rep(int.dense - int.sparse, times = length(eval.tibble))
    L.d.eval$Integral       = rep(int.dense, times = length(eval.tibble))

    
    # integral of sparse set with residual based estimator
    L.s.eval$Integral       = L.d.eval$Integral - int.delta.eval$ESTIMATE
    
    
    L.compare.eval$centered = L.compare.eval$ESTIMATE - L.compare.eval$Integral
    L.delta.eval$centered   = L.delta.eval$ESTIMATE   - int.delta.eval$ESTIMATE
    
    list(delta = L.delta.eval, mu.dense = L.d.eval, mu.sparse = L.s.eval, delta.int = int.delta.eval, delta.compare = L.compare.eval) 
    
  }, future.seed = T)
  
  result.delta    = bind_rows(result[1,])
  result.d        = bind_rows(result[2,])
  result.s        = bind_rows(result[3,])
  result.int      = bind_rows(result[4,])
  result.compare  = bind_rows(result[5,])
  
  
  return(list(delta = result.delta, sparse = result.s, dense = result.d, delta.int = result.int, delta.compare =  result.compare) )
}






#######################################################################################################
###### Multiplier Bootstrap    ########################################################################
#######################################################################################################


## For Multiplier Bootstrap if integral = T:
P.Cov = function(cov, eval.type = "full"){
  if(eval.type == "full"){cov.col.mean = matrix(rep(apply(cov,2,mean), each = dim(cov)[1]),ncol = dim(cov)[1], byrow = T)
                          cov.row.mean = matrix(rep(apply(cov,1,mean), each = dim(cov)[1]),ncol = dim(cov)[1])
                          cov.full.mean = rep(mean(cov), each = length(cov), ncol = dim(cov)[1])
                          cov.int =  - cov.col.mean - cov.row.mean + cov.full.mean 
  return(cov.int)}
  else if(eval.type == "diagonal"){
    return(-apply(cov,2,mean)-apply(cov,1,mean) + mean(cov))
    }
}



q.month = function(data.s, data.d, bandwidth, est, cov, from = 1, to = 12, alpha = NA, B = 1000, depend = F, int = F, constant = NA){
  
  if(any(is.na(alpha))){    alpha    = numeric(length(from:to)) + 0.9} # Default 90%
  if(any(is.na(constant))){ constant = numeric(length(from:to))}       # H0: delta = 0
  
  res = sapply(from:to,  function(m){
    
    S.d  = data.d |> filter(MONTH %in% month.name[m]) |> dplyr::select(4:dim(data.d)[2])
    S.d  = S.d[rowSums(is.na(S.d)) == 0,]
    S.s = data.s |> filter(MONTH %in% month.name[m])|> dplyr::select(4:dim(data.s)[2])
    S.s = S.s[rowSums(is.na(S.s)) == 0,]
    
    eval.grid =  0:96/96 
    
    list = q.MB(S.s, S.d, 
                est$sparse |> filter(MONTH %in% month.name[m]), 
                est$dense |> filter(MONTH %in% month.name[m]), 
                unlist(cov[[m]]), 
                eval.grid, bandwidth[,m], alpha = alpha[m], B = B, depend = depend, int = int, H0 = constant[m])
    
    print(paste0("done: month ", month.name[m], ": ", list$quantile ))
    
    if(int == F){  obs.cv = sqrt(dim(S.s)[1]) * max(abs((est$delta    |> filter(MONTH %in% month.name[m]))$ESTIMATE -  constant[m] )/ sqrt(unlist(cov[[m]]))) 
    }else{         obs.cv = sqrt(dim(S.s)[1]) * max(abs((est$delta    |> filter(MONTH %in% month.name[m]))$ESTIMATE - (est$delta.int |> filter(MONTH %in% month.name[m]))$ESTIMATE) / sqrt(unlist(cov[[m]])))}
    
    p.value = sum(list$sample >= obs.cv)/length(list$sample)
    
    c(list$quantile, p.value)})
  
  colnames(res) = month.name[from:to]
  rownames(res) = c( paste0(min(alpha),"-quantile") , "p.val" )
  
  return(res)
}


q.MB = function(Ys, Yd, Ys.est,  Yd.est, cov, x.eval, bandwidth, alpha = 0.9, B = 1000, depend = F, int = F, H0 = 0){
  
  grid.s  = 0:(dim(Ys)[2]-1)/(dim(Ys)[2]-1) # 0:24/24
  grid.d  = 0:(dim(Yd)[2]-1)/(dim(Yd)[2]-1) # 0:144/144
  
  bandwidth = round(bandwidth*34/24, digits = 2) #from 34h grid to 24h grid
  
  w.s = locPolWeights(x= grid.s, bw = bandwidth[2], deg = 2, xeval= x.eval, kernel=EpaK)$locWeig
  w.d = locPolWeights(x= grid.d, bw = bandwidth[1], deg = 2, xeval= grid.s, kernel=EpaK)$locWeig
  
  l.s = list(S.s = t(Ys), weight.s = w.s, est.s = Ys.est$ESTIMATE)
  l.d = list(S.d = t(Yd), weight.d = w.d, est.d = Yd.est$ESTIMATE) 
  
  if(int == F){
    l.s = modifyList(l.s, list(constant = H0))
    l.d = modifyList(l.d, list(constant = H0))
  }
  
  if(int == T){
    
    l.s = modifyList(l.s, list(int.s = Ys.est$Integral[1])) # easier to code (integral estimation in est.result)
    l.d = modifyList(l.d, list(int.d = Yd.est$Integral[1])) # easier to code (integral estimation in est.result)
    
  }
  
  delta.Bootstrap = unlist(future_lapply(1:B, function(i) {MB(l.s,l.d,cov, dependent = depend,int = int) },future.seed = T))
  
  q = quantile(delta.Bootstrap, probs = alpha, Type = 2, na.rm = T)
  
  return(list(quantile = q, sample = delta.Bootstrap))
} 


MB = function(list.s, list.d, cov, dependent = F, int = F){
  
  S.s   = list.s[[1]];         S.d   = list.d[[1]]
  w.s   = list.s[[2]];         w.d   = list.d[[2]]  
  est.s = unlist(list.s[[3]]); est.d = unlist(list.d[[3]])
  
  if(int == T){ int.s = list.s[[4]] ; int.d = list.d[[4]] }
  else{         H0.s  = list.s[[4]] ; H0.d  = list.d[[4]] }
  
  n  = dim(S.s)[2]; p  = dim(S.s)[1]
  nd = dim(S.d)[2]; pd = dim(S.d)[1]
  
  
  if(dependent == T){
    
    ln.func = function(n){floor(2*n^(1/3))}
    k1 = function(h, n, func) {
      L = func(n)
      ifelse(abs(h) < L, 1 / (2*L - 1), 0)
    }
    
    q.n = 1/(2*ln.func(n)-1)
    q.nd = 1/(2*ln.func(nd)-1)
    
    w.n  = rnorm(3*n, mean = 0, sd = 1/sqrt(q.n))
    w.nd = rnorm(3*nd, mean = 0, sd = 1/sqrt(q.nd))
    
    g.n  = numeric(n)
    g.nd = numeric(nd)
    
    for(j in 1:n){
      g.n[j] = sum(sapply((-ln.func(n)):ln.func(n), function(h) k1(h, n, ln.func)) * w.n[j:(j+2*ln.func(n))])
    }
    for(j in 1:nd){
      g.nd[j] = sum(sapply((-ln.func(nd)):ln.func(nd), function(h) k1(h, nd, ln.func)) * w.nd[j:(j+2*ln.func(nd))])
    }
    g.n  = g.n  - mean(g.n)   # Remark 1 of Bücher: -1
    g.nd = g.nd - mean(g.nd)
    
  }else{
    g.n  = rnorm(n)
    g.nd = rnorm(nd)
  }
  
  if(int == T){ S.s.int = sapply(1:n,  function(i){ mean(w.s %*% S.s[, i])})
                S.d.int = sapply(1:nd, function(i){ mean(w.s %*% drop(w.d %*% S.d[, i]))}) }
  
  
  if(int == F){return(max(abs((1/sqrt(n-1)*( (w.s %*% S.s - est.s) - H0.s )) %*% g.n  -
                                (n/(nd*sqrt(n-1))*( (w.s %*% drop(w.d %*% S.d) - est.d) - H0.d )) %*% g.nd)/sqrt(cov))) }
  
  else{        return(max(abs(1/sqrt(n-1)*(sweep(w.s %*% S.s, 2, S.s.int, "-") - est.s + int.s) %*% g.n  -
                                (n/(nd*sqrt(n-1))*(w.s %*% drop(sweep(w.d %*% S.d, 2, S.d.int, "-")) - est.d + int.d)) %*% g.nd)/sqrt(cov)))}
}





CB = function(data.s, estimation, cov, q.list, center = F){
  
  est.cv = estimation$delta
  if(center == T){ est.cv$ESTIMATE = est.cv$centered }
  
  UP = c()
  LO = c()
  
  for(m in 1:12){S.s = data.s |> filter(MONTH %in% month.name[m])|> dplyr::select(4:dim(data.s)[2])
                 S.s = S.s[rowSums(is.na(S.s)) == 0,]
                 
                 up = c(est.cv |> filter(MONTH ==  month.name[m]) |> select(ESTIMATE))$ESTIMATE 
                 lo = c(est.cv |> filter(MONTH ==  month.name[m]) |> select(ESTIMATE))$ESTIMATE

                 UP = c(UP, up  + q.list[m] * sqrt(cov[[m]]/dim(S.s)[1])) 
                 LO = c(LO, lo  - q.list[m] * sqrt(cov[[m]]/dim(S.s)[1])) 
  }
  
  
  UP = unlist(UP)
  UP = unname(UP)
  
  LO = unlist(LO)
  LO = unname(LO)
  
  est.cv$UP = UP
  est.cv$LO = LO
  
  return(est.cv |> select( c(-4) ) )
}





############################
### Integral delta Test ####
############################

test.int = function(data.s, data.d, int.est, var.s, var.d, from = 1, to = 12, test = "two-sided", alpha = NA, approx = "normal",...){
  
  if(any(is.na(alpha))){ alpha = numeric(length(from:to)) + 0.9} # Default 90%
  
  nd   = unlist(lapply(from:to, function(m) {S.d  = data.d |> filter(MONTH %in% month.name[m]) |> dplyr::select(4:dim(data.d)[2])
  S.d  = S.d[rowSums(is.na(S.d)) == 0,]
  return(dim(S.d)[1])}))
  
  n   = unlist(lapply(from:to, function(m) {S.s = data.s |> filter(MONTH %in% month.name[m])|> dplyr::select(4:dim(data.s)[2])
  S.s = S.s[rowSums(is.na(S.s)) == 0,]
  return(dim(S.s)[1])}))
  
  test.stat = unlist(lapply(from:to, function(m) { sqrt(n[m]) * int.est[m] / sqrt(var.s[m] + n[m]/nd[m] * var.d[m])}))
  
  if( !(approx %in% c("normal","MB"))){ stop("Quantile approximation not implemented.") }
  
  if(approx == "normal"){
    if(test == "two-sided"){  p.val = unlist(lapply(from:to, function(m) {2*pnorm(abs(test.stat[m]),0,1,lower.tail = F) }))}
    if(test == "right-sided"){p.val = unlist(lapply(from:to, function(m) {  pnorm(test.stat[m],0,1,lower.tail = F) }))}
    if(test == "left-sided"){ p.val = unlist(lapply(from:to, function(m) {  pnorm(test.stat[m],0,1,lower.tail = T) }))}
    
    if(test == "two-sided"){ confInterval = lapply(from:to, function(m) {c(int.est[m]-sqrt(var.s[m] + n[m]/nd[m] * var.d[m])*qnorm((1 + alpha[m])/2) / sqrt(n[m]),int.est[m] + sqrt(var.s[m] + n[m]/nd[m] * var.d[m])*qnorm((1 + alpha[m])/2) / sqrt(n[m]))})}
    else{                    confInterval = lapply(from:to, function(m) {c(int.est[m]-sqrt(var.s[m] + n[m]/nd[m] * var.d[m])*qnorm(alpha[m]) / sqrt(n[m]),int.est[m] + sqrt(var.s[m] + n[m]/nd[m] * var.d[m])*qnorm(alpha[m]) / sqrt(n[m]))})}
    
  }
  else{
    
    MB.int = function(l.s, l.d, cov, dependent = F){
      
      S.s = l.s[[1]]; w.s = l.s[[2]]; int.est.s = unlist(l.s[[3]]) 
      S.d = l.d[[1]]; w.d = l.d[[2]]; int.est.d = unlist(l.d[[3]]) 
      
      n  = dim(S.s)[2]; p  = dim(S.s)[1]
      nd = dim(S.d)[2]; pd = dim(S.d)[1]
      
      if(dependent == T){
        
        ln.func = function(n){floor(2*n^(1/3))}
        k1 = function(h, n, func) {
          L = func(n)
          ifelse(abs(h) < L, 1 / (2*L - 1), 0)
        }
        
        q.n = 1/(2*ln.func(n)-1);                      q.nd = 1/(2*ln.func(nd)-1)
        w.n  = rnorm(3*n, mean = 0, sd = 1/sqrt(q.n)); w.nd = rnorm(3*nd, mean = 0, sd = 1/sqrt(q.nd))
        g.n  = numeric(n);                             g.nd = numeric(nd)
        
        for(j in 1:n){   g.n[j] = sum(sapply((-ln.func(n)):ln.func(n),   function(h) k1(h, n, ln.func))  * w.n[j:(j+2*ln.func(n))])  }
        for(j in 1:nd){ g.nd[j] = sum(sapply((-ln.func(nd)):ln.func(nd), function(h) k1(h, nd, ln.func)) * w.nd[j:(j+2*ln.func(nd))])}
        
        g.n  = g.n  - mean(g.n) ; g.nd = g.nd - mean(g.nd) }
      else{ 
        g.n  = rnorm(n);          g.nd = rnorm(nd)}
      
      
      S.s.int = sapply(1:n,  function(i){ mean(w.s %*% S.s[, i])})
      S.d.int = sapply(1:nd, function(i){ mean(w.s %*% drop(w.d %*% S.d[, i]))}) 
      
      
      return(  (  (1/sqrt(n-1)* (S.s.int - int.est.s) %*% g.n)  - (n/(nd*sqrt(n-1))* (S.d.int - int.est.d)  %*% g.nd)  )/ sqrt(cov)  )
      
    }
    
    
    args.MB = as.list(...)
    est     = args.MB[[1]]; bandwidth = args.MB[[2]]; B = args.MB[[3]]; depend = args.MB[[4]]
    
    bandwidth = round(bandwidth*34/24, digits = 2) #from 34h grid to 24h grid
    
    res = sapply(from:to,  function(m){
      
      S.d  = data.d |> filter(MONTH %in% month.name[m]) |> dplyr::select(4:dim(data.d)[2])
      S.d  = S.d[rowSums(is.na(S.d)) == 0,]
      S.s  = data.s |> filter(MONTH %in% month.name[m])|> dplyr::select(4:dim(data.s)[2])
      S.s  = S.s[rowSums(is.na(S.s)) == 0,]
      
      x.eval =  0:96/96 
      
      grid.s  = 0:(dim(S.s)[2]-1)/(dim(S.s)[2]-1) # 0:24/24
      grid.d  = 0:(dim(S.d)[2]-1)/(dim(S.d)[2]-1) # 0:144/144
      
      w.s = locPolWeights(x= grid.s, bw = bandwidth[2,m], deg = 2, xeval= x.eval, kernel=EpaK)$locWeig
      w.d = locPolWeights(x= grid.d, bw = bandwidth[1,m], deg = 2, xeval= grid.s, kernel=EpaK)$locWeig
      
      l.s = list(S.s = t(S.s),  w.s = w.s, int.s = (est$sparse|> filter(MONTH %in% month.name[m]))$Integral[1])
      l.d = list(S.d = t(S.d),  w.d = w.d, int.d = (est$dense |> filter(MONTH %in% month.name[m]))$Integral[1]) 
      

      cov = var.s[m] + n[m]/nd[m] * var.d[m]
      
      int.Bootstrap = unlist(future_lapply(1:B, function(i) { MB.int(l.s,l.d, cov, dependent = depend)  },future.seed = T))
      
      if(test == "two-sided"){    p.value = sum(abs(int.Bootstrap) >= abs(test.stat)[m] )/length(int.Bootstrap)
                                  q.upper       = quantile(int.Bootstrap, probs = (alpha[m]+1)/2,     Type = 2, na.rm = T)
                                  q.lower       = quantile(int.Bootstrap, probs = (1-(alpha[m]+1)/2), Type = 2, na.rm = T)
                                  result = c(p.value = p.value, quantile.LO = q.lower, quantile.UP = q.upper)}
      if(test == "right-sided"){  p.value = sum( int.Bootstrap     >=     test.stat[m]  )/length(int.Bootstrap)
                                  q       = quantile(int.Bootstrap, probs = alpha[m], Type = 2, na.rm = T)
                                  result = c(p.value = p.value, quantile = q)}
      if(test == "left-sided"){   p.value = sum( int.Bootstrap     <=     test.stat[m]  )/length(int.Bootstrap)
                                  q       = quantile(int.Bootstrap, probs = (1-alpha[m]), Type = 2, na.rm = T)
                                  result = c(p.value = p.value, quantile = q)}
      
      result
      
    })
  
    if(test == "two-sided"){   confInterval = lapply(from:to, function(m) {c(int.est[m] - sqrt(var.s[m] + n[m]/nd[m] * var.d[m])* abs(res[2,m]) / sqrt(n[m]),
                                                                           int.est[m] + sqrt(var.s[m] + n[m]/nd[m] * var.d[m])* abs(res[3,m]) / sqrt(n[m]))})  }
    
    if(test == "right-sided"){ confInterval = lapply(from:to, function(m) {c(-Inf, int.est[m] + sqrt(var.s[m] + n[m]/nd[m] * var.d[m])* abs(res[2,m]) / sqrt(n[m]))})}
    if(test == "left-sided"){  confInterval = lapply(from:to, function(m) {c(int.est[m] + sqrt(var.s[m] + n[m]/nd[m] * var.d[m])* abs(res[2,m]) / sqrt(n[m]), Inf)})} 
    
    p.val = res[1,]
    
  }
    
    names(confInterval) = month.name[from:to]
    return(list(pval = round(p.val, digits = 5), confInterval = confInterval))
}


