



est.results = function(data.1,data.2,bandwidth,from = 1,to = 12, eval.grid = NA){
  
  p  = dim(data.1)[2]- 4
  pd = dim(data.2)[2] - 4
  
  if(any(is.na(eval.grid))) {eval.grid =  (20:(136-20))/136} # every 15 Minute from -19:00 till +05:00: 4*5 + (4*24+1) + 4*5 = 20 + 97 + 20
  
  eval.tibble = tibble(TIME = hms::as_hms(seq(from = as.POSIXct("1970-01-01 00:00:00"),to   = as.POSIXct("1970-01-01 23:45:00"),by   = "15 min")))
  eval.tibble = eval.tibble  |>  add_row(TIME = as_hms("23:59:59"))
  
  result = future_sapply(from:to, function(m){
    
    bw.2      = bandwidth[1,m]
    bw.delta  = bandwidth[2,m]
    bw.1      = bandwidth[3,m]
    
    sample.2  = data.2  |>  filter(MONTH == month.name[m]) |> dplyr::select(4:dim(data.2)[2])
    sample.2  = sample.2[rowSums(is.na(sample.2)) == 0,]
    
    L.2 = tibble(TIME = colnames(sample.2), MEAN = apply(sample.2,2,mean,na.rm = T))
    
    L.2.eval          = eval.tibble
    L.2.eval$ESTIMATE = locPolSmootherC(x = (0:pd)/pd, y = L.2$MEAN, xeval = eval.grid, bw = bw.2, deg = 2, EpaK)$beta0
    L.2.eval$MONTH    = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    sample.1  = data.1  |>  filter(MONTH == month.name[m]) |> dplyr::select(4:dim(data.1)[2])
    sample.1  = sample.1[rowSums(is.na(sample.1)) == 0,]
    
    L.1        = tibble(TIME  = colnames(sample.1), MEAN = apply(sample.1,2,mean,na.rm = T))
    part2.hat  = locPolSmootherC(x = (0:pd)/pd,y = L.2$MEAN, xeval = (0:p)/p, bw = bw.2,deg = 2, EpaK)$beta0
    res        = part2.hat - L.1$MEAN 
    
    L.1.eval           = eval.tibble
    part2.hat.eval     = locPolSmootherC(x = (0:pd)/pd, y = L.2$MEAN, xeval= eval.grid ,bw = bw.2,deg = 2, EpaK)$beta0
    L.1.eval$ESTIMATE  = part2.hat.eval - locPolSmootherC(x = (0:p)/p, y = res, xeval = eval.grid, bw = bw.delta,deg = 2, EpaK)$beta0
    L.1.eval$MONTH     = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    L.delta.eval            = eval.tibble
    L.delta.eval$ESTIMATE   = locPolSmootherC(x = (0:p)/p, y = res,xeval = eval.grid ,bw = bw.delta,deg = 2, EpaK)$beta0 
    L.delta.eval$MONTH      = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    int.delta.eval          = eval.tibble
    int.delta.hat           = integrate(function(x)locPolSmootherC(x = (0:p)/p, y = res,xeval = x ,bw = bw.delta, deg = 2, EpaK)$beta0,lower = eval.grid[1],upper =  eval.grid[length(eval.grid)], stop.on.error = FALSE)$value/(eval.grid[length(eval.grid)]-eval.grid[1])  
 
    int.delta.eval$ESTIMATE = rep(int.delta.hat, times = length(eval.tibble))
    int.delta.eval$MONTH    = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    
    L.compare = tibble(TIME = colnames(sample.1), MEAN = apply(sample.1,2,mean,na.rm = T))
    
    L.compare.eval          = eval.tibble
    L.compare.eval$ESTIMATE = L.2.eval$ESTIMATE - locPolSmootherC(x = (0:p)/p, y = L.compare$MEAN, xeval = eval.grid, bw = bw.1, deg = 2, EpaK)$beta0
    L.compare.eval$MONTH    = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    int.part1  = integrate(function(x)locPolSmootherC(x = (0:p)/p  , y = L.compare$MEAN, xeval = x, bw = bw.1, deg = 2, EpaK)$beta0,lower = eval.grid[1],upper =  eval.grid[length(eval.grid)], stop.on.error = FALSE)$value/(eval.grid[length(eval.grid)]-eval.grid[1])  
    int.part2  = integrate(function(x)locPolSmootherC(x = (0:pd)/pd, y = L.2$MEAN  , xeval = x,     bw = bw.2, deg = 2, EpaK)$beta0,lower = eval.grid[1],upper =  eval.grid[length(eval.grid)], stop.on.error = FALSE)$value/(eval.grid[length(eval.grid)]-eval.grid[1])  
    
    L.compare.eval$Integral = rep(int.part2 - int.part1, times = length(eval.tibble))
    L.2.eval$Integral       = rep(int.part2, times = length(eval.tibble))
    
    
    # integral of sparse set with residual based estimator
    L.1.eval$Integral       = L.2.eval$Integral - int.delta.eval$ESTIMATE
    
    
    L.compare.eval$centered = L.compare.eval$ESTIMATE - L.compare.eval$Integral
    L.delta.eval$centered   = L.delta.eval$ESTIMATE   - int.delta.eval$ESTIMATE
    
    list(delta = L.delta.eval, mu.part2 = L.2.eval, mu.part1 = L.1.eval, delta.int = int.delta.eval, delta.compare = L.compare.eval) 
    
  }, future.seed = T)
  
  result.delta   = bind_rows(result[1,])
  result.2       = bind_rows(result[2,])
  result.1       = bind_rows(result[3,])
  result.int     = bind_rows(result[4,])
  result.compare = bind_rows(result[5,])
  
  
  return(list(delta = result.delta, part1 = result.1, part2 = result.2, delta.int = result.int, delta.compare =  result.compare) )
}






q.month = function(sample1, sample2, bandwidth, est, cov, from = 1, to = 12,alpha = NA, B = 1000, depend = F, int = F, constant = NA){
  
  if(any(is.na(alpha))){    alpha    = numeric(length(from:to)) + 0.9} # Default 90%
  if(any(is.na(constant))){ constant = numeric(length(from:to))}       # H0: delta = 0
  
  res = sapply(from:to,  function(m){
    
    sample2  = sample2  |>  filter(MONTH %in% month.name[m]) |> dplyr::select(4:dim(sample2)[2])
    sample2  = sample2[rowSums(is.na(sample2)) == 0,]
    sample1 = sample1  |>  filter(MONTH %in% month.name[m])|> dplyr::select(4:dim(sample1)[2])
    sample1 = sample1[rowSums(is.na(sample1)) == 0,]
    
    eval.grid =  0:96/96 
    
    list = q.MB(sample1, sample2, 
                est$part1 |> filter(MONTH %in% month.name[m]), 
                est$part2 |> filter(MONTH %in% month.name[m]), 
                unlist(cov[[m]]), 
                eval.grid, bandwidth[,m], alpha = alpha[m], B = B, depend = depend, int = int, H0 = constant[m])
    
    print(paste0("done: month ", month.name[m], ": ", list$quantile ))
    
    if(int == F){obs = sqrt(dim(sample1)[1]) * max(abs((est$delta |> filter(MONTH %in% month.name[m]))$ESTIMATE -  constant[m] )/ sqrt(unlist(cov[[m]]))) }
    else{        obs = sqrt(dim(sample1)[1]) * max(abs((est$delta |> filter(MONTH %in% month.name[m]))$ESTIMATE - (est$delta.int |> filter(MONTH %in% month.name[m]))$ESTIMATE) / sqrt(unlist(cov[[m]]))) }
    
    p.value = sum(list$sample >= obs)/length(list$sample)
    c(list$quantile, p.value)})
  
  colnames(res) = month.name[from:to]
  rownames(res) = c( paste0(min(alpha),"-quantile") , "p.val" )
  
  return(res)
}


q.MB = function(Y1, Y2, Y1.est,  Y2.est, cov, x.eval, bandwidth, alpha = 0.9, B = 1000, depend = F, int = F, H0 = 0){
  
  grid.1  = 0:(dim(Y1)[2]-1)/(dim(Y1)[2]-1)
  grid.2  = 0:(dim(Y2)[2]-1)/(dim(Y2)[2]-1)
  
  bandwidth = round(bandwidth*34/24, digits = 2) #from 34h grid to 24h grid
  
  w.1 = locPolWeights(x= grid.1, bw = bandwidth[2], deg = 2, xeval= x.eval, kernel=EpaK)$locWeig
  w.2 = locPolWeights(x= grid.2, bw = bandwidth[1], deg = 2, xeval= grid.1,  kernel=EpaK)$locWeig
  
  l1 = list(sample1 = t(Y1), weight.s = w.1, sample1.mean = unlist(colMeans(Y1)), est1 = Y1.est$ESTIMATE)
  l2 = list(sample2 = t(Y2), weight.d = w.2, sample2.mean = unlist(colMeans(Y2)), est2 = Y2.est$ESTIMATE) 
  
  
  l.1 = list(S.1 = t(Y1), weight.1 = w.1, est.1 = Y1.est$ESTIMATE)
  l.2 = list(S.2 = t(Y2), weight.2 = w.2, est.2 = Y2.est$ESTIMATE) 
  
  
  if(int == F){
    l.1 = modifyList(l.1, list(constant = H0))
    l.2 = modifyList(l.2, list(constant = H0))
  }

  
  if(int == T){
    
    l.1 = modifyList(l.1, list(int.1 = Y1.est$Integral[1])) # easier to code (integral estimation in est.result)
    l.2 = modifyList(l.2, list(int.2 = Y2.est$Integral[1])) # easier to code (integral estimation in est.result)
    
  }
  
  delta.Bootstrap = unlist(future_lapply(1:B, function(i) {MB(l.1, l.2, cov, dependent = depend, int = int) },future.seed = T))
  
  q = quantile(delta.Bootstrap, probs = alpha, Type = 2, na.rm = T)
  
  return(list(quantile = q, sample = delta.Bootstrap))
} 


MB = function(list.1, list.2, cov, dependent = F, int = F){
  
  S.1   = list.1[[1]];         S.2   = list.2[[1]]
  w.1   = list.1[[2]];         w.2   = list.2[[2]]  
  est.1 = unlist(list.1[[3]]); est.2 = unlist(list.2[[3]])
  
  if(int == T){ int.1  = list.1[[4]] ; int.d = list.2[[4]] }
  else{         H0.1   = list.1[[4]] ; H0.2  = list.2[[4]] }
  
  n1 = dim(S.1)[2]; p  = dim(S.1)[1]
  n2 = dim(S.2)[2]; pd = dim(S.2)[1]
  
  
  if(dependent == T){
    
    ln.func = function(n){floor(2*n^(1/3))}
    k1 = function(h, n, func) {
      L = func(n)
      ifelse(abs(h) < L, 1 / (2*L - 1), 0)
    }
    
    q.n1 = 1/(2*ln.func(n1)-1)
    q.n2 = 1/(2*ln.func(n2)-1)
    
    w.n1  = rnorm(3*n1, mean = 0, sd = 1/sqrt(q.n1))
    w.n2  = rnorm(3*n2, mean = 0, sd = 1/sqrt(q.n2))
    
    g.n1  = numeric(n1)
    g.n2  = numeric(n2)
    
    for(j in 1:n1){
      g.n1[j] = sum(sapply((-ln.func(n1)):ln.func(n1), function(h) k1(h, n1, ln.func)) * w.n1[j:(j+2*ln.func(n1))])
    }
    for(j in 1:n2){
      g.n2[j] = sum(sapply((-ln.func(n2)):ln.func(n2), function(h) k1(h, n2, ln.func)) * w.n2[j:(j+2*ln.func(n2))])
    }
    g.n1  = g.n1  - mean(g.n1)   # Remark 1 of Bücher: -1
    g.n2  = g.n2  - mean(g.n2)
    
  }else{
    g.n1  = rnorm(n1)
    g.n2  = rnorm(n2)
  }
  
  if(int == T){ S.1.int = sapply(1:n1, function(i){ mean(w.1 %*% S.1[, i])})
                S.2.int = sapply(1:n2, function(i){ mean(w.1 %*% drop(w.2 %*% S.2[, i]))}) }
  
  
  if(int == F){return(max(abs((1/sqrt(n1-1)*( (w.1 %*% S.1 - est.1) - H0.1 )) %*% g.n1  -
                                (n1/(n2*sqrt(n1-1))*( (w.1 %*% drop(w.2 %*% S.2) - est.2) - H0.2 )) %*% g.n2)/sqrt(cov))) }
  
  else{        return(max(abs(1/sqrt(n1-1)*(sweep(w.1 %*% S.1, 2, S.1.int, "-") - est.1 + int.1 ) %*% g.n1  -
                                (n1/(n2*sqrt(n1-1))*(w.1 %*% drop(sweep(w.2 %*% S.2, 2, S.2.int, "-")) - est.2 + int.d)) %*% g.n2)/sqrt(cov)))}
}






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
      
      l.s = list(S.s = t(S.s),  w.s = w.s, int.s = (est$part1|> filter(MONTH %in% month.name[m]))$Integral[1])
      l.d = list(S.d = t(S.d),  w.d = w.d, int.d = (est$part2 |> filter(MONTH %in% month.name[m]))$Integral[1]) 
      
      
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
