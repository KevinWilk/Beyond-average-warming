



est.results = function(data.1,data.2,bandwidth,from = 1,to = 12, eval.grid = NA){
  
  p  = dim(data.1)[2]- 4
  pd = dim(data.2)[2] - 4
  
  if(any(is.na(eval.grid))) {eval.grid =  (20:(136-20))/136} # every 15 Minute from -19:00 till +05:00: 4*5 + (4*24+1) + 4*5 = 20 + 97 + 20
  
  eval.tibble = tibble(TIME = hms::as_hms(seq(from = as.POSIXct("1970-01-01 00:00:00"),to   = as.POSIXct("1970-01-01 23:45:00"),by   = "15 min")))
  eval.tibble = eval.tibble %>% add_row(TIME = as_hms("23:59:59"))
  
  result = future_sapply(from:to, function(m){
    
    bw.2      = bandwidth[1,m]
    bw.delta  = bandwidth[2,m]
    bw.1      = bandwidth[3,m]
    
    sample.2  = data.2 %>% filter(MONTH == month.name[m]) |> dplyr::select(4:dim(data.2)[2])
    sample.2  = sample.2[rowSums(is.na(sample.2)) == 0,]
    
    L.2 = tibble(TIME = colnames(sample.2), MEAN = apply(sample.2,2,mean,na.rm = T))
    
    L.2.eval          = eval.tibble
    L.2.eval$ESTIMATE = locPolSmootherC(x = (0:pd)/pd, y = L.2$MEAN, xeval = eval.grid, bw = bw.2, deg = 2, EpaK)$beta0
    L.2.eval$MONTH    = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    sample.1  = data.1 %>% filter(MONTH == month.name[m]) |> dplyr::select(4:dim(data.1)[2])
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
    int_delta_hat           = integrate(function(x)locPolSmootherC(x = (0:p)/p, y = res,xeval = x ,bw = bw.delta, deg = 2, EpaK)$beta0,lower = eval.grid[1],upper =  eval.grid[length(eval.grid)], stop.on.error = FALSE)$value/(eval.grid[length(eval.grid)]-eval.grid[1])  
 
    int.delta.eval$ESTIMATE = rep(int_delta_hat, times = length(eval.tibble))
    int.delta.eval$MONTH    = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    
    L.compare = tibble(TIME = colnames(sample.1), MEAN = apply(sample.1,2,mean,na.rm = T))
    
    L.compare.eval          = eval.tibble
    L.compare.eval$ESTIMATE = L.2.eval$ESTIMATE - locPolSmootherC(x = (0:p)/p, y = L.compare$MEAN, xeval = eval.grid, bw = bw.1, deg = 2, EpaK)$beta0
    L.compare.eval$MONTH    = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    int_part1  = integrate(function(x)locPolSmootherC(x = (0:p)/p  , y = L.compare$MEAN, xeval = x, bw = bw.1, deg = 2, EpaK)$beta0,lower = eval.grid[1],upper =  eval.grid[length(eval.grid)], stop.on.error = FALSE)$value/(eval.grid[length(eval.grid)]-eval.grid[1])  
    int_part2  = integrate(function(x)locPolSmootherC(x = (0:pd)/pd, y = L.2$MEAN  , xeval = x,     bw = bw.2, deg = 2, EpaK)$beta0,lower = eval.grid[1],upper =  eval.grid[length(eval.grid)], stop.on.error = FALSE)$value/(eval.grid[length(eval.grid)]-eval.grid[1])  
    
    L.compare.eval$Integral = rep(int_part2 - int_part1, times = length(eval.tibble))
    
    L.compare.eval$centered = L.compare.eval$ESTIMATE - L.compare.eval$Integral
    L.delta.eval$centered   = L.delta.eval$ESTIMATE   - int.delta.eval$ESTIMATE
    
    list(delta = L.delta.eval, mu_part2 = L.2.eval, mu_part1 = L.1.eval, delta_int = int.delta.eval, delta_compare = L.compare.eval) 
    
  }, future.seed = T)
  
  result.delta   = bind_rows(result[1,])
  result.2       = bind_rows(result[2,])
  result.1       = bind_rows(result[3,])
  result.int     = bind_rows(result[4,])
  result.compare = bind_rows(result[5,])
  
  
  return(list(delta = result.delta, part1 = result.1, part2 = result.2, delta_int = result.int, delta_compare =  result.compare) )
}






q.month = function(sample1, sample2, bandwidth, est, cov, from = 1, to = 12,alpha = NA, B = 1000, depend = F, int = F, constant = NA, grid = "sparse to dense"){
  
  if(any(is.na(alpha))){    alpha    = numeric(length(from:to)) + 0.9} # Default 90%
  if(any(is.na(constant))){ constant = numeric(length(from:to))}       # H0: delta = 0
  
  res = sapply(from:to,  function(m){
    
    sample2  = sample2 %>% filter(MONTH %in% month.name[m]) |> dplyr::select(4:dim(sample2)[2])
    sample2  = sample2[rowSums(is.na(sample2)) == 0,]
    sample1 = sample1 %>% filter(MONTH %in% month.name[m])|> dplyr::select(4:dim(sample1)[2])
    sample1 = sample1[rowSums(is.na(sample1)) == 0,]
    
    eval.grid =  (20:(136-20))/136
    
    list = q.MB(sample1, sample2, est$part1%>%filter(MONTH %in% month.name[m]), est$part2%>%filter(MONTH %in% month.name[m]), 
                unlist(cov[[m]]), eval.grid, bandwidth[,m], alpha = alpha[m], B = B, depend = depend, int = int, H0 = constant[m], grid = grid)
    
    print(paste0("done: month ", month.name[m], ": ", list$quantile ))
    
    if(int == F){obs = sqrt(dim(sample1)[1]) * max(abs((est$delta%>%filter(MONTH %in% month.name[m]))$ESTIMATE -  constant[m] )/ sqrt(unlist(cov[[m]]))) }
    else{        obs = sqrt(dim(sample1)[1]) * max(abs((est$delta%>%filter(MONTH %in% month.name[m]))$ESTIMATE - (est$delta_int%>%filter(MONTH %in% month.name[m]))$ESTIMATE) / sqrt(unlist(cov[[m]]))) }
    
    p.value = sum(list$sample >= obs)/length(list$sample)
    c(list$quantile, p.value)})
  
  colnames(res) = month.name[from:to]
  rownames(res) = c( paste0(min(alpha),"-quantile") , "p.val" )
  
  return(res)
}


q.MB = function(Y1, Y2, Y1.est,  Y2.est, cov, x.eval, bandwidth, alpha = 0.9, B = 1000, depend = F, int = F, H0 = 0, grid = "sparse to dense"){
  
  if(grid == "sparse")      {          grid1  = c((0:(dim(Y1)[2]-1))/(dim(Y1)[2]-1))[6:(dim(Y1)[2]-5)]
  grid2  = c((0:(dim(Y1)[2]-1))/(dim(Y1)[2]-1))[6:(dim(Y1)[2]-5)]
  }else if(grid == "dense"){           grid1  = c((0:(dim(Y2)[2]-1))/(dim(Y2)[2]-1))[31:(dim(Y2)[2]-30)]
  grid2  = c((0:(dim(Y2)[2]-1))/(dim(Y2)[2]-1))[31:(dim(Y2)[2]-30)]
  }else if(grid == "sparse to dense"){ grid1  = c((0:(dim(Y1)[2]-1))/(dim(Y1)[2]-1))[6:(dim(Y1)[2]-5)]
  grid2  = c((0:(dim(Y2)[2]-1))/(dim(Y2)[2]-1))[31:(dim(Y2)[2]-30)]}
  
  
  w1 = locPolWeights(x= grid1, bw = bandwidth[2], deg = 2, xeval= x.eval, kernel=EpaK)$locWeig
  w2 = locPolWeights(x= grid2, bw = bandwidth[1], deg = 2, xeval= grid1, kernel=EpaK)$locWeig
  
  col.1.25h = !grepl("Before|After", colnames(Y1)) | colnames(Y1) == "After 00:00:00"
  col.2.25h = !grepl("Before|After", colnames(Y2)) | colnames(Y2) == "After 00:00:00"
  
  l1 = list(sample1 = t(Y1[,col.1.25h]), weight.s = w1, sample1.mean = unlist(colMeans(Y1[,col.1.25h])), est1 = Y1.est$ESTIMATE)
  l2 = list(sample2 = t(Y2[,col.2.25h]), weight.d = w2, sample2.mean = unlist(colMeans(Y2[,col.2.25h])), est2 = Y2.est$ESTIMATE) 
  
  if(int == F){
    l1 = modifyList(l1, list(constant = H0))
    l2 = modifyList(l2, list(constant = H0))
  }
  
  if(int == T){
    
    grid1  = (0:(dim(Y1)[2]-1))/(dim(Y1)[2]-1)
    grid2  = (0:(dim(Y2)[2]-1))/(dim(Y2)[2]-1)
    
    integral.2 = integrate(function(x)locPolSmootherC(x = grid2, y = colMeans(Y2), xeval = x ,bw = bandwidth[1],deg = 2, kernel = EpaK)$beta0,lower = x.eval[1],upper = x.eval[length(x.eval)])$value/(x.eval[length(x.eval)]-x.eval[1])
    
    res        = colMeans(Y1) - locPolSmootherC(x = grid2, y = colMeans(Y2), xeval = grid1, bw = bandwidth[1],deg = 2, EpaK)$beta0
    integral.1 = integral.2   + integrate(function(x)locPolSmootherC(x = grid1, y = res, xeval = x ,bw = bandwidth[2],deg = 2, kernel = EpaK)$beta0,lower = x.eval[1],upper =  x.eval[length(x.eval)])$value/(x.eval[length(x.eval)]-x.eval[1])
    
    l1 = modifyList(l1, list(int.s = integral.1))
    l2 = modifyList(l2, list(int.d = integral.2))
    
    l1 = modifyList(l1, list(bw1 = bandwidth[2]))    
    l2 = modifyList(l2, list(bw2 = bandwidth[1]))   
  }
  
  delta_Bootstrap = unlist(future_lapply(1:B, function(i) {MB(l1, l2, cov, dependent = depend,int = int) },future.seed = T))
  
  q = quantile(delta_Bootstrap, probs = alpha, Type = 2, na.rm = T)
  
  return(list(quantile = q, sample = delta_Bootstrap))
} 


MB = function(list1, list2, cov, dependent = F, int = F){
  
  sample1 = list1[[1]] 
  sample2 = list2[[1]]
  weights1 = list1[[2]]
  weights2 = list2[[2]]  
  Mean1    = list1[[3]] 
  Mean2    = list2[[3]] 
  T_val1   = unlist(list1[[4]]) 
  T_val2   = unlist(list2[[4]])
  
  int1 = list1[[5]] 
  int2 = list2[[5]]
  
  n1 = dim(sample1)[2]
  p1 = dim(sample1)[1]
  n2 = dim(sample2)[2]
  p2 = dim(sample2)[1]
  
  if(int == T){
    h1 = list1[[6]]
    h2 = list2[[6]]
  }
  
  if(dependent == T){
    
    
    ln_func = function(n){floor(2*n^(1/3))}
    k1 = function(h, n, func) {
      L = func(n)
      ifelse(abs(h) < L, 1 / (2*L - 1), 0)
    }
    
    
    q_n1 = 1/(2*ln_func(n1)-1)
    q_n2 = 1/(2*ln_func(n2)-1)
    
    w_n1 = rnorm(3*n1, mean = 0, sd = 1/sqrt(q_n1))
    w_n2 = rnorm(3*n2, mean = 0, sd = 1/sqrt(q_n2))
    
    g_n1 = numeric(n1)
    g_n2 = numeric(n2)
    
    for(j in 1:n1){
      g_n1[j] = sum(sapply((-ln_func(n1)):ln_func(n1), function(h) k1(h, n1, ln_func)) * w_n1[j:(j+2*ln_func(n1))])
    }
    for(j in 1:n2){
      g_n2[j] = sum(sapply((-ln_func(n2)):ln_func(n2), function(h) k1(h, n2, ln_func)) * w_n2[j:(j+2*ln_func(n2))])
    }
    g_n1 = g_n1 - mean(g_n1)   # Remark 1 of Bücher: -1
    g_n2 = g_n2 - mean(g_n2) }
  else{
    g_n1 = rnorm(n1)
    g_n2 = rnorm(n2)}
  
  if(int == T){f1_int = sapply(1:n1,  function(i){mean(locPolSmootherC((0:(p1-1))/(p1-1), sample1[, i], seq(0, 1, length.out = 1000), h1, 2, EpaK)$beta0)})
  f2_int = sapply(1:n2, function(i){mean(locPolSmootherC((0:(p2-1))/(p2-1), sample2[, i], seq(0, 1, length.out = 1000), h2, 2, EpaK)$beta0)})}
  
  if(int == F){return(max(abs((1/sqrt(n1-1)*( (weights1 %*% sample1 - T_val1) - int1)) %*% g_n1  -
                                (n1/(n2*sqrt(n1-1))*( (weights1 %*% drop(weights2 %*% sample2) - T_val2) - int2)) %*% g_n2)/sqrt(cov))) }
  
  else{        return(max(abs(1/sqrt(n1-1)*(sweep(weights1 %*% sample1, 2, f1_int, "-") - T_val1 + int1) %*% g_n1  -
                                (n1/(n2*sqrt(n1-1))*(weights1 %*% drop(sweep(weights2 %*% sample2, 2, f2_int, "-")) - T_val2 + int2)) %*% g_n2)/sqrt(cov)))}
}

