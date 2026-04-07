########################################################################################################################
################ Mean and Difference function: Cross Validation ########################################################
########################################################################################################################

k.fold.hc.cv = function(sample, h.seq, deg = 2, diff = F, K = 5,...){
  
  p = dim(sample)[2] - 2
  n = dim(sample)[1] 
  
  n.years = length(unique(sample$Year))
  grp.ind = rep(seq_len(K), times = c(rep(n.years %/% K, K - 1), n.years %/% K + n.years %% K))
  grp     = split(unique(unique(sample$Year)),grp.ind)
  
  x.design = (0:p)/p
  
  if(diff == F){
    help = function(h,...){
      w_h = locPolWeights(x = x.design, bw = h, deg = deg, xeval= x.design, kernel=EpaK)$locWeig
    
        max_diff = do.call(cbind,future_lapply(1:K, function(kk) {
        test     = which(sample$Year %in% grp[[kk]])
        train    = setdiff(1:n, test)
      
        test_grp     = as.vector(apply(sample[test,-1],2,mean, na.rm = T))
        pred         = w_h %*% as.vector(apply(sample[train,-1],2,mean, na.rm = T))
        max(abs(test_grp - pred))}))
    
      mean(max_diff)}
  
  mean_sup = sapply(h.seq, function(x){help(x,...)})
  return(h.seq[max(which(mean_sup == min(mean_sup)))])
  }else{help = function(h,...){
    
        arg = list(...)
        sample.d = arg[[1]]
        h.d      = arg[[2]]
    
        pd = dim(sample.d)[2] - 2
        nd = dim(sample.d)[1] 
    
        w.d = locPolWeights(x = (0:pd)/pd, bw = h.d, deg = deg, xeval= x.design, kernel=EpaK)$locWeig
        w_h = locPolWeights(x = x.design, bw = h, deg = deg, xeval= x.design, kernel=EpaK)$locWeig
    
        res = (-1)*sample + as.vector(w.d %*% apply(sample.d[,-1],2,mean, na.rm = T))
      
        max_diff = do.call(cbind,future_lapply(1:K, function(kk) {
            test     = which(sample$Year %in% grp[[kk]])
            train    = setdiff(1:n, test)
      
            test_grp     = as.vector(apply(res[test,-1],2,mean, na.rm = T))
            pred         = w_h %*% as.vector(apply(res[train,-1],2,mean, na.rm = T))
            max(abs(test_grp - pred))}))
    
      mean(max_diff)}
  
  mean_sup = sapply(h.seq, function(x){help(x,...)})
  return(h.seq[max(which(mean_sup == min(mean_sup)))])
  }
} 

bw_month = function(data.sparse,data.dense){
  
  bw_list = matrix(NA, nrow = 2, ncol = 12,dimnames = list(c("dense", "delta"), month.name)) |> as.data.frame()
  
  results = future_sapply(c(1,4,7,10), function(m){
    m.neighbor = c((m - 2) %% 12 + 1,m,m %% 12 + 1) 
    sample.dense  = data.dense %>% filter(MONTH %in% month.name[m.neighbor]) |> dplyr::select(1,4:dim(data.dense)[2])
    sample.dense  = sample.dense[rowSums(is.na(sample.dense)) == 0,]
    K.dense  = length(unique(sample.dense$Year))
    sample.sparse = data.sparse %>% filter(MONTH %in% month.name[m.neighbor])|> dplyr::select(1,4:dim(data.sparse)[2])
    sample.sparse = sample.sparse[rowSums(is.na(sample.sparse)) == 0,]
    K.sparse = length(unique(sample.sparse$Year))
    
    bw.dense      = k.fold.hc.cv(sample = sample.dense, h.seq = seq(round(0.05*(30/34), digits = 2), 0.12,0.002), deg = 2, K = K.dense)  
    bw.delta      = k.fold.hc.cv(sample = sample.sparse,h.seq = seq(round(0.095*(30/34), digits = 2),0.18,0.002), deg = 2, diff = T,K = K.sparse, sample.dense, bw.dense) 
    
    sample.dense  = data.dense %>% filter(MONTH %in% month.name[m]) |> dplyr::select(1,4:dim(data.dense)[2])
    sample.dense  = sample.dense[rowSums(is.na(sample.dense)) == 0,]
    sample.sparse = data.sparse %>% filter(MONTH %in% month.name[m])|> dplyr::select(1,4:dim(data.sparse)[2])
    sample.sparse = sample.sparse[rowSums(is.na(sample.sparse)) == 0,]
    print(c(dim(sample.dense)[1],dim(sample.sparse)[1]))
    c(bw.dense,bw.delta)
  }, future.seed = T)

  bw_list[1, c(12,1:11)] = rep(results[1,],each = 3)
  bw_list[2, c(12,1:11)] = rep(results[2,],each = 3)
  
  return(bw_list)
}

est.results = function(data.sparse,data.dense,bandwidth,from = 1,to = 12){
  
  p  = dim(data.sparse)[2]- 4
  pd = dim(data.dense)[2] - 4
  
  #eval.grid =  (12:(120-12))/120  # every 15 Minute from -21:00 till +03:00: 4*3 + (4*24+1) + 4*3 = 12 + 97 + 12
  eval.grid =  (20:(136-20))/136  # every 15 Minute from -19:00 till +05:00: 4*5 + (4*24+1) + 4*5 = 20 + 97 + 20
  
  eval.tibble = tibble(TIME = hms::as_hms(seq(from = as.POSIXct("1970-01-01 00:00:00"),to   = as.POSIXct("1970-01-01 23:45:00"),by   = "15 min")))
  eval.tibble = eval.tibble %>% add_row(TIME = as_hms("23:59:59"))
  
  result = future_sapply(from:to, function(m){
    
    bw.dense  = bandwidth[1,m]
    bw.delta  = bandwidth[2,m]
    
    sample.dense  = data.dense %>% filter(MONTH == month.name[m]) |> dplyr::select(4:dim(data.dense)[2])
    sample.dense  = sample.dense[rowSums(is.na(sample.dense)) == 0,]
    
    L.dense = tibble(TIME = colnames(sample.dense), MEAN = apply(sample.dense,2,mean,na.rm = T))
    
    L.dense.eval          = eval.tibble
    L.dense.eval$ESTIMATE = locPolSmootherC(x = (0:pd)/pd, y = L.dense$MEAN, xeval = eval.grid, bw = bw.dense, deg = 2, EpaK)$beta0
    L.dense.eval$MONTH    = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    sample.sparse  = data.sparse %>% filter(MONTH == month.name[m]) |> dplyr::select(4:dim(data.sparse)[2])
    sample.sparse  = sample.sparse[rowSums(is.na(sample.sparse)) == 0,]
    
    L.sparse = tibble(TIME  = colnames(sample.sparse), MEAN = apply(sample.sparse,2,mean,na.rm = T))
    dense.hat               = locPolSmootherC(x = (0:pd)/pd,y = L.dense$MEAN, xeval = (0:p)/p, bw = bw.dense,deg = 2, EpaK)$beta0
    res                     = dense.hat - L.sparse$MEAN 
    
    L.sparse.eval           = eval.tibble
    dense.hat.eval          = locPolSmootherC(x = (0:pd)/pd, y = L.dense$MEAN, xeval= eval.grid ,bw = bw.dense,deg = 2, EpaK)$beta0
    L.sparse.eval$ESTIMATE  = dense.hat.eval - locPolSmootherC(x = (0:p)/p, y = res, xeval = eval.grid, bw = bw.delta,deg = 2, EpaK)$beta0
    L.sparse.eval$MONTH     = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    L.delta.eval            = eval.tibble
    L.delta.eval$ESTIMATE   = locPolSmootherC(x = (0:p)/p, y = res,xeval = eval.grid ,bw = bw.delta,deg = 2, EpaK)$beta0 
    L.delta.eval$MONTH      = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    int.delta.eval          = eval.tibble
    int_delta_hat           = integrate(function(x)locPolSmootherC(x = (0:p)/p, y = res,xeval = x ,bw = bw.delta,deg = 2, EpaK)$beta0,lower = eval.grid[1],upper =  eval.grid[length(eval.grid)])$value/(eval.grid[length(eval.grid)]-eval.grid[1])  
 
    int.delta.eval$ESTIMATE = rep(int_delta_hat, times = length(eval.tibble))
    int.delta.eval$MONTH    = factor(rep(month.name[m], each = length(eval.grid)), levels = month.name)
    
    list(delta = L.delta.eval, mu_dense = L.dense.eval, mu_sparse = L.sparse.eval, delta_int = int.delta.eval) 
    
  }, future.seed = T)
  
  result.delta  = bind_rows(result[1,])
  result.dense  = bind_rows(result[2,])
  result.sparse = bind_rows(result[3,])
  result.int    = bind_rows(result[4,])
  
  
  return(list(delta = result.delta, sparse = result.sparse, dense = result.dense, delta_int = result.int) )
}




########################################################################################################################
################ Lagged kernel: Cross Validation #######################################################################
########################################################################################################################



cov.k.fold.hc.cv = function(sample, h.seq, deg = 1, K = 5, w.parallel = F){
  
  p = dim(sample)[2] - 1
  n = dim(sample)[1] 
  
  n.years = length(unique(sample$Year))
  grp.ind = rep(seq_len(K), times = c(rep(n.years %/% K, K - 1), n.years %/% K + n.years %% K))
  grp     = split(unique(unique(sample$Year)),grp.ind)
  
  help = function(h){
    w_h = local.polynomial.weights(p, h, p.eval = p, m = deg, parallel = w.parallel, parallel.environment = F)
    
    max_diff = do.call(cbind,future_lapply(1:K, function(kk) {
      test     = which(sample$Year %in% grp[[kk]])
      train    = setdiff(1:n, test)
      
      test_grp     = matrix(observation.transformation(sample[test,-1], grid.type = "full"), p, p)
      train_grp    = observation.transformation(sample[train,-1], grid.type = "less")
      pred         = eval_weights(w_h, train_grp)
      max(abs((test_grp - pred)[!as.logical(diag(p))]))}))
    
    mean(max_diff)}  
  mean_sup = mean_sup = sapply(h.seq, help)
  print("done")
  
  return(h.seq[max(which(mean_sup == min(mean_sup)))])
}

cov.bw.month = function(data, h.seq){
  
  bw_list = matrix(NA, nrow = 1 , ncol = 4) |> as.data.frame()
  results = sapply(c(1,4,7,10), function(m){
    
    m.neighbor = c((m - 2) %% 12 + 1,m,m %% 12 + 1) 
    data       = data %>% filter(MONTH %in% month.name[m.neighbor]) |> dplyr::select(1,4:dim(data)[2])
    data       = data[rowSums(is.na(data)) == 0,] 
    K.data     = length(unique(data$Year))
    bw         = cov.k.fold.hc.cv(data,h.seq = h.seq, w.parallel = T, deg = 1, K = K.data)
    bw
  })
  bw_list[1, c(12,1:11)] = rep(results,each = 3)
  dimnames(bw_list)  = list(paste0("lag ", 0), month.name)
  return(bw_list)
}


cov.weights = function(data.sparse,data.dense,bandwidth.s,bandwidth.d,from = 1, to = 12, p.eval , eval.type = "full", deg = 1, max.lag = 10){
  
  pd = dim(data.dense)[2] - 3
  p = dim(data.sparse)[2] - 3
  
  for(m in from:to){
    
    file.s = paste0("Application/Data examples/",data.example[[l]],"/weights/",month.name[m],"/sparse/",eval.type,"_w_lag0.rds")
    #file.d = paste0("Application/Data examples/",data.example[[l]],"/weights/",month.name[m],"/dense/", eval.type,"_w_lag0.rds")
    
    w  = local.polynomial.weights(p,  bandwidth.s[1,m], p.eval = p.eval, parallel = T, m = deg, del = 0, grid.type = "less", eval.type = eval.type, parallel.environment = F)
    #wd = local.polynomial.weights(pd, bandwidth.d[1,m], p.eval = p.eval, parallel = T, m = deg, del = 0, grid.type = "less", eval.type = eval.type, parallel.environment = F)    
    saveRDS(w, file = file.s)
    #saveRDS(wd, file = file.d)
    rm(w); gc()
    #rm(wd); gc()
    
    for(k in 1:max.lag){
      
     file.s = paste0("Application/Data examples/",data.example[[l]],"/weights/",month.name[m],"/sparse/",eval.type,"_w_lag",k,".rds")
     #file.d = paste0("Application/Data examples/",data.example[[l]],"/weights/",month.name[m],"/dense/", eval.type,"_w_lag",k,".rds")
    
     bw.lag   = round(bandwidth.s[1,m]*1.1^k, digits = 2)
     #bw.d.lag = round(bandwidth.d[1,m]*1.1^k, digits = 2)
      
      w.lag   = lag.local.polynomial.weights(p, bw.lag, p.eval = p.eval, parallel = T, m = deg, del = 0, grid.type = "lesseq", eval.type = eval.type, parallel.environment = F)
      #wd.lag  = lag.local.polynomial.weights(pd,bw.d.lag, p.eval = p.eval, parallel = T, m = deg, del = 0, grid.type = "lesseq", eval.type = eval.type, parallel.environment = F)
      saveRDS(w.lag, file = file.s)
      #saveRDS(wd.lag, file = file.d)
      
      print(paste0("done: lag ", k))
      
      rm(w.lag); gc()
      #rm(wd.lag); gc()
    } 
    print(paste0("done: ", month.name[m]))
  }
}  



### Integral delta Test

test.int = function(data.sparse, data.dense, int.est, var.s, var.d, from = 1, to = 12, test = "two-sided", alpha = 0.9){
  
  nd = unlist(lapply(from:to, function(m) {sample.dense  = data.dense %>% filter(MONTH %in% month.name[m]) |> dplyr::select(4:dim(data.dense)[2])
                                           sample.dense  = sample.dense[rowSums(is.na(sample.dense)) == 0,]
                                           return(dim(sample.dense)[1])}))
  n  = unlist(lapply(from:to, function(m) {sample.sparse = data.sparse %>% filter(MONTH %in% month.name[m])|> dplyr::select(4:dim(data.sparse)[2])
                                           sample.sparse = sample.sparse[rowSums(is.na(sample.sparse)) == 0,]
                                           return(dim(sample.sparse)[1])}))
  test.statistic = unlist(lapply(from:to, function(m) { sqrt(n[m]) * int.est[m] / sqrt(var.s[m] + n[m]/nd[m] * var.d[m])}))
  
  if(test == "two-sided"){  p.val = unlist(lapply(from:to, function(m) {2*pnorm(abs(test.statistic[m]),0,1,lower.tail = F) }))}
  if(test == "right-sided"){p.val = unlist(lapply(from:to, function(m) {  pnorm(test.statistic[m],0,1,lower.tail = F) }))}
  if(test == "left-sided"){ p.val = unlist(lapply(from:to, function(m) {  pnorm(test.statistic[m],0,1,lower.tail = T) }))}
  
  confInterval = lapply(from:to, function(m) {c(int.est[m]-sqrt(var.s[m] + n[m]/nd[m] * var.d[m])*qnorm(alpha) / sqrt(n[m]),int.est[m] + sqrt(var.s[m] + n[m]/nd[m] * var.d[m])*qnorm(alpha) / sqrt(n[m]))})
  names(confInterval) = month.name[from:to]
  return(list(pval = round(p.val, digits = 5), confInterval = confInterval))
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


q.month = function(data.sparse, data.dense, bandwidth, est, cov, from = 1, to = 12,alpha = 0.9, B = 1000, depend = F, int = F, constant = NA){
  
  res = sapply(from:to,  function(m){
    
    sample.dense  = data.dense %>% filter(MONTH %in% month.name[m]) |> dplyr::select(4:dim(data.dense)[2])
    sample.dense  = sample.dense[rowSums(is.na(sample.dense)) == 0,]
    sample.sparse = data.sparse %>% filter(MONTH %in% month.name[m])|> dplyr::select(4:dim(data.sparse)[2])
    sample.sparse = sample.sparse[rowSums(is.na(sample.sparse)) == 0,]
    
    eval.grid =  (20:(136-20))/136
    
    if(any(is.na(constant))){ constant = numeric(12)} #H0: delta = 0
    
    list = q.MB(sample.sparse, sample.dense, est$sparse%>%filter(MONTH %in% month.name[m]), est$dense%>%filter(MONTH %in% month.name[m]), unlist(cov[[m]]), eval.grid, bandwidth[,m], alpha = alpha, B = B, depend = depend, int = int, H0 = constant[m])
    
    print(paste0("done: month ", month.name[m], ": ", list$quantile ))
    
    if(int == F){obs = sqrt(dim(sample.sparse)[1]) * max(abs((est$delta%>%filter(MONTH %in% month.name[m]))$ESTIMATE -  constant[m] )/ sqrt(unlist(cov[[m]]))) }
    else{        obs = sqrt(dim(sample.sparse)[1]) * max(abs((est$delta%>%filter(MONTH %in% month.name[m]))$ESTIMATE - (est$delta_int%>%filter(MONTH %in% month.name[m]))$ESTIMATE) / sqrt(unlist(cov[[m]]))) }
    
    p.value = sum(list$sample >= obs)/length(list$sample)
    c(list$quantile, p.value)})
  
  colnames(res) = month.name[from:to]
  rownames(res) = c( paste0(alpha,"-quantile") , "p.val" )
  
  return(res)
}
  

q.MB = function(Ys, Yd, Ys.est,  Yd.est, cov, x.eval, bandwidth, alpha = 0.9, B = 1000, depend = F, int = F, H0 = 0){
  
  grid.s  = c((0:(dim(Ys)[2]-1))/(dim(Ys)[2]-1))[6:(dim(Ys)[2]-5)]
  grid.d  = c((0:(dim(Yd)[2]-1))/(dim(Yd)[2]-1))[31:(dim(Yd)[2]-30)]

  w.s = locPolWeights(x= grid.s, bw = bandwidth[2], deg = 2, xeval= x.eval, kernel=EpaK)$locWeig
  w.d = locPolWeights(x= grid.d, bw = bandwidth[1], deg = 2, xeval= grid.s, kernel=EpaK)$locWeig
  
  col.s.25h = !grepl("Before|After", colnames(Ys)) | colnames(Ys) == "After 00:00:00"
  col.d.25h = !grepl("Before|After", colnames(Yd)) | colnames(Yd) == "After 00:00:00"
  
  ls = list(sample.s = t(Ys[,col.s.25h]), weight.s = w.s, sample.s.mean = unlist(colMeans(Ys[,col.s.25h])), est.s = Ys.est$ESTIMATE)
  ld = list(sample.d = t(Yd[,col.d.25h]), weight.d = w.d, sample.d.mean = unlist(colMeans(Yd[,col.d.25h])), est.d = Yd.est$ESTIMATE) 
  
  if(int == F){
    ls = modifyList(ls, list(constant = H0))
    ld = modifyList(ld, list(constant = H0))
  }
  
  if(int == T){
    
    grid.s  = (0:(dim(Ys)[2]-1))/(dim(Ys)[2]-1)
    grid.d  = (0:(dim(Yd)[2]-1))/(dim(Yd)[2]-1)
    
    ls = modifyList(ls, list(int.s = integrate(function(x)locPolSmootherC(x = grid.s, y = colMeans(Ys), xeval = x ,bw = bandwidth[2],deg = 2, kernel = EpaK)$beta0,lower = x.eval[1],upper =  x.eval[length(x.eval)])$value/(x.eval[length(x.eval)]-x.eval[1])))
    ld = modifyList(ld, list(int.d = integrate(function(x)locPolSmootherC(x = grid.d, y = colMeans(Yd), xeval = x ,bw = bandwidth[1],deg = 2, kernel = EpaK)$beta0,lower = x.eval[1],upper =  x.eval[length(x.eval)])$value/(x.eval[length(x.eval)]-x.eval[1])))
    
    
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

  if(int == T){f1_int = sapply(1:n,  function(i){mean(locPolSmootherC((0:(p-1))/(p-1),   sample[, i], seq(0, 1, length.out = 1000), max(h,0.09), 2, EpaK)$beta0)})
               f2_int = sapply(1:nd, function(i){mean(locPolSmootherC((0:(pd-1))/(pd-1), sample_d[, i], seq(0, 1, length.out = 1000), h_d, 2, EpaK)$beta0)})}

  if(int == F){return(max(abs((1/sqrt(n-1)*( (weights %*% sample - T_val1) - int1)) %*% g_n  -
                     (n/(nd*sqrt(n-1))*( (weights %*% drop(weights_d %*% sample_d) - T_val2) - int2)) %*% g_nd)/sqrt(cov))) }
  
  else{        return(max(abs(1/sqrt(n-1)*(sweep(weights %*% sample, 2, f1_int, "-") - T_val1 + int1) %*% g_n  -
                     (n/(nd*sqrt(n-1))*(weights %*% drop(sweep(weights_d %*% sample_d, 2, f2_int, "-")) - T_val2 + int2)) %*% g_nd)/sqrt(cov)))}
}





CB = function(data.sparse, estimation, cov, q.list, center = F){
  
  new_est = estimation$delta
  if(center == T){new_est$ESTIMATE = estimation$delta$ESTIMATE  - estimation$delta_int$ESTIMATE}
  
  UP = c()
  LO = c()
  
  for(m in 1:12){sample.sparse = data.sparse %>% filter(MONTH %in% month.name[m])|> dplyr::select(4:dim(data.sparse)[2])
                    sample.sparse = sample.sparse[rowSums(is.na(sample.sparse)) == 0,]
                          
                    UP = c(UP,new_est |> filter(MONTH ==  month.name[m]) |> select(ESTIMATE)  + q.list[m] * sqrt(cov[[m]]/dim(sample.sparse)[1])) 
                    LO = c(LO,new_est |> filter(MONTH ==  month.name[m]) |> select(ESTIMATE)  - q.list[m] * sqrt(cov[[m]]/dim(sample.sparse)[1])) }
  
  
  UP = unlist(UP)
  UP = unname(UP)
  
  LO = unlist(LO)
  LO = unname(LO)
  
  new_est$UP = UP
  new_est$LO = LO
  return(new_est)
}


