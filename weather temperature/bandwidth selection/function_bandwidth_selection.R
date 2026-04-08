

########################################################################################################################
################ Mean and Difference function: Cross Validation ########################################################
########################################################################################################################

# K-Fold hv-block cross validation for univariate local polynomial estimator

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


# runs k.fold.hc.cv() on densely (2000-2025) observed data sets and sparsly (1952-1972) observed residuals for each season 
# with K = length of time period (1952-1972: 21 years, 2000-2025: 26 years)

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




########################################################################################################################
################ Lagged kernel: Cross Validation #######################################################################
########################################################################################################################

# K-Fold hv-block cross validation for bivariate local polynomial estimator (lag 0)

cov.k.fold.hc.cv = function(sample, h.seq, deg = 1, K = 5, w.parallel = F, period.n = NA){
  
  p = dim(sample)[2] - 1
  n = dim(sample)[1] 
  
  n.years = length(unique(sample$Year))
  grp.ind = rep(seq_len(K), times = c(rep(n.years %/% K, K - 1), n.years %/% K + n.years %% K))
  grp     = split(unique(unique(sample$Year)),grp.ind)
  
  help = function(h){
    w_h = local_polynomial_weights(p, h, p.eval = p, m = deg, parallel = w.parallel, parallel.environment = F)

    max_diff = do.call(cbind,future_lapply(1:K, function(kk) {
      test     = which(sample$Year %in% grp[[kk]])
      train    = setdiff(1:n, test)
      
      test_grp     = matrix(observation.transformation(sample[test,-1], grid.type = "full", periodic = T, m = period.n), p, p)
      train_grp    = observation.transformation(sample[train,-1], grid.type = "less", periodic = T, m = period.n)
      pred         = eval_weights(w_h, train_grp)
      max(abs((test_grp - pred)[!as.logical(diag(p))]))}))
    
    mean(max_diff)}  
  mean_sup = mean_sup = sapply(h.seq, help)
  print("done")
    
  return(h.seq[max(which(mean_sup == min(mean_sup)))])
}



# runs cov.k.fold.hc.cv() for each season with K = length of time period (1952-1972: 21 years, 2000-2025: 26 years)

cov.bw.month = function(data, h.seq, period.n = NA){
  
  bw_list = matrix(NA, nrow = 1 , ncol = 4) |> as.data.frame()
  results = sapply(c(1,4,7,10), function(m){     # for each season
    
    m.neighbor = c((m - 2) %% 12 + 1,m,m %% 12 + 1) 
    data       = data %>% filter(MONTH %in% month.name[m.neighbor]) |> dplyr::select(1,4:dim(data)[2])
    data       = data[rowSums(is.na(data)) == 0,] 
    K.data     = length(unique(data$Year))
    bw         = cov.k.fold.hc.cv(data,h.seq = h.seq, w.parallel = T, deg = 1, K = K.data, period.n = period.n)
    bw
  })
  bw_list[1, c(12,1:11)] = rep(results,each = 3)
  dimnames(bw_list)  = list(paste0("lag ", 0), month.name)
  return(bw_list)
}



