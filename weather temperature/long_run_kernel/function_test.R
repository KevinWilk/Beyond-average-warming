inference.autocovariance.test = function(sample, m, alpha, max.lag = 1, procedure = "random"){ 
  
  q_weighted_chisq_imhof = function(lam, alpha = 0.9){
    lam = as.numeric(lam)
    
    tailp = function(x) imhof(x, lam)$Qq   # upper tail P(S >= x)
    
    mu  = sum(lam)
    sig = sqrt(2 * sum(lam^2))
    lower = 0
    upper = max(mu + 20 * sig, 1)  
    
    while (tailp(upper) > alpha) upper = upper * 2
    
    uniroot(function(x) tailp(x) - alpha,
            lower = lower, upper = upper, tol = 1e-8)$root
  }
  
  Xi_hat_MC <- function(lam, alpha = 0.9, B = 100000){
    lam = lam[lam > 0]
    Z = matrix(rnorm(B * length(lam)), B)
    S = Z^2 %*% lam
    as.numeric(quantile(S, probs = alpha))
  }
  
  res.q = future_lapply(1:max.lag,function(h) {
    
    sample  = sample %>% filter(MONTH %in% month.name[m])|> dplyr::select(c(1,3:dim(sample)[2]))
    sample  = sample[rowSums(is.na(sample)) == 0,]
    n.days  = length(unique(sample$DAY))
    n.Years = length(unique(sample$Year))
    Y = do.call(rbind,lapply(1:n.days,  function(kk){apply(sample %>% filter(DAY == unique(sample$DAY)[kk])|> dplyr::select(-c(1,2)),2,mean)}))
    
    p    = ncol(Y) 
    n.T  = nrow(Y)
    
    barY = colMeans(Y)
    
    D    = sweep(Y, 2, barY, "-")
    A    = D[1:(n.T-h), , drop=FALSE]            
    B    = D[(1+h):n.T, , drop=FALSE]  
    
    c.hat = array(0, dim = c(p, p, p, p))
    for(r in 1:(n.T-h)){M    = outer(A[r, ], B[r, ])        
    c.hat = c.hat + array( tcrossprod(as.vector(M)), dim = c(p, p, p, p) )}
    c.hat = c.hat/(p^2 * n.T)
    C.hat = as.tensor(c.hat)  
    sv   = svd.tensor(C.hat,i = c(1, 2),j = c(3, 4))
    lambda_hat = as.vector(sv$d)
    
    lam = lambda_hat[lambda_hat > 1e-10]
    
    Xi_hat_safe = function(lam, alpha){tryCatch(q_weighted_chisq_imhof(lam, alpha = 1-alpha),
                                                error = function(e) Xi_hat_MC(lam, alpha = alpha))}
    
    if(procedure == "-MC"){   Xi_hat = Xi_hat_MC(lam, alpha = alpha)}
    if(procedure == "-imhof"){Xi_hat = q_weighted_chisq_imhof(lam, alpha = alpha)}
    else{                     Xi_hat = Xi_hat_safe(lam, alpha = alpha)}
    
    q = sqrt(Xi_hat/n.T) / mean(diag(crossprod(D) / (n.T-1))) 
    q                                         
  }, future.seed = TRUE)
  return(unlist(res.q))
}






