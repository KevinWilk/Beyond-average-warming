#########################################################
#########################################################
###                                                   ###
###  Adaptation and expansion of R Package BiLocPol   ###
###                                                   ###
#########################################################
#########################################################




####################################
#### Calculation of the weights ####
####################################

#' Expanded and adapted local_polynomial_weights() to local.polynomial.weights()
#'   Features: - Grid.type -lesseq- is implemented
#'             - x.design apdapted to 0:(p-1)/(p-1) 
#'
#' Adapted observation_grid() to observation.grid()
#'   Feature: - grid adapted to 0:(p-1)/(p-1)

local.polynomial.weights = function(p, h, p.eval, parallel = F, m = 1,
                                    del = 0, x.design.grid = NULL,
                                    grid.type = "less", eval.type = "full", parallel.environment = T, ...){
  
  
  observation.grid = function(p = NULL, x = NULL, comp = "less") {
    if(is.null(x) & is.null(p))
      stop("p or grid needs to be supplied")
    else if(is.null(x))
      x = (0:(p-1))/(p-1)
    x.grid = expand.grid(x, x)
    
    if (comp == "less")
      b = x.grid[,1] < x.grid[,2]
    if (comp == "lesseq")
      b = x.grid[,1] <= x.grid[,2]
    if (comp == "gtr")
      b = x.grid[,1] > x.grid[,2]
    if (comp == "gtreq")
      b = x.grid[,1] >= x.grid[,2]
    if (comp == "without diagonal")
      b = !(x.grid[,1] == x.grid[,2])
    if (comp == "full")
      b = T
    x.grid[b, ]
  }
  
  
  if( !(grid.type %in% c("less", "lesseq", "without diagonal")) ){
    stop("grid type is not feasible - choose -less-, -lesseq- or -without diagonal-.")
  }
  
  
  if(!is.null(x.design.grid) & is.vector(x.design.grid)){
    x.design.grid = observation.grid(x = x.design.grid, comp = grid.type)
  }else{
    x.design.grid = observation.grid(p = p, comp = grid.type)
  }
  
  if( !(eval.type %in% c("full", "diagonal")) ){
    stop("evaluation is only possible for lower diagonal -full- or only on the -diagonal-.")
  }
  
  
  x.eval.grid = switch(eval.type,
                       full = observation.grid(p = p.eval, comp = "lesseq"),
                       diagonal = matrix( (0:(p.eval-1))/(0:(p.eval-1))/(p.eval-1), p.eval, 2)   )
  
  if (del > m) {
    m = 2; del = 2
    warning("del > m; calculations executed for del = m = 2.")
  }
  if (parallel) {
    if (parallel.environment) {
      cl = parallel::makeCluster(parallel::detectCores( ) - 1)
      future::plan(future::multisession)
    }
    w = future.apply::future_apply(x.eval.grid, 1, FUN = weights_point,
                                   x.design.grid = x.design.grid, h = h, m = m, del = del, ...,
                                   future.seed = T)
    if(parallel.environment){
      parallel::stopCluster(cl)
    }
  }else{
    w = apply(x.eval.grid, 1, weights_point,
              x.design.grid = x.design.grid,
              h = h, m = m, del = del, ...)
  }
  
  if(del == 0){
    weights = w
  }else{
    if (grid.type == "less") {k = 2} else {k = 1}
    weights = slice_matrix(w, p*(p-1)/k, switch(del, 3, 6))
  }
  rm(w)
  L = list(design = x.design.grid, x.design = (0:(p-1))/(p-1), p = p, grid.type = grid.type,
           eval = x.eval.grid, x.eval = (0:(p.eval-1))/(p.eval-1), p.eval = p.eval, eval.type = eval.type,
           bandwidth = h, m = m, del = del, weights = weights)
}



###################################
#### Evaluation of the weights ####
###################################

#' Extended eval_weights() to eval.weights()
#' Feature: - Evaluates the weights from local_polynomial_weights with |lag| > 0 according to Wilk/Holzmann (2026) 
#'  @param lag: lag = 0 refers to eval_weights()
#'  @param lag: |lag| > 0  refers to evaluating the weights of lagged kernels

eval.weights = function(W, Z,lag = 0){
  
  if(lag == 0){
    
    eval_weights(W, Z) # function of R Package BiLocPol 
    
  }else{
  
    if(W$del == 0){ 
      if (W$grid.type %in% c("lesseq","less")){
        if (W$eval.type == "full") {
          M1       = matrix(0, W$p.eval, W$p.eval)
          M2       = matrix(0, W$p.eval, W$p.eval)
          M.up     = upper.tri(M1, T)
          M1[M.up] = if(Z$lag >= 0) crossprod(W$weights, Z$Z1) else crossprod(W$weights, Z$Z2)
          M2[M.up] = if(Z$lag >= 0) crossprod(W$weights, Z$Z2) else crossprod(W$weights, Z$Z1)
          M        = M1 + t(M2)
          diag(M)  = diag(M) / 2
        } else {
          if(W$eval.type != "diagonal") {stop("eval.type needs to be diagonal or full")}
          M = (crossprod(W$weights, Z$Z1) + crossprod(W$weights, Z$Z2))/2
        }
      } else {
        if(W$grid.type != "full") {stop("grid type not implemented")}
        if(W$p*(W$p-1) != length(Z)) {stop("Z is not in correct representation")}
        else {
          if(W$eval.type != "diagonal") {stop("eval.type needs to be diagonal or full")}
          M = crossprod(W$weights, Z)
        }
      }
      if(Z$lag >= 0){return(M)}
      else{return(t(M))}
    }
    else{stop("grid type not implemented")}
  }
}


#################################
#### Computational functions ####
#################################

#' Adapted and extended observation_transformation() to observation.transformation()
#'
#'  Transforms observations according to Wilk/Holzmann (2026) based on Berger/Holzmann (2024)
#'
#'   Features: - returns vector with periodic transformed obseravtions
#'   
#'  @param periodic: periodic = F refers to no periodic transformed obseravtions
#'                   periodic = T refers to periodic transformed obseravtions
#'  @param m: periodic length

observation.transformation = function(Y, grid.type = "less", na.rm = F, periodic = F, m = NA, lag = 0){
  
  p = length(Y[1, ])
  n = length(Y[, 1])
  
  if(lag == 0){
  
    if(periodic == T){
      
      if(is.na(m)){m = p - 1}
      map = ((1:p) - 1) %% m + 1
      
      Y.means.base = vapply(1:m,
                              function(g) mean(unlist(Y[, which(map == g)]), 
                                na.rm = TRUE),numeric(1)
                           )
      Y.means = Y.means.base[map]
    }else{
      Y.means = colMeans(Y, na.rm = na.rm) 
    }
  
    Y.2 = apply(Y, 1, tcrossprod)
    Z = Y.2 - as.vector(tcrossprod(Y.means))
    Z.mean =  Z |> 
                apply(1, 
                    function(x){sum(x, na.rm = T)/(n - 1 - sum(is.na(x)))}
                     )
  
    if(periodic == T){
      
      col = rep(1:p, each = p)
      row = rep(1:p, times = p)
      group = (map[row] - 1) * m + map[col]
      Z.mean.base = tapply(Z.mean, group, mean, na.rm = TRUE)
      Z.mean = as.numeric(Z.mean.base[group])
      
    }
    if (Z.mean |> is.na() |> any()) { warning("There are NA values in the empirical covariance") }
    if (grid.type == "less") {
     M2 = as.vector(upper.tri(matrix(0, p, p)))  
    } else if (grid.type == "without diagonal") {
      M2 = as.vector(!diag(T, p, p))
    } else if (grid.type == "full") {
      M2 = T
    } else if (grid.type == "diagonal"){
      M2 = as.vector(diag(T, p, p))
    } else {stop("grid type not implemented")}
    return(Z.mean[M2])
  
  }else{
    
    if(periodic == T){
      if(is.na(m)){m = p - 1}
      map = ((1:p) - 1) %% m + 1
      Y.means.base = vapply(1:m,
                            function(g) mean(unlist(Y[, which(map == g)]), 
                                             na.rm = na.rm),numeric(1)
                            )
      
      Y.means = Y.means.base[map]
      
    }else{
      Y.means = colMeans(Y, na.rm = na.rm) 
    }
    
    
    if(lag >= 0){
      
      Y1.lag = do.call(cbind, lapply(lapply(
                                      1:(n - abs(lag)), 
                                      function(i) tcrossprod(unlist(Y[i,]), unlist(Y[i+abs(lag),]))
                                      ), as.vector))
      
      Y2.lag = do.call(cbind, lapply(lapply(
                                      (abs(lag)+1):n,   
                                      function(i) tcrossprod(unlist(Y[i,]), unlist(Y[i-abs(lag),]))
                                      ), as.vector))
      
      }else{        
      
      Y1.lag = do.call(cbind, lapply(lapply(
                                      (abs(lag)+1):n,   
                                      function(i) tcrossprod(unlist(Y[i,]), unlist(Y[i-abs(lag),]))
                                      ), as.vector))
      
      Y2.lag = do.call(cbind, lapply(lapply(
                                      1:(n - abs(lag)), 
                                      function(i) tcrossprod(unlist(Y[i,]), unlist(Y[i+abs(lag),]))
                                      ), as.vector))
      
      }
    
    Z1 = Y1.lag - as.vector(tcrossprod(Y.means))       
    Z2 = Y2.lag - as.vector(tcrossprod(Y.means))
    
    Z1.mean =  Z1 |>
      apply(1, function(x){
        sum(x, na.rm = na.rm)/(n - 1 - sum(is.na(x)))}
      )                                        
    Z2.mean =  Z2 |>
      apply(1, function(x){
        sum(x, na.rm = na.rm)/(n - 1 - sum(is.na(x)))}
      )     
    
    if(periodic == T){
      col = rep(1:p, each = p)
      row = rep(1:p, times = p)
      group = (map[row] - 1) * m + map[col]
      
      Z1.mean.base = tapply(Z1.mean, group, mean, na.rm = TRUE)
      Z1.mean =  as.numeric(Z1.mean.base[group])
      
      Z2.mean.base = tapply(Z2.mean, group, mean, na.rm = TRUE)
      Z2.mean =  as.numeric(Z2.mean.base[group])
    }
    
    if (Z1.mean |> is.na() |> any()) { warning("There are NA values in the empirical lagged covariance") }
    if (Z2.mean |> is.na() |> any()) { warning("There are NA values in the empirical lagged covariance") }
    
    if (grid.type == "lesseq") {
      M = as.vector(upper.tri(matrix(0, p, p)) | diag(T, p, p))  
      return(list(Z1 = Z1.mean[M],Z2 = Z2.mean[M], lag = lag))
    }else if (grid.type == "diagonal"){
      M = as.vector(diag(T, p, p))
      return((Z1.mean[M]+Z2.mean[M])/2)
    }else if (grid.type == "full") {
      
      M.diag = as.vector(diag(T, p, p))
      Z1.mean[M.diag] = Z1.mean[M.diag]/2
      Z2.mean[M.diag] = Z2.mean[M.diag]/2
      
      upr.Z1 = matrix(Z1.mean,p,p)
      upr.Z1[lower.tri(upr.Z1, diag = F)] = 0
      upr.Z2 = matrix(Z2.mean,p,p)
      upr.Z2[lower.tri(upr.Z2, diag = F)] = 0
      
      return(as.vector(upr.Z1+t(upr.Z2)))
    } else {stop("grid type not implemented")}
    
  }
}




