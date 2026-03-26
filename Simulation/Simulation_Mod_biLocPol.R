#########################################################
#########################################################
###                                                   ###
###  Adaptation and expansion of R Package BiLocPol   ###
###                                                   ###
#########################################################
#########################################################


# Grid setting: t_j = (j -0.5)/p
# like in Berger and Holzmann (2025)


####################################
#### Calculation of the weights ####
####################################

#' Expanded and adapted local_polynomial_weights() to local.polynomial.weights()
#'   Features: - Grid.type -lesseq- is implemented
#'


local.polynomial.weights = function(p, h, p.eval, parallel = F, m = 1,
                                    del = 0, x.design.grid = NULL,
                                    grid.type = "less", eval.type = "full", parallel.environment = T, ...){
  
  observation.grid = function(p = NULL, x = NULL, comp = "less") {
    if(is.null(x) & is.null(p))
      stop("p or grid needs to be supplied")
    else if(is.null(x))
      x = (1:p - 0.5)/p
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
                       diagonal = matrix( (1:p.eval - 0.5)/p.eval , p.eval, 2)   )
  
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
  L = list(design = x.design.grid, x.design = (1:p - 0.5)/p , p = p, grid.type = grid.type,
           eval = x.eval.grid, x.eval = (1:p.eval - 0.5)/p.eval , p.eval = p.eval, eval.type = eval.type,
           bandwidth = h, m = m, del = del, weights = weights)
}



