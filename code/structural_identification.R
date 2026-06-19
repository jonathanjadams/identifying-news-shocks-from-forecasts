##########################################################################################
#
# structural_identification.R
#
# Code to define structural identification in the "indentification" from news project
# Philip Barrett, Washington DC
# First version: 10apr2023
#
##########################################################################################

library(tidyverse)
library(vars)
library(expm)
library(lubridate)
library(magrittr)
library(combinat)
library(clue)
select <- dplyr::select
# Otherwise picks up the MASS usage

struct_id_rf <- function(l.rf.sim, sign.var=NA, n.match.pds=NULL, pt.est.only = FALSE ){
# Wrapper function for converting the bias-corrected reduced form estimates into
# structural estimates
  
  ### 1. Structural identification of the point estimate
  rf.pt <- l.rf.sim$rf
  sign.var.idx <- which( l.rf.sim$rf$args$var.names == sign.var )
  struct.pt <- structural_identification( rf.pt$l.B, rf.pt$Sigma, sign.var.idx )
  pt.est <- c( rf.pt, struct.pt )
  resid.pt <- structural_residuals( pt.est )
  pt.est <- c( rf.pt, struct.pt, resid.pt )
    # Add the structural residuals
  
  ### 2. Construct the loss function for ordering
  pt.est.irf <- irf_var_fcast( pt.est, n.pds = n.match.pds )
  loss.target <- pt.est.irf %>% 
    arrange( dimension, type, pd ) %>% 
    select(dimension, type, pd, contains('x.')) %>% 
    gather( variable, value, -dimension, -type, -pd ) %>%
    spread( dimension, value ) %>%
    select( -type, -pd, -variable ) %>%
    as.matrix
      # The IRF target
  if( pt.est.only ) return( list( pt=pt.est, irf=pt.est.irf ) )
  
  ### 3. Apply to the bootstrap
  loss.fn <- function( cand.irf ){
    err <- ( cand.irf - loss.target )
    return( err^2 %>% colSums() )
    # return( err^2 %>% sum() )
  }
  n.boot <- l.rf.sim$boot %>% length()
  message('Structural identification of the bootstrap')
  l.struct.id <- lapply( 1:n.boot, function(i){
    if(i %% 50 == 0) message('  n=', i)
    out <- structural_identification( l.rf.sim$boot[[i]]$l.B, l.rf.sim$boot[[i]]$Sigma,
                                      ordering.loss.fn=loss.fn, n.match.pds=n.match.pds )
    return(out)
  } )
  
  ### 4. Form up the results and return
  l.out <- list( pt=pt.est, boot=lapply( 1:n.boot, function(i) c( l.rf.sim$boot[[i]], l.struct.id[[i]] ) ) )
  return(l.out)
  
}

structural_identification <- function(l.B.j, m.Sigma, sign.var.idx=NA, 
                                      ordering.loss.fn=NULL, n.match.pds=NULL,
                                      full.rf=NA ){
  # Structural identification of the model from the reduced form parameters
  
  B.1 <- l.B.j[[1]]
  # Extract the first element
  n.var <- B.1 %>% nrow
  # Useful
  m.Sigma.11 <- m.Sigma[1:n.var,1:n.var]
  m.Sigma.21 <- m.Sigma[n.var+1:n.var,1:n.var]
  m.Sigma.22 <- m.Sigma[n.var+1:n.var,n.var+1:n.var]
  # The parts of Sigma
  m.phi <- m.Sigma.11 - B.1 %*% m.Sigma.21 - t(m.Sigma.21) %*% t(B.1) + 
    B.1 %*% m.Sigma.22 %*% t(B.1)
  m.phi.temp <- m.Sigma.21 - m.Sigma.22 %*% t(B.1)
  m.psi <- m.Sigma.22 - m.phi.temp %*% solve(m.phi) %*% t(m.phi.temp)
      # The matrices phi and psi
  phi.psi.svd <- svd( m.phi + m.psi )
  m.U <- phi.psi.svd$u
  m.lambda <- phi.psi.svd$d %>% sqrt %>% diag
  # The SVD of the sum of phi and psi
  svd.err.1 <- (m.phi + m.psi - m.U %*% m.lambda^2 %*% t(m.U)) %>% abs %>% max
  if(svd.err.1 > 1e-08) message('Warning: SVD decomposition error #1=', svd.err.1)
  # The error on the first SVD.  Fails if m.phi+m.psi is not PD
  m.lambda.inv <- (1/phi.psi.svd$d) %>% sqrt %>% diag
  m.kappa <- m.lambda.inv %*% t(m.U) %*% m.phi %*% m.U %*% m.lambda.inv
  kappa.svd <- m.kappa %>% svd
  m.V.prime <- kappa.svd$u
  m.D.v.2 <- kappa.svd$d %>% diag
  m.D.u.2 <- diag(n.var) - m.D.v.2
  # The variance matrices
  svd.err.2 <- (m.kappa - m.V.prime %*% m.D.v.2 %*% t(m.V.prime)) %>% abs %>% max
  if(svd.err.2 > 1e-08) message('Warning: SVD decomposition error #2=', svd.err.2)
  # The error on the second SVD
  
  # browser()
  m.A <- m.U %*% m.lambda %*% m.V.prime
  # m.A.sign.correction <- m.A %>% diag %>% sign %>% diag
  #   # Sign correction via positive diagonals of A
  if(!is.na(sign.var.idx)){
    m.A.sign.correction <- m.A[sign.var.idx,] %>% sign %>% diag
    # Sign correction via positive shocks to one variable
    m.A <- m.A %*% m.A.sign.correction
    # Imposes the sign correction
  }
  
  m.C <- m.phi.temp %*% solve( m.D.v.2 %*% t(m.A) )
  # The structural shock matrices
  
  # Reorder the shocks if there is a target ordering function
  # browser()
  if( !is.null(ordering.loss.fn) ){
    # browser()
    shk.order <- shk_sign_order( ordering.loss.fn, m.A, m.C, m.D.v.2, m.D.u.2, l.B.j, n.match.pds )
    m.A <- (m.A %*% shk.order$sign)[,shk.order$order]
    m.C <- (m.C %*% shk.order$sign )[,shk.order$order]
    m.D.v.2 <- m.D.v.2[shk.order$order,shk.order$order]
    m.D.u.2 <- m.D.u.2[shk.order$order,shk.order$order]
  }
  
  # Decomposition checks
  A.check.1 <- ( ( m.A %*% m.D.v.2 %*% t(m.A) ) - m.phi) %>%
    abs %>% max
  A.check.2 <- ( ( m.A %*% m.D.u.2 %*% t(m.A) ) - m.psi) %>%
    abs %>% max
  A.check.3 <- ( ( m.A %*% t(m.A) ) - m.phi - m.psi ) %>%
    abs %>% max
  
  out <- list( A=m.A, C=m.C, D.v.2=m.D.v.2, D.u.2=m.D.u.2,
               A.checks=c(A.check.1, A.check.2, A.check.3) )
  
  return( out ) 
}


shk_sign_order <- function( fn, m.A, m.C, m.D.v.2, m.D.u.2, l.B.j, n.match.pds ){
# Uses a loss function to evaluate alternate orderings and signs of the shocks
  
  # browser()
  n.shk <- m.A %>% ncol
  n.signs <- 2^n.shk
  m.signs <- (-1)^sapply( 1:n.signs, function(i) intToBits(i) %>% as.numeric() )[1:n.shk,]
  # The matrix of possible signs
  ##### TO ADD HERE: COMPUTE THE IRFS AND COMPARE ACROSS ORDERINGS AND SIGNAGES ####
  this.est <- list( A=m.A, l.B=l.B.j, C=m.C, D.v.2=m.D.v.2, D.u.2=m.D.u.2 )
  this.irf <- irf_var_fcast( this.est, n.pds = n.match.pds )
  this.irf.vals <- this.irf %>%
    arrange( dimension, type, pd ) %>% 
    select(dimension, type, pd, contains('x.')) %>% 
    gather( variable, value, -dimension, -type, -pd ) %>%
    spread( dimension, value ) %>%
    select( -type, -pd, -variable ) %>%
    as.matrix
  # browser()
  m.orderings <- apply( m.signs, 2, 
                        function(x) shk_order( fn, this.irf.vals * (rep(1,nrow(this.irf.vals)) %*% t(x) ) ) )
  # Compute the loss-minimizing orderings conditional on each signing of all the shocks
  # browser()
  
  v.loss <- sapply( 1:n.signs, 
                    function(i){
                      this.sign <- m.signs[,i]
                      this.order <- m.orderings[,i]
                      this.irf.sign <- this.irf.vals * (rep(1,nrow(this.irf.vals)) %*% t(this.sign) )
                      fn( this.irf.sign[,this.order] ) %>% sum()
                      # NB: l.B.j is just propagation so does not need permuting
                    } )
  min.idx <- which.min(v.loss)
    # The loss-minimizing sign index
  return( list( sign=m.signs[,min.idx] %>% diag, order=m.orderings[ , min.idx ] ) )
}


shk_order <- function( fn, this.irf.vals.in ){
  # Uses a loss function to evaluate alternate orderings of the shocks.
  
  n.shk <- this.irf.vals.in %>% ncol
  m.cost <- apply( this.irf.vals.in, 2, fn )
  m.opt.order <- solve_LSAP(m.cost) %>% c
  return(m.opt.order)
  # fn.old <- function(x) sum(fn(x))
  # l.orderings <- permn(1:n.shk)
  # v.loss <- sapply( l.orderings, function(x) fn.old( this.irf.vals.in[,x] ) )
  # return( l.orderings[[ which.min(v.loss) ]] )
    # Commented version does this by hand.
}

structural_residuals <- function( l.rf ){
# Computes the structural residuals
  
  # Extract elements of the reduced form
  B.1 <- l.rf$l.B[[1]]
  m.A <- l.rf$A
  m.C <- l.rf$C
  big.resid <- l.rf$big.resid
  n.var <- l.rf$n.var
  p <- l.rf$lags
  big.B.x <- l.rf$big.B.x
  
  # Compute the structural residuals
  big.A <- rbind(
    cbind( B.1 %*% m.C + m.A, B.1 %*% m.A ),
    cbind( m.C, m.A )
  )
  u.v.fitted <- solve( big.A, t(big.resid) ) %>% t()
    # Inverting equation (4)
  v.fitted <- u.v.fitted[,1:n.var]
  u.fitted <- u.v.fitted[,n.var+1:n.var]
    # The fitted values
  epsilon.fitted <- rbind( NA, u.fitted[-1,] + v.fitted[-nrow(v.fitted),] )
  
  # Estimation checks
  m.x <- l.rf$data[,n.var+1:n.var]
    # The non-forecast data
  big.X <- lapply( 0:(p-1), function(i) m.x[(p:nrow(m.x))-i,] ) %>%
    reduce(cbind)
  big.X <- big.X[-nrow(big.X),]
  if(p==1){
    resid.check <- (m.x[-(1:(1+p)),] - 
                      (t(big.B.x[1:n.var,] %*% t(big.X)) + t(m.A %*% t(epsilon.fitted)) + t(m.C %*% t(v.fitted)))[-(1:p),] ) %>% 
      abs %>% max
  }else{
    resid.check <- (m.x[-(1:p),] - 
                      (t(big.B.x[1:n.var,] %*% t(big.X)) + t(m.A %*% t(epsilon.fitted[-1,])) + t(m.C %*% t(v.fitted[-1,]))) ) %>% 
      abs %>% max
  } 
  
  # The output
  df.resid <- data.frame( date=( l.rf$df %>% select(date))[-(1:(p-1)),], 
                          v.fitted %>% set_colnames(paste0('vhat.', 1:n.var)),
                          u.fitted %>% set_colnames(paste0('uhat.', 1:n.var)),
                          epsilon.fitted %>% set_colnames(paste0('epshat.', 1:n.var))) %>%
    as_tibble()
  out <- list( resid.check=resid.check, df.resid=df.resid )
  
  return(out)
    
}


var.decomp <- function( l.est, n.pds=20 ){
  # The variance decomposition
  
  l.var.ratio <- l.var.y <- l.var.decomp <- l.covar.decomp <- list()
  # l.var.decomp.ratio.sumy <- l.var.decomp.ratio <- 
  n.vars <- l.est$n.var
  # Number of variables
  n.shks <- l.est$A %>% ncol()
  p <- l.est$lags
  # The number of lags
  B <- l.est$big.B.x
  # Extract the VAR elements
  var.u <- l.est$D.u.2 %>% diag
  var.v <- l.est$D.v.2 %>% diag
  
  A.sq.extended <- lapply( 1:n.shks, function(x){ 
    out <- 0*diag( n.vars * p )
    out[1:n.vars,1:n.vars] <- l.est$A[,x] %*% t(l.est$A[,x])
    return(out)
  } )
  C.sq.extended <- lapply( 1:n.shks, function(x){ 
    out <- 0*diag( n.vars * p )
    out[1:n.vars,1:n.vars] <- l.est$C[,x] %*% t(l.est$C[,x])
    return(out)
  } )
  A.C.prod.extended <- lapply( 1:n.shks, function(x){ 
    out <- 0*diag( n.vars * p )
    out[1:n.vars,1:n.vars] <- l.est$A[,x] %*% t(l.est$C[,x])
    return(out)
  } )
  
  # browser()
  cov.names <- sapply( 1:p, function(i) paste0( l.est$args$var.names, '.l', i-1 ) ) %>% c
  news.increment <- lapply( 1:n.shks, function(x) (var.v[x] * C.sq.extended[[x]] ) %>%
                              set_rownames(cov.names) %>% set_colnames(cov.names) )
  for( i in 1:n.pds ){
    if( i==1 ){
      l.covar.decomp[[i]] <- list( contemp = lapply( 1:n.shks, function(x) (var.u[x] * ( A.sq.extended[[x]] ) ) %>%
                                                       set_rownames(cov.names) %>% set_colnames(cov.names) ),
                                   news.rolling = lapply( 1:n.shks, function(x) 0 * C.sq.extended[[x]] ),
                                   news = news.increment )
      l.covar.decomp[[i]]$epsilon <- lapply( 1:n.shks, function(x) l.covar.decomp[[i]]$contemp[[x]] + 
                                               l.covar.decomp[[i]]$news[[x]] )
    }else{
      l.covar.decomp[[i]] <- list( contemp = lapply( 1:n.shks, function(x) (B %*% l.covar.decomp[[i-1]]$contemp[[x]] %*% t(B)) + 
                                                       var.u[x] * A.sq.extended[[x]] ),
                                   news.rolling=lapply( 1:n.shks, function(x) (B %*% l.covar.decomp[[i-1]]$news.rolling[[x]] %*% t(B)) +
                                                          var.v[x] * (A.sq.extended[[x]] + B %*% C.sq.extended[[x]] %*% t(B) + 
                                                                        A.C.prod.extended[[x]] %*% t(B) + B %*% t(A.C.prod.extended[[x]]) ) ) )
      l.covar.decomp[[i]]$news <- lapply( 1:n.shks, function(x) l.covar.decomp[[i]]$news.rolling[[x]] + 
                                            C.sq.extended[[x]] * var.v[x] )
      l.covar.decomp[[i]]$epsilon <- lapply( 1:n.shks, function(x) l.covar.decomp[[i]]$contemp[[x]] + l.covar.decomp[[i]]$news[[x]] )
    }
    # big.sigma.denominator <- ( sapply(l.var.decomp[[i]]$contemp, diag) %>% rowSums() + 
    #                              sapply(l.var.decomp[[i]]$news, diag) %>% rowSums() )
    #   # The total variance of each variable at horizon i
    # l.var.decomp.ratio[[i]] <- lapply( c(contemp='contemp', news='news'), function(nn) sapply( 1:n.vars, 
    #                               function(x) diag(l.var.decomp[[i]][[nn]][[x]]) / big.sigma.denominator ) )
    # l.var.decomp.ratio.sumy[[i]] <- sapply( l.var.decomp.ratio[[i]], rowSums )
    l.var.decomp[[i]] <- lapply( l.covar.decomp[[i]], function(x) sapply(x, diag)[1:n.vars,] )
    l.var.y[[i]] <- lapply( l.var.decomp[[i]], rowSums )
    # l.var.shk[[i]] <- l.var.decomp[[i]]$epsilon %>% colSums()
    l.var.ratio[[i]] <- lapply( c('contemp', 'news', 'epsilon') %>% set_names(.,.), 
                                function(x) l.var.decomp[[i]][[x]] / ( l.var.y[[i]][[x]] %*% t(rep(1, n.shks)) ) ) 
  }
  
  df.var.decomp.ratio <- suppressMessages(
    lapply( 1:n.pds, function(i) lapply( c(contemp='contemp', news='news', 'epsilon'),
                                         function(nn) l.var.ratio[[i]][[nn]] %>%
                                           as.data.frame %>% rownames_to_column('outcome') %>%
                                           gather(dimension, variance.ratio, -outcome) %>%
                                           mutate(pd=i, type=case_when(nn=='contemp'~'u',
                                                                       nn=='news'~'v',
                                                                       nn=='epsilon'~'epsilon'),
                                                  dimension=dimension %>% gsub('V','',.) %>% as.numeric()) ) %>%
              reduce(full_join) ) %>%
      reduce(full_join) )
  
  return( list( var.decomp=l.var.decomp, var.decomp.ratio=l.var.ratio, covar.decomp=l.covar.decomp,
                df=df.var.decomp.ratio ) )
  
}

