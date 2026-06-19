##########################################################################################
#
# var_sim_fns.R
#
# Code define the simulation functions for the "Identification from 
# Rational Expectations" project
#
# Philip Barrett, Washington DC
# First version: 26jun2022
#
##########################################################################################

#### Libraries
library(vars)
library(tsDyn)
library(tidyverse)
library(mvtnorm)
library(expm)

## Simulate the VAR
sim_var_fcast <- function( x.0, l.B.j, A, C, D.u.2, D.v.2, n.pds, u=NULL, v=NULL ){
# Simulate the VAR and return the rational forecasts as well
  
  n.x <- x.0 %>% ncol
  n.lags <- l.B.j %>% length() # 
  n.lags.0 <- x.0 %>% nrow
  if( n.lags != n.lags.0 ) stop('Dimension of x.0 and j.B.j do not match')
    # Dimensions
  m.B.1 <- l.B.j %>% reduce(cbind)
  lagmat <- if( n.lags > 1 ) cbind( diag( n.x*(n.lags-1) ), 
                                    matrix(0, n.x*(n.lags-1), n.x ) ) else NULL
  m.B <- rbind( m.B.1, lagmat )
    # The reduced form VAR matrix
  m.A <- rbind( A, matrix(0, nrow=n.x*(n.lags-1), ncol=n.x) )
  m.C <- rbind( C, matrix(0, nrow=n.x*(n.lags-1), ncol=n.x) )
    # The extended versions of A and C for the stacked VAR
  if( is.null(u) ) u <- mvrnorm( n.pds, rep(0,n.x), D.u.2 )
  if( is.null(v) ) v <- mvrnorm( 1+n.pds, rep(0,n.x), D.v.2 )
    # The primitive innovations
  epsilon <- u + v[1:n.pds,]
    # The reduced form innovations
  # browser()
  x <- matrix( NA, n.pds+1, n.x*n.lags )
  x[1,] <- x.0 %>% t %>% c
    # Fill in the top row
  fcast <- matrix( NA, n.pds+1, n.x*n.lags )
    # The matrix of one-period forecasts
  
  
  for( this.pd in 1:n.pds ){
    fcast[this.pd,] <- m.B %*% c(x[this.pd,]) + 
      m.A %*% c(v[this.pd,])
    x[this.pd+1,] <- m.B %*% c(x[this.pd,]) +
      m.A %*% c(epsilon[this.pd,]) + m.C %*% c(v[this.pd+1,])
        # NB. Are computing the x outcome in the *next* period, so this is correct???
  }
  fcast[n.pds+1,] <- m.B %*% c(x[n.pds+1,]) + 
    m.A %*% c(v[n.pds+1,])  
  
  df <- data.frame( pd=0:n.pds, x=x[,1:n.x], f=fcast[,1:n.x], 
                    v=v, u=rbind(NA,u), epsilon=rbind(NA,epsilon) ) %>%
    as_tibble()
  return(df)
}

## Simulate the VARMA
sim_varma_fcast <- function( x.0, l.B.j, A, C, G, H, D.u.2, D.v.2, n.pds, u=NULL, v=NULL ){
# Simulate the VARMA and return the rational forecasts
  
  # Create the state variables
  df.var <- sim_var_fcast( x.0, l.B.j, A, C, D.u.2, D.v.2, n.pds, u, v )
  
  # Create the observations
  df.out <- state_to_varma_fcast( df.var, G, H )
  return(df.out)
  
}

state_to_varma_fcast <- function( df.var, G, H ){
  # Simulate the VARMA and return the rational forecasts
  
  # Extract the relevant parts as matrices
  m.x <- df.var %>% select( contains('x') ) %>% as.matrix()
  m.f.x <- df.var %>% select( contains('f') ) %>% as.matrix()
  m.u <- df.var %>% select( contains('u') ) %>% as.matrix()
  m.v <- df.var %>% select( contains('v') ) %>% as.matrix()
  m.v.l1 <- rbind( NA, m.v[-nrow(m.v),])
  
  # Create the observations and the forecast
  m.y <- ( m.x + (m.u + m.v.l1) %*% t(G) + m.v %*% t(H) ) %>%
    set_colnames( paste0( 'y.', 1:ncol(m.x) ) )
  m.f.y <- ( m.f.x + m.v %*% t(G) ) %>%
    set_colnames( paste0( 'f.y.', 1:ncol(m.x) ) )
  
  # Amalgamate and return
  df.out <- cbind( df.var, m.y, m.f.y ) %>%
    as_tibble()
  return(df.out)
  
}

irf_var_fcast <- function( l.est, var.names=NULL, n.pds=20, unit='sd', shk.size=NA ){
# Computes the impulse responses
  
  # Unpack the estimated VAR
  A <- l.est$A
  l.B.j <- l.est$l.B
  C <- l.est$C
  sd.u <- l.est$D.u.2 %>% diag %>% sqrt
  sd.v <- l.est$D.v.2 %>% diag %>% sqrt
  
  # Create inputs to the simulation
  n.vars <- C %>% nrow()
  n.lags <- l.B.j %>% length
  x.0 <- matrix(0, n.lags, n.vars)
  u.0 <- matrix( 0, n.pds, n.vars )
  v.0 <- matrix( 0, n.pds+1, n.vars )
  
  # Create the IRFs
  irf.fn <- function(i){
    this.u <- u.0 ; this.v <- v.0
    this.u[1,i] <- if(unit=='sd') sd.u[i] else if(is.na(shk.size)) 1 else shk.size$u[i]
    this.v[2,i] <- if(unit=='sd') sd.v[i] else if(is.na(shk.size)) 1 else shk.size$v[i] 
      # Set up the shocks
    dr.this.irf <- list( u=sim_var_fcast( x.0, l.B.j, A, C, NULL, NULL, n.pds, this.u, v.0 ) %>%
                           mutate(type='u'),
                      v=sim_var_fcast( x.0, l.B.j, A, C, NULL, NULL, n.pds, u.0, this.v ) %>%
                        mutate(type='v') ) %>%
      reduce(rbind) %>%
      mutate(dimension=i)
    return(dr.this.irf)
  }
  fcast.names <- if(is.null(var.names)) paste0('f.', 1:n.vars) else paste0( var.names, '_f' )
  if(is.null(var.names)) var.names <- paste0('x.', 1:n.vars)
  df.irfs <- suppressMessages( lapply( 1:n.vars, irf.fn ) %>% 
                                 reduce(rbind) %>%
                                 { if(n.vars>1) 
                                   {.} %>% rename_at(vars(c(paste0('x.', 1:n.vars), paste0('f.', 1:n.vars))),
                                           ~ c(var.names, fcast.names) ) else {.} }
                               )
    # Rstudio thinks that this is wrong but it is not
  return(df.irfs)
}

irf_var_epsilon <- function(df.irfs){
# Computes the non-causal epsilon shocks
  
  # browser()
  if('pc' %in% names(df.irfs) ) df.irfs <- df.irfs %>% select(-pc)
  
  df.wts <- df.irfs %>% 
    # filter(pc=='pt.est') %>% 
    select( type, pd, dimension, contains(c('v.', 'u.'))) %>% 
    filter( pd==1 ) %>%
    gather( variable, value, -type, -dimension, -pd ) %>% 
    filter( value > 0 ) %>%
    rename( shk.wt=value ) %>%
    select( type, dimension, shk.wt )
      # The weights on the types of shock. These are already 1SD shocks. So need
      # to multiply the IRFs by the magnitude of the shock SD once more. This
      # gives variance-weighted shocks
  
  # browser()
  df.irfs.adj <- df.irfs %>%
    mutate( pd=ifelse(type=='v', pd-1, pd) ) %>%
    # Change the timing so the epsilon shock occurs simultaneously across u
    # and v.
    filter( !(pd %in% range(pd)) ) %>%
    left_join(df.wts, by=c('type', 'dimension')) 
  
  df.out <- df.irfs.adj %>%
    gather( variable, value, -pd, -type, -dimension, -shk.wt, -shk.name, -irf.sign ) %>%
      # Drop the final period and period -1 due to the shift in timing.
    mutate( eps.irf=shk.wt*value ) %>%
    select( -value ) %>%
    full_join( {.} %>% group_by( dimension, pd, variable, shk.name, irf.sign ) %>%
                 summarise( eps.irf=sum(eps.irf, na.rm=T) ) %>%
                    # Add across u and v 
                 mutate( type='epsilon' ) ) %>%
    spread( variable, eps.irf ) %>%
    arrange( dimension, pd ) 
      # The impulse responses
  
  return( df.out )
}


# irf_var_rf <- function( l.est, var.names=NULL, n.pds=20 ){
# # Computes the impulse responses
#   
#   # Unpack the estimated VAR
#   A <- l.est$A
#   l.B.j <- l.est$l.B
#   C <- l.est$C
#   
#   # Create inputs to the simulation
#   n.vars <- C %>% nrow()
#   n.lags <- l.B.j %>% length
#   x.0 <- matrix(0, n.lags, n.vars)
#   m.B.1 <- l.B.j %>% reduce(cbind)
#   lagmat <- if( n.lags > 1 ) cbind( diag( n.vars*(n.lags-1) ), 
#                                     matrix(0, n.vars*(n.lags-1), n.vars ) ) else NULL
#   m.B <- rbind( m.B.1, lagmat )
#   # The reduced form VAR matrix
#   
#   # Create the IRFs
#   if(is.null(var.names)) var.names <- paste0('x.', 1:n.vars)
#   m.irf.select <- rbind( diag(n.vars), matrix(0, n.vars*(n.lags-1), n.vars ) )
#   df.irf <- suppressMessages(
#     lapply( (-1):n.pds, 
#          function(i) t(m.irf.select) %*% ( if(i<0) 0*m.B else m.B %^% i) %*% m.irf.select %>%
#            set_colnames(var.names) %>% 
#            set_rownames(paste0('shk.', 1:n.vars) ) %>%
#            as.data.frame %>%
#            rownames_to_column('shk') %>%
#            mutate(pd=i+1) ) %>%
#     reduce(rbind) %>%
#     select(shk, pd, everything()) %>%
#     arrange(shk, pd)
#   )
#   
#   return(df.irf)
# }
