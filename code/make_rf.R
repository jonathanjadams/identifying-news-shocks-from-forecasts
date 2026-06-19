##########################################################################################
#
# make_rf.R
#
# Code to define the functions that make the reduced form VARs
# Philip Barrett, Washington DC
# First version: 01apr2023
#
##########################################################################################

library(tidyverse)
library(vars)
library(expm)   
library(magrittr)
library(RcppEigen)  # Use to speed up the linear regression in here.
select <- dplyr::select
# Otherwise picks up the MASS usage

make.rf.var <- function( df, var.names, fcast.names, pd.name='pd', lags=1, 
                         m.Y.in=NA, restrict=TRUE, big.B.only=FALSE, do.aic=TRUE ){
# Makes a VAR and returns the reduced form coefficients. No bias correction.
  
  if( is.na(m.Y.in)){
    m.Y.in <- df %>%
      select(-!!pd.name) %>%
      # select(!!fcast.names, !!var.names) %>%
      as.matrix
    # The matrix of data for the VAR
  }
  n.var <- var.names %>% length()
  y.means <- m.Y.in %>% colMeans()
  # if( (y.means %>% abs %>% max) > 1e-10 ) message( 'Maximum mean bigger than consistent with deseasonalizing')
    # Should be v close to zero if de-seasonalizing correctly
  p <- lags
    # Because using the name "lags" causes errors later.  Weird.
  
  ## Augment the forecasts with data in the lag structure
  # browser()
  m.f.aug <- (VAR(m.Y.in, p=p ) %>% fitted())[,var.names]
    # Make the augmented forecasts from the fitted values for the VAR
  df.aug <- df
  df.aug[,fcast.names] <- rbind( matrix(NA, p-1, n.var), m.f.aug, NA )
    # The final "NA" shifts the timing to make sure we put forecasts in before
    # the priod that there is a realization.
  df.aug <- df.aug %>% na.omit()
    # The augmented forecast, dropping the requisite lags
  m.Y <- df.aug %>%
    select(!!fcast.names, !!var.names) %>%
    as.matrix
    # The new data Y
  m.f <- df.aug %>%
    select(!!fcast.names) %>%
    as.matrix
  m.x <- df.aug %>%
    select(!!var.names) %>%
    as.matrix
  
  ### Estimate the coefficients with the theoretical restrictions
  if(restrict | p==1 ){
    
    # browser()
    if((!restrict) & p==1) message('Cannot run unrestricted estimation for the one-lag case. Doing restricted instead.')
    
    # Compute the estimated B_i from the top half of the VAR(m-1) formulation
    st.fla.partial <- paste0( 'lag(', var.names, ',' )
    st.fla.rhs <- (if(p>1) sapply( 1:(p-1), function(i) paste0( st.fla.partial, i, ')' ) ) %>% 
                     paste( collapse=' + ' ) %>%
                     unlist()  %>% paste0( ' + ', . ) else c() ) %>%
      paste0( paste0( 'lag(', fcast.names, ')', collapse = ' + '), . ) %>%
      paste0( ' ~ 0 + ', . )
    # The RHS of the regression formula.  Includes variable lags for all
    # variables but only one lag of the forecast
    l.fcast.fla <- lapply( fcast.names, function(this.fcast) 
      paste0(this.fcast, st.fla.rhs ) %>% as.formula() )
    # browser()
    # fcast.lm <- lapply( l.fcast.fla, fastLm, data=df.aug )
    fcast.lm <- lapply( l.fcast.fla, lm, data=df.aug )
    big.B.est <- sapply( fcast.lm, coef ) %>% t() %>%
      set_colnames( sapply( 1:p, function(i) paste0( var.names, '.l', i ) ) %>% c() )
    l.B.j.est <- lapply( 1:p, function(i) big.B.est[,(i-1)*n.var+1:n.var] )
    B.1 <- l.B.j.est[[1]]
    # The estimated coefficients
    
    ## Make the VAR-1 form of the data generating process for x
    lagmat.x <- if( p > 1 ) cbind( diag( n.var*(p-1) ), 
                                   matrix(0, n.var*(p-1), n.var ) ) else NULL
    big.B.x <- rbind( big.B.est, lagmat.x ) %>% 
      set_rownames( c( var.names,  if( p > 1 ) sapply( 1:(p-1), function(i) paste0( var.names, '.l', i )) else NULL ) %>% c ) 
    
    ## Make the VAR-1 form of the joint process for x and f
    lagmat.f <- if( p > 2 ) cbind( diag( n.var*2*(p-2) ), 
                                   matrix(0, n.var*2*(p-2), 2*n.var ) ) %>% 
      set_rownames(colnames(big.B.x)[1:(n.var*2*(p-2))]) else NULL
    big.B.f <- matrix(0, n.var*2*(max(1,p-1)), n.var*2*(max(1,p-1)) )
    big.B.f[ 1:n.var, 1:n.var ] <- B.1
    big.B.f[ n.var+(1:n.var), (1:n.var) ] <- diag(n.var)
    if(p>1) for(i in 1:(p-1)) big.B.f[ 1:n.var, 2*(i-1)*n.var + n.var+1:n.var ] <- l.B.j.est[[i+1]]
    big.B.f[ -(1:(2*n.var)),  ] <- lagmat.f
    # The reduced form VAR matrix
    
    ## Estimate the error covariance
    f.resid <- sapply( fcast.lm, resid )
    x.resid <- (m.x[-1,] - m.f[-nrow(m.f),]) %>% 
      tail( f.resid %>% nrow )
    all.resid <- cbind(f.resid,x.resid)
    m.Sigma <- cov(all.resid)
    # The reduced form error covariance
    
  }
  
  ### Estimated via an unrestricted VAR
  if(!restrict & p != 1){
    
    # Compute the estimated B_i in eq (2)
    est.var <- VAR(m.Y, p = p-1, type='none')
    m.B.var <- sapply( est.var %>% coef, function(x) x[,'Estimate'] ) %>% t
    B.1 <- m.B.var[1:n.var, 1:n.var]
    l.B.j.est <- c( list(B.1), lapply(2:p, function(k) m.B.var[ 1:n.var, 2*(k-2)*n.var+n.var+1:n.var] ) ) 
    big.B.est <- l.B.j.est %>% reduce(cbind) %>%
      set_colnames( sapply( 1:p, function(i) paste0( var.names, '.l', i ) ) %>% c() )
    
    ## Make the VAR-1 form of the  <-  generating process for x
    lagmat.x <- if( p > 1 ) cbind( diag( n.var*(p-1) ), 
                                   matrix(0, n.var*(p-1), n.var ) ) else NULL
    big.B.x <- rbind( big.B.est, lagmat.x ) %>% 
      set_rownames( c( var.names,  if( p > 1 ) sapply( 1:(p-1), function(i) paste0( var.names, '.l', i )) else NULL ) %>% c ) 
    
    ## Make the VAR-1 form of the joint process for x and f
    lagmat.f <- if( p > 2 ) cbind( diag( n.var*2*(p-2) ), 
                                   matrix(0, n.var*2*(p-2), 2*n.var ) ) %>% 
      set_rownames(colnames(m.B.var)[1:(n.var*2*(p-2))]) else NULL
    big.B.f <- rbind( m.B.var, lagmat.f )
    # The reduced form VAR matrix
    
    ## Estimate the error covariance
    all.resid <- est.var %>% resid()
    f.resid <- all.resid[,1:n.var]
    x.resid <- all.resid[,-(1:n.var)]
    m.Sigma <- cov(all.resid)
    # The reduced form error covariance
  }
  if( big.B.only ) return(big.B.f)
    # Returns only the estimated B from the joint process
  
  # Compute the likelihood and AIC
  m.Omega <- cov(x.resid)
  big.T <- x.resid %>% nrow()
  if(do.aic){
    
  ll.x.form <- - big.T * n.var / 2 * log( 2 * pi ) - big.T / 2 * det(m.Omega, log=TRUE) -
    1/2 * (x.resid %*% solve( m.Omega, t(x.resid))) %>% diag() %>% sum()
  ll.f.form <- - big.T * n.var * log( 2 * pi ) - big.T / 2 * det(m.Sigma, log=TRUE) -
    1/2 * (all.resid %*% solve( m.Sigma, t(all.resid))) %>% diag() %>% sum()
  ll <- c( x.form=ll.x.form, f.form=ll.f.form )
    # The log likelihood (See 11.3.29 in Hamilton for details)
  k.x <- n.var * n.var * p + n.var * (1+n.var) / 2
  k.f <- (if(restrict) n.var * n.var * p else (2*n.var)^2 *p ) + 2*n.var * (1+2*n.var) / 2
  aic <- c( x.form= 2* k.x - 2*ll.x.form,
            f.form= 2* k.f - 2*ll.f.form )
  
  # Initialize the output
  out <- list( B=big.B.est, l.B=l.B.j.est, Sigma=m.Sigma, n.var=n.var, lags=p, 
               mu=y.means, data=m.Y, df=df.aug, big.B.x=big.B.x, big.B.f=big.B.f, big.resid=all.resid,
               loglik=ll, aic=aic,
               args=list( var.names=var.names, fcast.names=fcast.names, restrict=restrict,
                          pd.name=pd.name, lags=lags) )  
  }else{
    out <- list( B=big.B.est, l.B=l.B.j.est, Sigma=m.Sigma, n.var=n.var, lags=p, 
                 mu=y.means, data=m.Y, df=df.aug, big.B.x=big.B.x, big.B.f=big.B.f, big.resid=all.resid,
                 args=list( var.names=var.names, fcast.names=fcast.names, restrict=restrict,
                            pd.name=pd.name, lags=lags) )  
  }
  return( out )
}