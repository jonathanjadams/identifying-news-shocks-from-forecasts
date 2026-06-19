##########################################################################################
#
# policy_counterfactuals.R
#
# Code to compute the policy counterfactuals from the "Identification from 
# Rational Expectations" project
#
# Philip Barrett, Washington DC
# First version: 28jan2023
#
##########################################################################################

#### Libraries
library(tidyverse)
library(broom)

policy_counterfactual_obj <- 
  function( df.irfs, var.names, policy.shk.idces, v.F ){
# Policy counterfactual for the explicit objective function

    other.shk.idx <- setdiff( 1:length(var.names), policy.shk.idces )
    other.shks <- paste0( rep(c('u.','v.'), each=length(other.shk.idx) ), other.shk.idx )
      # The non-policy shocks
    df.F.mult <- df.irfs %>% select(pd, type, dimension, data.names) %>%
      mutate( data.F=as.matrix( {.} %>% select(data.names) ) %*% v.F %>% c )
      # The F-multiplied data
    df.lhs <- df.F.mult %>% select( pd, type, dimension, lhs=data.F ) %>% 
      filter(dimension %in% other.shk.idx ) %>%
      mutate( shk = paste0( type, '.', dimension) ) %>%
      select( pd, shk, lhs )
      # The LHS data
    df.rhs <- df.F.mult %>% 
      mutate( rhs = - data.F ) %>% 
        # Flip the signon the RHS variables
      select( pd, type, dimension, rhs ) %>% 
      filter(dimension %in% policy.shk.idces ) %>%
      mutate( shk.name = paste0( 'pol.shk.', type, '.', dimension) ) %>%
      select( pd, rhs, shk.name ) %>%
      spread( shk.name, rhs )
      # The RHS data
    df.reg <- left_join( df.lhs, df.rhs )
      # The regression dataset
    reg.fla <- paste0( 'lhs ~ 0 + ', paste( names(df.rhs)[-1], collapse=' + ' ) ) %>%
      as.formula()
      # The regression formula
    df.coef.lm <- df.reg %>% group_by(shk) %>% 
      do(tidy(lm( reg.fla, data=. )))
      # Do the main regression
    df.alpha <- df.coef.lm %>% select( shk, term, estimate ) %>%
      spread(term, estimate ) %>% ungroup() %>%
      set_names( gsub('pol.shk', 'alpha', names({.})))
    alpha.u <- df.alpha %>% 
      column_to_rownames('shk') %>% 
      select( contains( 'alpha.u' ) ) %>% 
      as.matrix
    alpha.v <- df.alpha %>% 
      column_to_rownames('shk') %>% 
      select( contains( 'alpha.v' ) ) %>% 
      as.matrix
      # Shape and extract the regression coefficients
    df.fit <- df.reg %>%
      full_join(df.alpha) %>%
      full_join( df.reg %>% group_by(shk) %>% do(augment(lm( reg.fla, data=. ))) ) %>%
      select( -(.hat:.std.resid) ) %>%
      rename( resid=.resid, fitted=.fitted )
      # The fitted regressions
    
    out <- list( alpha.u=alpha.u, alpha.v=alpha.v, df.alpha=df.alpha,
                 df.fit=df.fit, shk.idx=policy.shk.idces )
    return(out)
    
}



policy_counterfactual_rule <- 
  function( df.irfs, var.names, policy.var, policy.shk.idx, alpha.u=NULL, 
            alpha.v=NULL, phi.j=NULL ){
# Computes the policy counterfactual from the dataset of impulse responses. NB:
# Elements of phi.j *should* be labelled with the corresponding variables,
# otherwise  will be automatically assigned in order

# This is set up to minimize the error relative to a given *policy rule*, given
# by:
#      policy var = phi.j ' X
# The mapping to the the form:
#      min ||F*X||
# is given by:
#      F = phi.j - 1_{policy var}
# In practice this is implemented by defining the systematic and policy-driven 
# parts for the IRFs
#      Systematic part = phi.j ' X
#      Policy part = policy var
# Then we regress the difference between the systematic and policy parts for 
# the non-policy shocks on the same for the policy shock.  This gives:
#      F * [ IRFs of other shocks ] = - F * [ Policy shock IRFs ] * alpha + epsilson
# Just as in the paper. Notes:
#  1. To minimize the variance of a given variable, we only need set phi.j=0
#  2. To minimize the variance of a combination of variables, just apply a minus
#     sign to the relevant entry in phi.j, i.e. to minimize pi + k*y, set the 
#     policy variable to pi and phi.j to -k in the y entry.
  
  ## Switch for alpha.u and alpha.v provided
  if( !is.null(phi.j) ){
    
    #### TO DO HERE: CHECK IF THE CODE WILL WORK WITH MULTIPLE POLICY VARIABLES
    # Need to think about the policy.var thing
    # browser()
    
    # We have to compute the optimal alpha.u and alpha.v
    other.vars <- var.names[var.names!=policy.var]
      # The non-policy variables
    if( phi.j %>% names %>% is.null ) names(phi.j) <- other.vars
      # Set the names of phi.j if required
    other.shk.idx <- setdiff( 1:length(var.names), policy.shk.idx )
    other.shks <- paste0( rep(c('u.','v.'), each=length(other.shk.idx) ), other.shk.idx )
      # The non-policy shocks
    
    ## Extract the relevant objects
    df.sys <- df.irfs %>% 
      select(type, dimension, pd, all_of(other.vars)) %>%
      gather( variable, value, -pd, -type, -dimension ) %>%
      mutate( phi.j=phi.j[variable], shk=paste0(type, '.', dimension) ) %>%
      group_by(type, dimension, shk, pd) %>%
      summarise( sys = sum(phi.j * value)) %>%
      ungroup()
        # The systematic part of the policy rule
    df.reg <- df.sys %>%
      full_join( df.irfs %>% select( type, dimension, pd, 'policy'=policy.var ) ) %>%
      filter( !(dimension==policy.shk.idx) ) %>%
      full_join( df.sys %>% filter( dimension==policy.shk.idx, type=='u') %>%
                   select( pd, sys.u=sys) ) %>%
      full_join( df.irfs %>% filter( dimension==policy.shk.idx, type=='u') %>% 
                   select( pd, 'policy.u'=policy.var ) ) %>%
      full_join( df.sys %>% filter( dimension==policy.shk.idx, type=='v') %>%
                   select( pd, sys.v=sys) ) %>%
      full_join( df.irfs %>% filter( dimension==policy.shk.idx, type=='v') %>% 
                   select( pd, 'policy.v'=policy.var ) ) %>%
      mutate( Y=sys-policy, X.u=policy.u-sys.u, X.v=policy.v-sys.v )
        # The Dataset for the regression
    df.coef.lm <- df.reg %>% group_by(shk) %>% 
      do(broom::tidy(lm( Y~0+X.u+X.v, data=. )))
        # Do the main regression
    df.alpha <- df.coef.lm %>% select( shk, term, estimate ) %>%
      spread(term, estimate ) %>%
      rename(alpha.u=X.u, alpha.v=X.v )
    alpha.u <- df.alpha$alpha.u %>% set_names(df.alpha$shk)
    alpha.v <- df.alpha$alpha.v %>% set_names(df.alpha$shk)
      # Shape and extract the regression coefficients
    df.target <- df.reg %>%
      full_join(df.alpha) %>%
      mutate( fit = alpha.u * X.u + alpha.v * X.v,
              reg.err = Y - fit,
              target = sys + alpha.u * sys.u + alpha.v * sys.v,
              rule.fit = policy + alpha.u * policy.u + alpha.v * policy.v,
              rule.err = target - rule.fit,
              err.err=reg.err-rule.err )
    
    # m.alpha <- sapply( l.df.reg, function(df) lm( Y~0+X.u+X.v, df ) %>% coef() )
    # alpha.u <- m.alpha['X.u',]
    # alpha.v <- m.alpha['X.v',]
    
  }else{
    if( is.null(alpha.v) | is.null(alpha.u ) ){
      stop('Must supply alpha.v and alpha.u if no phi.j in policy counterfactual')
    }
  }
  
  # Now create some outputs
  
  out <- list( alpha.u=alpha.u, alpha.v=alpha.v, df.alpha=df.alpha,
               df.target=df.target, policy.var=policy.var, shk.idx=policy.shk.idx )
  return(out)
}

policy_irfs <- function(l.est, policy.rule, var.names=NULL, n.pds=20 ){
# Makes the impulse responses to the "other shocks when there is a policy rule
  
  policy.shk.idx <- policy.rule$shk.idx
    # The policy shock index
  policy.alpha <- policy.rule$df.alpha
    # The alphas on the policy rule
  
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
  non.policy.shks <- setdiff( 1:n.vars, policy.shk.idx )
  
  # Create the IRFs
  irf.fn <- function(i){
    
    this.u.u <- u.0 ; this.u.v <- u.0 ; 
    this.v.u <- v.0 ; this.v.v <- v.0
    
    this.pol.u <- policy.alpha %>% filter( shk==paste0('u.',i)) %>%
      select( -shk ) %>% unlist()
    this.u.u[1,i] <- sd.u[i]
    this.u.u[1,policy.shk.idx] <- this.pol.u[ grepl( 'alpha.u', names(this.pol.u) ) ] * sd.u[policy.shk.idx]
    this.v.u[2,policy.shk.idx] <- this.pol.u[ grepl( 'alpha.v', names(this.pol.u) ) ] * sd.v[policy.shk.idx]
      # The endogenous response to a surprise
    
    this.pol.v <- policy.alpha %>% filter( shk==paste0('v.',i)) %>% 
      select( -shk ) %>% unlist()
    this.v.v[2,i] <- sd.v[i]
    this.u.v[1,policy.shk.idx] <- this.pol.v[ grepl( 'alpha.u', names(this.pol.u) ) ] * sd.u[policy.shk.idx]
    this.v.v[2,policy.shk.idx] <- this.pol.v[ grepl( 'alpha.v', names(this.pol.u) ) ] * sd.v[policy.shk.idx]
      # The endogenous response to a news shock
    
    df.this.irf <- list( u=sim_var_fcast( x.0, l.B.j, A, C, NULL, NULL, n.pds, this.u.u, this.v.u ) %>%
                           mutate(type='u'),
                         v=sim_var_fcast( x.0, l.B.j, A, C, NULL, NULL, n.pds, this.u.v, this.v.v ) %>%
                           mutate(type='v') ) %>%
      reduce(rbind) %>%
      mutate(dimension=i)
    return(df.this.irf)
  }
  fcast.names <- if(is.null(var.names)) paste0('f.', 1:n.vars) else paste0( var.names, '_f' )
  if(is.null(var.names)) var.names <- paste0('x.', 1:n.vars)
  df.irfs <- suppressMessages( lapply( non.policy.shks, irf.fn ) %>% 
                                 reduce(rbind) %>%
                                 { if(n.vars>1) 
                                 {.} %>% rename_at(vars(c(paste0('x.', 1:n.vars), paste0('f.', 1:n.vars))),
                                                   ~ c(var.names, fcast.names) ) else {.} }
  )
  # Rstudio thinks that this is wrong but it is not
  
  return(df.irfs)
}


struct_convert <- function( this.struct, this.cf, policy.shk.idx ){
# Creates a new structural form of the DGP when the counterfactual policy is applied.
  
  out <- this.struct
    # Set up the output
  sd.u <- this.struct$D.u.2 %>% diag %>% sqrt
  sd.v <- this.struct$D.v.2 %>% diag %>% sqrt
    # Standard deviations.  Need these for scaling.
  
  this.pol.u.u <- this.cf$df.alpha %>% filter( grepl( 'u', shk) ) %>%
    select(contains('alpha.u')) %>% as.matrix
  this.pol.u.v <- this.cf$df.alpha %>% filter( grepl( 'u', shk) ) %>%
    select(contains('alpha.v')) %>% as.matrix
    # Responses to a surprise
  
  new.A <- out$A[,-policy.shk.idx] + 
    out$A[,policy.shk.idx] %*% (t(this.pol.u.u) * sd.u[policy.shk.idx]) %*% 
      diag( 1 / sd.u[-policy.shk.idx ] ) +
    out$C[,policy.shk.idx] %*%(t(this.pol.u.v) * sd.v[policy.shk.idx]) %*% 
      diag( 1 / sd.u[-policy.shk.idx ] )
    # The new A matrix -- add together the policy surprise and news.
  
  this.pol.v.u <- this.cf$df.alpha %>% filter( grepl( 'v', shk) ) %>%
    select(contains('alpha.u')) %>% as.matrix
  this.pol.v.v <- this.cf$df.alpha %>% filter( grepl( 'v', shk) ) %>%
    select(contains('alpha.v')) %>% as.matrix
    # Responses to news shocks
  
  new.C <- out$C[,-policy.shk.idx] + 
    out$A[,policy.shk.idx] %*% (t(this.pol.v.u) * sd.u[policy.shk.idx]) %*% 
      diag( 1 / sd.v[-policy.shk.idx ] ) +
    out$C[,policy.shk.idx] %*%(t(this.pol.v.v) * sd.v[policy.shk.idx]) %*% 
      diag( 1 / sd.v[-policy.shk.idx ] )
    # The new C matrix -- add together the policy surprise and news.
  out$A <- new.A ; out$C <- new.C
      # Replace the A and C matrices in the standard form
  out$D.u.2 <- this.struct$D.u.2[-policy.shk.idx,-policy.shk.idx]
  out$D.v.2 <- this.struct$D.v.2[-policy.shk.idx,-policy.shk.idx]
      # The variance vectors
  return(out)
}

cf_policy_pc <- function( df.irfs.pt, l.struct, df.labels, var.names, policy.shk.idx, 
                          p.vals, policy.var=NULL, alpha.u=NULL, alpha.v=NULL, 
                          phi.j=NULL, v.F=NULL ){
# Wrapper for the policy counterfactual
  
  n.boot <- l.struct$boot %>% length()
  n.pds <- df.irfs.pt$pd %>% max
  
  if( is.null(phi.j)){
    l.pol.cf <- policy_counterfactual_obj( df.irfs.pt, var.names, policy.shk.idx, v.F )
  }else{
    l.pol.cf <- policy_counterfactual( df.irfs.pt, var.names, policy.var, policy.shk.idx,
                                  alpha.u=alpha.u, alpha.v=alpha.v, phi.j=phi.j )
  }
  
  df.pol.pt <- policy_irfs( l.struct$pt, l.pol.cf, 
                            var.names, n.pds = n.pds.irf ) %>%
    left_join( df.labels %>% select( dimension, shk.name, irf.sign ) )
      # The point IRFs
  df.pol.boot <- lapply( 1:n.boot, function(this.boot) 
    policy_irfs( l.struct$boot[[this.boot]], l.pol.cf, 
                 data.names, n.pds = n.pds.irf ) %>% mutate(i.boot=this.boot) ) %>%
    rbindlist %>%
    left_join( df.labels %>% select( dimension, shk.name, irf.sign ) )
      # The bootstrapped IRFs
  df.pol.pc <-
    lapply( c( p.vals, rev(1-p.vals ) ), function(x)
      df.pol.boot %>% group_by( type, dimension, shk.name, irf.sign, pd ) %>%
        summarise_if( is.numeric, quantile, probs=x, na.rm=T ) %>%
        mutate( pc=paste0('pc.', x ) ) %>%
        select(-i.boot)
    ) %>%
    reduce(rbind) %>%
    select( type, dimension, pd, pc, everything() )
  df.joint <- full_join( df.pol.pc, df.pol.pt %>% mutate(pc='pt'))
  
  l.struct.cf.pt <- struct_convert( l.struct$pt, l.pol.cf, policy.shk.idx )
    # The VAR form
  m.var.pt <- (l.struct$pt %>% var.decomp(n.pds=n.pds))$var.decomp %>%
    sapply( ., function(x) x$epsilon %>% rowSums )
  m.var.cf.pt <- (l.struct.cf.pt %>% var.decomp(n.pds=n.pds))$var.decomp %>%
    sapply( ., function(x) x$epsilon %>% rowSums )
  m.var.ratio <- m.var.cf.pt / m.var.pt
  l.var.decomp <- list( orig=m.var.pt, cf=m.var.cf.pt, ratio=m.var.ratio )
    # The variance decompositions
  
  l.irf.norm <- lapply( list( orig=df.irfs.pt, cf=df.pol.pt ),
                         function(x) x %>% 
                           select( dimension, shk.name, pd, var.names ) %>% 
                           gather( variable, irf, var.names ) %>% 
                           group_by(dimension, shk.name, variable ) %>% 
                           summarise( irf_norm=sum(irf^2) ) %>% 
                           spread( variable, irf_norm ) %>% 
                           ungroup() )
  l.irf.norm$ratio <- l.irf.norm$cf %>% select(var.names) %>% colSums() / 
    l.irf.norm$orig %>% select(var.names) %>% colSums()
      # The simpler IRF norm
  
  return( list( cf=l.pol.cf, df=df.joint, 
                var.decomp=l.var.decomp, irf.norm=l.irf.norm, F=v.F ) ) 
}


policy_irfs_sd <- function(df.irfs.pc, policy.rule, timing.adj=FALSE ){
# Computes the policy IRFs and the associated standard deviation
  
  # browser()
  df.irfs.sd <- df.irfs.pc %>% 
    filter( pc %in% c('sd', 'pt.est' ), type!='epsilon' ) %>%
    select( -i.boot, -shk.wt ) %>% 
    mutate( shk=paste0(type, '.', dimension )) 
    # Select the point estimate and standard deviation
  if(timing.adj){
    df.irfs.sd <- df.irfs.sd %>%
      mutate( pd=ifelse(type=='v', pd+1, pd))
  }
    
  
  pol.shk.idx <- policy.rule$shk.idx
  df.pol.shk <- df.irfs.sd %>% 
    filter( dimension==pol.shk.idx ) %>% 
    gather( variable, value, -pc, -shk, -shk.name, 
            -type, -dimension, -pd, -irf.sign ) %>% 
    ungroup() %>%
    select( pd, pc, type, variable, value ) %>%
    spread(type, value) %>%
    rename(u.pol=u, v.pol=v)
    # The policy shock in long form
  df.other.shks <- df.irfs.sd %>% 
    filter( dimension!=pol.shk.idx ) %>% 
    gather( variable, value, -pc, -shk, -shk.name, 
            -type, -dimension, -pd, -irf.sign )
    # The other shocks
  df.pol.cf <- df.other.shks %>%
    left_join( df.pol.shk ) %>%
    left_join( policy.rule$df.alpha ) %>%
    mutate( alpha.u=ifelse(pc=='sd', abs(alpha.u), alpha.u),
            alpha.v=ifelse(pc=='sd', abs(alpha.v), alpha.v) ) %>%
    mutate(pol.cf=value+alpha.u*u.pol+alpha.v*v.pol)
  df.out <- df.pol.cf %>%
    select( -value, -u.pol, -v.pol, -alpha.u, -alpha.v ) %>%
    spread( variable, pol.cf ) %>%
    select( pd, everything(), type, dimension, shk.name, irf.sign) %>%
    arrange( pc, type, dimension, pd)
  
  return(df.out)
  
}

policy_history <- function(){
# Makes the alternate history when there is a policy rule
  
}
