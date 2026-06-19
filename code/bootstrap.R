##########################################################################################
#
# bootstrap.R
#
# Code to run bootstrap for the "Identification from Rational Expectations" project
#
# Philip Barrett, Washington DC
# First version: 19jul2022
#
##########################################################################################

rf_debias <- function( l.rf, method='resample', n.boot=1000, burn=100, pt.est.only=FALSE ){
# De-biases the reduced form estimates
  
  ## 1. Compute the bias in the point estimates ##
  message('Stage 1: Correcting bias in point estimates')
  rf.B.debiased <- rf_debias_1st_stage( l.rf, method, n.boot, burn )
    # The debiased mean
  message( '** Maximum absolute root of the reduced form VAR: ** ', 
           '\n  inital B             = ', round(rf.B.debiased$df.bias.persist$B.in[1], 2),
           '\n  bias, bootstrap avg. = ', round(rf.B.debiased$df.bias.persist$bias.avg[1], 4),
           '\n  rmse, bootstrap avg. = ', round(rf.B.debiased$df.bias.persist$rmse[1], 5),
           '\n  corrected B          = ', round(rf.B.debiased$df.bias.persist$B.adj[1], 2),
           '\n  correction error     = ', round(rf.B.debiased$df.bias.persist$adj.err[1], 5) )
  
  ## 2. Fill in the other elements of l.rf.adj ##
  n.var <- l.rf$n.var ; lags <- l.rf$lags
  nonzero.col.idces <- c(1:n.var, if(lags>1) 1:n.var + n.var*rep( 2*1:(lags-1)-1, each=n.var ) )
    # Set up
  make.rf.est <- function( this.B, this.Sigma, this.resid, this.data ){
  # Assemble everything into a single list
    B.out <- this.B[1:n.var, nonzero.col.idces]
    l.B.out <- lapply( 1:lags, function(i) B.out[ , n.var*(i-1)+1:n.var ])
    big.B.x.out <- l.rf$big.B.x ; big.B.x.out[1:n.var,] <- B.out
      # The various formulations of th B matrices
    x.resid <- this.resid[,n.var+1:n.var]
    m.Omega <- cov(x.resid)
      # Just the x residual covariance
    big.T <- this.resid %>% nrow()
    ll.x.form <- - big.T * n.var / 2 * log( 2 * pi ) - big.T / 2 * det(m.Omega, log=TRUE) -
      1/2 * (x.resid %*% solve( m.Omega, t(x.resid))) %>% diag() %>% sum()
    ll.f.form <- - big.T * n.var * log( 2 * pi ) - big.T / 2 * det(rf.B.debiased$Sigma, log=TRUE) -
      1/2 * (this.resid %*% solve( rf.B.debiased$Sigma, t(this.resid))) %>% diag() %>% sum()
    ll.out <- c( x.form=ll.x.form, f.form=ll.f.form )
      # The log likelihood (See 11.3.29 in Hamilton for details)
    k.x <- 2 * (n.var * n.var * lags) + n.var * (1+n.var) / 2
    k.f <- (if(l.rf$args$restrict) 2 * (n.var * n.var * lags) else (2*n.var)^2 * lags ) + 2*n.var * (1+2*n.var) / 2
      # Need to double the number of parameters from the lag structure to
      # account for the forecast estimation
    aic.out <- c( x.form= 2* k.x - 2*ll.x.form,
                  f.form= 2* k.f - 2*ll.f.form )
      # The AIC
    lambda <- eigen(this.B)$values %>% abs()
    
    l.out <-  list( B=B.out, l.B=l.B.out, Sigma=this.Sigma, n.var=n.var,
                    lags=lags, data=this.data, big.B.x=big.B.x.out, 
                    big.B.f=this.B, big.resid=this.resid, args=l.rf$args,
                    lambda=lambda, ll=ll.out, aic=aic.out )
  }
  l.rf.adj <- make.rf.est( rf.B.debiased$B, rf.B.debiased$Sigma, rf.B.debiased$resid, l.rf$data )
    # Create the reduced form estimate list for the debiased estimate using the data
  if( pt.est.only ) return(list( rf=l.rf.adj, stage.1=rf.B.debiased$df.bias.persist ) )
  
  ## 3. Compute the bootstrap for the de-biased estimates for B ##
  message('Stage 2: Resampling from unbiased distribution')
  l.boot.B.Sigma <- rf_debias_2nd_stage( l.rf.adj, method, n.boot, burn )
    # Re-bootstrap to get a distribution for B and Sigma.
  stage.2.err <- rf.B.debiased$B-l.boot.B.Sigma$B.avg
  stage.2.err.check <- (stage.2.err[stage.2.err!=0]) %>% 
    (function(x) c( mean=mean(x), mean.abs=mean(abs(x)), max=max(abs(x)) ) )
  l.rf.boot <- lapply( 1:n.boot, 
      function(i) make.rf.est( l.boot.B.Sigma$l.B[[i]], l.boot.B.Sigma$l.Sigma[[i]], 
                               l.boot.B.Sigma$resid[[i]], l.boot.B.Sigma$data[[i]] ) )
  
  # Return output
  l.out <- list( rf=l.rf.adj, boot=l.rf.boot, 
                 stage.1=rf.B.debiased$df.bias.persist,
                 stage.2=l.boot.B.Sigma$df.lambda,
                 stage.2.err.check=stage.2.err.check)
  return(l.out)  
  
}

rf_debias_1st_stage <- function( l.rf, method='resample', n.boot=1000, burn=500 ){
# Computes the bias in the reduced form VAR by bootstrap method.  First stage:
# take the incoming RF model as biased and correct it using the bootstrap
  
  ## Set up
  rf.B <- l.rf$big.B.f ; lags <- l.rf$lags ; n.var <- l.rf$n.var
  
  # Check for stationarity
  rf.roots <- 1 / abs(eigen(rf.B)$values)
  if(any(rf.roots<1)){
    return( list( big.B.f=rf.B, Sigma=l.rf$Sigma ) )
  }
  
  if( method=='analytic' ){
  # Uses Pope's (1989) formula for analytical bias-correction
    if(lags<=2){
      m.data.lags <- l.rf$data
    }else{
      m.data.lags <- l.rf$data
      for(i in 1:(lags-2) ) m.data.lags <- cbind( m.data.lags, 
                                                  rbind( matrix(NA, nrow=i, ncol=2*n.var ), l.rf$data[-(nrow(l.rf$data)+1-(1:i)),]  ) )
    }
    cov.data <- m.data.lags[!apply(is.na(m.data.lags),1,any),] %>% cov
        # The covariance of the data
    TT <- m.data.lags %>% nrow()
    m.Sigma <- if( lags > 2 ) cbind( l.rf$big.resid, 
                   matrix( 0, nrow=nrow(l.rf$big.resid), ncol=2*n.var*(lags-2) )) %>% cov() else l.rf$Sigma
    k <- rf.B %>% nrow() ; B.t <- t(rf.B)
    eigenvalues <- ( rf.B %>% eigen )$values
      # Other preliminaries
    ev.sum <- lapply( eigenvalues, function(x) x * solve( diag(k) - x * B.t ) ) %>%
      Reduce( '+', . ) %>% Re
      # The the eigenvalue-weighted sum.  The imaginary part is zero due to
      # complex conjugates in the eigenvalues.
    m.bias <- - 1 / TT * 
      m.Sigma %*% ( solve( diag(k) - B.t ) + B.t %*% 
                      solve( diag(k) - B.t %*% B.t ) + ev.sum ) %*% solve(cov.data)
      # This doesn't seem to work.  I think because of the restrictions on our
      # estimation.  No further usage of this.
  }else{
    # Create the bootstrapped B estimates
    l.B.boot <- rf_bootstrap_B( l.rf, method, n.boot, burn )
    
    ## Compute the bias and associated statistics
    B.avg.boot <- ( l.B.boot %>% Reduce('+',.) ) / n.boot
    m.bias <- B.avg.boot - rf.B
    sd.bias <- sapply( l.B.boot, function(x) c( x - rf.B ) ) %>%
      apply( ., 1, sd ) / sqrt(n.boot)
        # The standard deviation of the average bias
    v.persist.avg.B <- abs(eigen(B.avg.boot)$values) # - 1 / rf.roots
    v.avg.persist <- sapply( l.B.boot, function(x) abs(eigen(x)$values) ) %>%
      apply( ., 1, mean )
    se.avg.persist <- sapply( l.B.boot, function(x) abs(eigen(x)$values) ) %>%
      apply( ., 1, sd ) / sqrt(n.boot)
  }
  
  
  B.adj <- rf_killian_stationarize(rf.B, m.bias)
  df.bias.persist <- data.frame( B.in = 1/rf.roots, 
                                 B.boot.avg=v.persist.avg.B, bias.avg=v.persist.avg.B-1/rf.roots,
                                 B.adj=abs(eigen(B.adj)$values), rmse=se.avg.persist ) %>%
    mutate( adj=B.adj-B.in,
            adj.err= bias.avg + adj )
  
  ## Compute the corresponding covariance matrix
  if(lags<=2){
    m.data.lags <- l.rf$data
  }else{
    m.data.lags <- l.rf$data
    for(i in 1:(lags-2) ) m.data.lags <- cbind( m.data.lags, 
                              rbind( matrix(NA, nrow=i, ncol=2*n.var ), l.rf$data[-(nrow(l.rf$data)+1-(1:i)),]  ) )
  }
  m.pred <- t( B.adj %*% t( m.data.lags ) )
  m.resid <- l.rf$data[-1,] - m.pred[-nrow(l.rf$data),1:(2*n.var)]
  m.resid <- m.resid[!apply(is.na(m.resid), 1, any), ]
  m.Sigma <- cov( m.resid )
  
  ## Return the output
  return( list( B=B.adj, bias=m.bias, adj=B.adj-rf.B, sd.bias=sd.bias, rf.roots=rf.roots, 
                df.bias.persist=df.bias.persist, Sigma=m.Sigma, resid=m.resid ) )
  
}

rf_debias_2nd_stage <- function( l.rf, method='resample', n.boot=1000, burn=500 ){
# Computes the bias in the reduced form VAR by bootstrap method.  Second stage:
# take the incoming RF model as de-biased and then adjust the bootstrap.
  
  ## Set up
  rf.B <- l.rf$big.B.f ; rf.roots <- 1 / abs(eigen(rf.B)$values)
  n.var <- l.rf$n.var
  
  # Create the bootstrapped B estimates
  l.boot <- rf_bootstrap_B( l.rf, method, n.boot, burn, include.sim=TRUE )
    # Need to return the simulations to get the post-correction Sigmas
  
  ## Compute the bias in the estimation procedure
  message(' Stationarizing bootstrapped estimates...')
  B.avg.boot <- ( l.boot$B %>% Reduce('+',.) ) / n.boot
  m.bias <- B.avg.boot - rf.B
  lambda.boot <- sapply( l.boot$B, function(x) eigen(x)$values %>% abs ) %>%
    apply( ., 1, mean )
  l.B.adj <- lapply( l.boot$B, function(x) rf_killian_stationarize( x, m.bias ) )
    # Correct the bootstrap to match the reduced form estimates. NB: This may
    # still fail if a draw from the 
  m.B.adj.avg <- ( l.B.adj %>% Reduce('+',.) ) / n.boot
  m.err <- m.B.adj.avg - rf.B
  message('    ...complete')
  
  ## Calculate impact on the persistence
  lambda.B.adj.avg <-  abs(eigen(m.B.adj.avg)$values)
  lambda.adj.avg <- sapply( l.B.adj, function(x) abs(eigen(x)$values) ) %>%
    apply( ., 1, mean ) 
  se.lambda.adj.avg <- sapply( l.B.adj, function(x) abs(eigen(x)$values) ) %>%
    apply( ., 1, sd ) / sqrt(n.boot)
  df.lambda <- data.frame( lambda.dgp = 1/rf.roots,  lambda.boot=lambda.boot,
                           err.boot=lambda.boot-1/rf.roots, 
                           lambda.B.adj.avg = lambda.B.adj.avg, err.B.adj=lambda.B.adj.avg-1/rf.roots,
                           lambda.adj.avg=lambda.adj.avg, err.adj.avg=lambda.adj.avg-1/rf.roots, 
                           se.adj.avg=se.lambda.adj.avg )
                                # err.avg=v.err.persist.avg, 
                                # avg.err=v.avg.err.persist, se=se.err.persist )
  
  ## Compute the simulated Sigmas
  l.resid <- lapply( 1:n.boot, function(i.boot){
    this.sim <- l.boot$sim[[i.boot]]
    this.pred <- t( l.B.adj[[i.boot]] %*% t( this.sim ) )
    this.resid <- (this.sim[-1,] - this.pred[-nrow(this.sim),])[,1:(2*n.var)]
    return(this.resid)
  } )
  l.Sigma <- lapply( l.resid, cov )
  
  return( list( l.B=l.B.adj, bias=m.bias, err=m.err, B.avg=m.B.adj.avg, 
                l.Sigma=l.Sigma, resid=l.resid, data=l.boot$sim, df.lambda=df.lambda ) )
  
}

rf_bootstrap_B <- function( l.rf, method='resample', n.boot=1000, burn=500, 
                            include.sim=FALSE ){
# Makes the bootstrapped persistence matrix (and returns the sim if requested)
  
  set.seed(42)
  
  ## Set up
  rf.B <- l.rf$big.B.f ; rf.resid <- l.rf$big.resid ; rf.Sigma <- l.rf$Sigma
  dim.B <- rf.B %>% nrow ; n.pds <- rf.resid %>% nrow
  lags <- l.rf$lags ; n.var <- l.rf$n.var
  
  ## Main loop
  l.B <- list()
  if(include.sim) l.sim <- list()
  
  message('  Bootstrapping reduced form... ')
  for( i.boot in 1:n.boot ){
    
    if( i.boot %% 200 == 0 ) message( '   Iteration ', i.boot)
    
    # Construct the residuals
    this.resid <- if(method=='resample') rf.resid[ 
      sample( 1:n.pds, n.pds+burn, TRUE ), ] else rmvnorm( n.pds+burn, sigma=rf.Sigma )
    # The sample of residuals
    if( lags > 2 ) this.resid <- 
      cbind( this.resid, matrix( 0, nrow=n.pds+burn, ncol=2*n.var*(lags-2) ))
    # Add empty residuals for lags if needed
    
    # Create the simulation
    this.sim <- rf_sim( rf.B, this.resid, burn )
    this.sim.no.lags <- this.sim[,1:(2*n.var)]
    # So this drops the extra columns as redundant if lag>1
    
    # Estimate and store the reduced form coefficients
    this.df <- as.data.frame( cbind( 1, this.sim.no.lags) ) %>%
      set_names( c(l.rf$args$pd.name, l.rf$args$fcast.names, l.rf$args$var.names) )
    l.B[[i.boot]] <- suppressMessages(
      make.rf.var( df=this.df, var.names=l.rf$args$var.names,
                   fcast.names=l.rf$args$fcast.names, pd.name=l.rf$args$pd.name, 
                   lags=l.rf$args$lags, restrict=l.rf$args$restrict, big.B.only=TRUE ) )
    if(include.sim) l.sim[[i.boot]] <- this.sim

  }
  
  l.out <- if(include.sim) list( B=l.B, sim=l.sim ) else l.B
    # Construct the output
  return(l.out)
}

rf_sim <- function( big.B, resid, burn=100 ){
# Simulates the reduced form VAR.  NB: returns forecast first

  dim.b <- big.B %>% ncol
  n.pd.burn <- resid %>% nrow
  sim <- matrix( NA, dim.b, n.pd.burn )
  sim[,1] <- c(resid[1,])
    # Flip sim for speed
  for( i in 2:n.pd.burn ) sim[,i] <- big.B %*% sim[,i-1] + c(resid[i,])
    # The main iteration
  out <- t( sim[,-(1:burn)])
  return(out)
}

rf_killian_stationarize <- function( B.orig, bias, delta.k=.05){
# Killian's stationarizing algorithm.
  
  orig.roots <- 1 / abs(eigen(B.orig)$values)
  if( any(orig.roots<1) ) return(B.orig)
  
  jj <- -1 ; stab <- FALSE ; jj.max <- 10000
  while( !stab & (jj < jj.max) ){
    jj <- jj + 1
    B.adj <- B.orig - (1-delta.k) ^ jj * bias
    lambda <- abs(eigen(B.adj)$values)
    adj.roots <- 1 / lambda
    stab <- all( adj.roots > 1 )
  }
  return( B.adj )
}























# 
# 
# 
# make_bootstrap_irf <- function( l.est, n.pds=20, n.boot=100, n.sim.length=NA, seed=42, unit='sd',
#                                 pc=c(.025, .05, .5, .95, .975 ), #struct.match='impact', 
#                                 delta.method=FALSE, return.all=FALSE, sim.coef=FALSE, 
#                                 n.match.pds=4, news.match=FALSE, sim.scale=1, fix.size=FALSE){
# # Computes bootstrapped distribution of impulse responses from an estimated VAR
#   
#   set.seed(seed)
#   
#   # Set up
#   if(is.na(n.sim.length)) n.sim.length <- l.est$data %>% nrow
#   n.x <- l.est$A %>% nrow
#   n.lags <- l.est$l.B %>% length()
#   x.0 <- matrix(0, n.lags, n.x)
#   
#   # Make the loss function for ordering the shocks
#   # if(struct.match=='impact'){
#   #   
#   #   if(news.match){
#   #     loss.fn <- function( A, C, m.D.v.2, m.D.u.2, l.B ){
#   #       # NB: Here the B matrices are redundant
#   #       if( any( (m.D.v.2 > 1-1e-10) != (l.est$D.v.2 > 1-1e-10) ) ){
#   #         return(Inf)
#   #       }else{
#   #         return( sum( (A %*% m.D.u.2 - l.est$A %*% l.est$D.u.2 )^2 ) +
#   #                      sum( (C %*% m.D.v.2 - l.est$C %*% l.est$D.v.2 )^2 ) )
#   #       }
#   #     }      
#   #   }else{
#   #     loss.fn <- function( A, C, m.D.v.2, m.D.u.2, l.B ){
#   #       # NB: HEre the D.u.2 and the B matrices are redundant
#   #       if( any( (m.D.v.2 > 1-1e-10) != (l.est$D.v.2 > 1-1e-10) ) ){
#   #         return(Inf)
#   #       }else{
#   #         # return( sum( (A-l.est$A)^2 )  )
#   #         return( sum( (A %*% m.D.u.2 - l.est$A %*% l.est$D.u.2 )^2 )  )
#   #       }
#   #     }
#   #   }
#   # }else if(struct.match=='irf'){
#     
#   pt.est.irf <- irf_var_fcast( l.est, n.pds = n.match.pds )
#   loss.target <- pt.est.irf %>% 
#     arrange( dimension, type, pd ) %>% 
#     select(dimension, type, pd, contains('x.')) %>% 
#     gather( variable, value, -dimension, -type, -pd ) %>%
#     spread( dimension, value ) %>%
#     select( -type, -pd, -variable ) %>%
#     as.matrix
#     # The IRF target
#   loss.fn <- function( cand.irf ){
#     err <- ( cand.irf - loss.target )
#     return( err^2 %>% sum() )
#   }
#   
#   # loss.fn <- function( A, C, m.D.v.2, m.D.u.2, l.B ){
#   #   
#   #   if( any( (m.D.v.2 > 1-1e-10) != (l.est$D.v.2 > 1-1e-10) ) ){
#   #     return(Inf)
#   #   }else{
#   #     this.est <- list( A=A, l.B=l.B, C=C, D.v.2=m.D.v.2, D.u.2=m.D.u.2 )
#   #     this.irf <- irf_var_fcast( this.est, n.pds = n.match.pds )
#   #       # The IRF for the given parameters
#   #     err <- ( this.irf %>% filter(type=='v') %>% 
#   #                select( contains('x.') ) %>% as.matrix ) - loss.target
#   #     return( err^2 %>% sum )
#   #     
#   #   }
#     # }
#     
#   #   
#   # }else{
#   #   loss.fn <- NULL
#   # }
#   
#   
#   # Compute the point estimate of the IRF
#   df.pt.irf <- 
#   irf_var_fcast( l.est, var.names=l.est$args$var.names, n.pds=n.pds, unit=unit ) %>%
#     mutate( type=as.factor(type) ) %>%
#     mutate_if(is.numeric, funs(ifelse(is.nan(.), 0, .))) %>%
#     mutate( pc='pt.est' ) 
#   
#   #### Switch based on whether to simulate the RF coefficients or not ###
#   if( sim.coef | delta.method ){
#     
#     # browser()
#     
#     ## Simulate the coefficients from the reduced form asymptotic variance.
#     
#     # Make the bootstrapped coefficient matrices
#     n.var <- l.est$n.var
#     v.B.sim <- mvtnorm::rmvnorm(n.boot, l.est$v.B$pt, l.est$v.B$vcov )
#     l.B.sim <- lapply( 1:n.boot, function(i) matrix(v.B.sim[i,], nrow=n.var, ncol=n.var*l.est$lags, byrow = TRUE ) )
#     l.l.B.sim <- lapply(l.B.sim, function(this.B) lapply( 1:l.est$lags,
#                           function(i) this.B[,(i-1)*n.var+1:n.var] ) )
#     l.Sigma.sim <- lapply( l.l.B.sim, function(this.l.B) make.resid( l.est$data, this.l.B) %>% cov() )
#     
#     if(delta.method){
#       # browser()
#       # Use the delta method to work out the covariance of the structural parameters
#       l.B.Sigma.sim <- cbind( v.B.sim, sapply(l.Sigma.sim, c ) %>% t )
#       vcov.b.Sigma <- cov(l.B.Sigma.sim)
#         # Compute the joint coefficients
#       struct.grad <- suppressMessages( structural_identification_grad( l.est$l.B, l.est$Sigma, ordering.loss.fn = loss.fn ) )
#       pt.est.struct <- vectorize_structural_identification( l.est$l.B, l.est$Sigma, ordering.loss.fn = loss.fn )
#       vcov.struct <- t(struct.grad) %*% vcov.b.Sigma %*% struct.grad
#       if( fix.size ){
#         vcov.struct[2*n.var^2 + 1:(2*n.var),] <- 
#           vcov.struct[,2*n.var^2 + 1:(2*n.var)] <- 0
#             # Set the shock *size* to a constant
#       }
#       # browser()
#       m.struct.sim <- mvtnorm::rmvnorm(n.boot, pt.est.struct, vcov.struct )
#       l.struct.id <- lapply( 1:n.boot, function(i) devectorize_struct( m.struct.sim[i,], l.est$n.var, l.est$lags ) )
#     }else{
#       # Produce the identified parameters
#       message('Computing strutural identification for simulated coefficients')
#       l.struct.id <- lapply( 1:n.boot, function(i){
#         if(i %% 20 == 0) message('n=', i)
#         out <- structural_identification( l.l.B.sim[[i]], l.Sigma.sim[[i]], ordering.loss.fn = loss.fn )
#         # out$l.B.j <- lapply( 1:l.est$lags, function(k) l.B.sim[[i]][,(k-1)*n.var + 1:n.var ] )
#         return(out)
#       } )
#       if( fix.size ){
#         # browser()
#         for( i in 1:n.boot ){
#           l.struct.id[[i]]$D.u.2 <- l.est$D.u.2
#           l.struct.id[[i]]$D.v.2 <- l.est$D.v.2
#             # Fix the size of the shocks
#         }
#       }
#     }
#     
#     ##### ALTERNATE IDEA: CREATE THE RESIDUALS AND MAKE A NUMERICAL COVARIANCE ######
#     ##### TO DO: FIX THIS #####
#     # v.Sigma.sqrt.sim <- mvtnorm::rmvnorm(n.boot, l.est$v.Sigma.sqrt$pt, l.est$v.Sigma.sqrt$vcov )
#     # l.Sigma.sqrt.sim <- lapply( 1:n.boot, function(i) matrix( v.Sigma.sqrt.sim[i,], 2*l.est$n.var, 2*l.est$n.var ) )
#     # l.Sigma.sim <- lapply( 1:n.boot, function(i) l.Sigma.sqrt.sim[[i]] %*% t(l.Sigma.sqrt.sim[[i]]) )
#     # Can't make the resampling of the shock variance :(
#     
#     # Make the full set of bootstrapped IRFs
#     message('\nCreating IRFs')
#     l.irfs <-
#       lapply( 1:n.boot, function(i){
#         if(i %% 20 == 0) message('n=', i)
#         irf_var_fcast( l.struct.id[[i]], var.names=l.est$args$var.names, n.pds=n.pds, unit=unit ) %>%
#           mutate( i.boot=i, type=as.factor(type) ) %>%
#           mutate_if(is.numeric, funs(ifelse(is.nan(.), 0, .)))
#       } )
#     # Use the simulation function -> Estimation function -> IRF function to
#     # produce the distribution of IRFs
#   }else{
#     
#     ## Simulate the new dataset from the point estimate parameters
#     
#     # Make the full set of bootstrapped IRFs
#     shk.size <- if(fix.size) list( u=diag(l.est$D.u.2)%>% sqrt, v=diag(l.est$D.u.2)%>% sqrt ) else NA
#     if( fix.size ) unit <- 'custom'
#     # browser()
#     
#     l.irfs <-
#       lapply( 1:n.boot, function(i){
#         message('n=', i)
#         sim_var_fcast( x.0=x.0, l.B.j = l.est$l.B, A=l.est$A, C=l.est$C, 
#                      D.u.2=l.est$D.u.2, D.v.2=l.est$D.v.2, n.pds=n.sim.length ) %>%
#           rename_with( ~l.est$args$var.names, contains('x.') ) %>%
#           rename_with( ~l.est$args$fcast.names, contains('f.') ) %>%
#           rename_with( ~l.est$args$pd.name, pd) %>%
#           select(c(l.est$args$pd.name, l.est$args$var.names, l.est$args$fcast.names)) %>%
#           make.struct.var( ., 
#                            var.names=l.est$args$var.names, fcast.names=l.est$args$fcast.names,
#                            pd.name=l.est$args$pd.name, lags=l.est$args$lags, #lag.max=l.est$args$lag.max,
#                            lag.select=l.est$args$lag.select, ordering.loss.fn=loss.fn, 
#                            n.match.pds = n.match.pds, sign.var = NA ) %>%
#             # NB: Need sign.var=NA to prevent unwanted flipping of IRFs
#         # lags = lags, restrict=TRUE,  n.pds.var.decomp = n.pds.irf )
#           irf_var_fcast( ., var.names=l.est$args$var.names, n.pds=n.pds, unit=unit, shk.size=shk.size ) %>%
#           mutate( i.boot=i, type=as.factor(type) ) %>%
#           mutate_if(is.numeric, funs(ifelse(is.nan(.), 0, .)))
#       } )
#   }
#   
#   df.irfs <- suppressMessages( l.irfs %>% reduce( rbind ) )
#     # Use the simulation function -> Estimation function -> IRF function to
#     # produce the distribution of IRFs
#   
#   if(return.all) return(df.irfs )
#   
#   # Summarize the IRFs by percentile
#   # browser()
#   df.irf.pc <- make.irf.pc( df.irfs, pc ) %>%
#       rbind( df.pt.irf )
#     
#   #   suppressMessages(
#   #   lapply( pc, function(x)
#   #     df.irfs %>% group_by( type, dimension, pd ) %>%
#   #       summarise_if( is.numeric, quantile, probs=x, na.rm=T ) %>%
#   #       mutate( pc=paste0('pc.', x ) )
#   #     ) %>%
#   #     reduce(rbind) %>%
#   #     rbind( df.irfs %>% group_by( type, dimension, pd ) %>%
#   #              summarise_if( is.numeric, sd, na.rm=T ) %>%
#   #              mutate( pc='sd' ) ) %>%
#   #     rbind( df.irfs %>% group_by( type, dimension, pd ) %>%
#   #              summarise_if( is.numeric, sd_trim, na.rm=T ) %>%
#   #              mutate( pc='sd.trim' ) )
#   # )
#   
#   return( df.irf.pc )
# }

# make_irf_pc <- function( df.irfs, pc=c(.025, .05, .5, .95, .975 ) ){
# # Makes the percentiles of the IRFs
#   
#   df.irf.pc <- suppressMessages(
#     lapply( pc, function(x)
#       df.irfs %>% group_by( type, dimension, pd ) %>%
#         summarise_if( is.numeric, quantile, probs=x, na.rm=T ) %>%
#         mutate( pc=paste0('pc.', x ) )
#     ) %>%
#       reduce(rbind) %>%
#       rbind( df.irfs %>% group_by( type, dimension, pd ) %>%
#                summarise_if( is.numeric, sd, na.rm=T ) %>%
#                mutate( pc='sd' ) ) %>%
#       rbind( df.irfs %>% group_by( type, dimension, pd ) %>%
#                summarise_if( is.numeric, sd_trim, na.rm=T ) %>%
#                mutate( pc='sd.trim' ) )
#   )
#   
#   return(df.irf.pc)
# }

# make_bootstrap_irf_coef <- function( l.est, n.pds=20, n.boot=100, n.sim.length=NA, seed=42, unit='sd',
#                                 pc=c(.025, .05, .5, .95, .975 ), struct.match=TRUE ){
# # Computes bootstrapped distribution of impulse responses from an estimated VAR
# # by drawing alternate coefficient matrices directly.
#   
#   set.seed(seed)
#   
#   # Make the loss function for ordering the shocks
#   browser()
#   if(struct.match=='impact'){
#     loss.fn <- function( A, C, m.D.v.2 ){
#       
#       if( any( (m.D.v.2 > 1-1e-10) != (l.est$D.v.2 > 1-1e-10) ) ){
#         return(Inf)
#       }else{
#         return( sum( (A-l.est$A)^2 )  ) # + sum( (C-l.est$C)^2 ) )
#       }
#     }
#   }else if(struct.match=='irf'){
#     # TO COMPLETE
#   }else{
#     loss.fn <- NULL
#   }
#   
#   
#   # Make the bootstrapped coefficient matrices
#   v.B.sim <- mvtnorm::rmvnorm(n.boot, l.est$v.B$pt, l.est$v.B$vcov )
#   l.B.sim <- lapply( 1:n.boot, function(i) matrix(v.B.sim[i,], nrow=l.est$n.var, ncol=l.est$n.var*l.est$lags, byrow = TRUE ) )
#   # v.Sigma.sqrt.sim <- mvtnorm::rmvnorm(n.boot, l.est$v.Sigma.sqrt$pt, l.est$v.Sigma.sqrt$vcov )
#   # l.Sigma.sqrt.sim <- lapply( 1:n.boot, function(i) matrix( v.Sigma.sqrt.sim[i,], 2*l.est$n.var, 2*l.est$n.var ) )
#   # l.Sigma.sim <- lapply( 1:n.boot, function(i) l.Sigma.sqrt.sim[[i]] %*% t(l.Sigma.sqrt.sim[[i]]) )
#       # Can't make the resampling of the shock variance :(
#   
#   # Produce the identified parameters
#   n.var <- l.est$n.var
#   l.B.1 <- lapply( l.B.sim, function(x) x[,1:n.var] )
#   l.struct.id <- lapply( 1:n.boot, function(i){
#     out <- structural_identification( l.B.1[[i]], l.est$Sigma, loss.fn )
#     out$l.B.j <- lapply( 1:l.est$lags, function(k) l.B.sim[[i]][,(k-1)*n.var + 1:n.var ] )
#     return(out)
#     } )
#   
#   # Make the full set of bootstrapped IRFs
#   l.irfs <-
#     lapply( 1:n.boot, function(i){
#         message('n=', i)
#         irf_var_fcast( l.struct.id[[i]], var.names=l.est$args$var.names, n.pds=n.pds, unit=unit ) %>%
#         mutate( i.boot=i, type=as.factor(type) ) %>%
#         mutate_if(is.numeric, funs(ifelse(is.nan(.), 0, .)))
#     } )
#   df.irfs <- suppressMessages( l.irfs %>% reduce( full_join ) )
#   # Use the simulation function -> Estimation function -> IRF function to
#   # produce the distribution of IRFs
#   
#   # Summarize the IRFs by percentile
#   df.irf.pc <- suppressMessages(
#     lapply( pc, function(x)
#       df.irfs %>% group_by( type, dimension, pd ) %>%
#         summarise_if( is.numeric, quantile, probs=x, na.rm=T ) %>%
#         mutate( pc=paste0('pc.', x ) )
#     ) %>%
#       reduce(full_join) %>%
#       full_join( irf_var_fcast( l.est, var.names=l.est$args$var.names, n.pds=n.pds, unit=unit ) %>%
#                    mutate( type=as.factor(type) ) %>%
#                    mutate_if(is.numeric, funs(ifelse(is.nan(.), 0, .))) %>%
#                    mutate( pc='pt.est' ) )
#   )
#   
#   return( df.irf.pc )
# }
# 
# make_irf_labels <- function(df.irf.in, var.names, fcast.names, sign.var='ir03', 
#                             shk.type.only=FALSE ){
# # Re-signs and labels the IRFs
#   
#   browser()
# 
#   # Compute the signs of the interest rate IRFs and correct
#   df.irf.sign <- df.irf.in %>% ungroup %>%
#     filter(pd==1, type=='u', pc=='pt.est') %>%
#     mutate( irf.sign=sign(get(sign.var))) %>%
#     select( dimension, irf.sign  )
#       # The sign of the sign variable
#   df.irf.rescaled <- df.irf.in %>%
#     ungroup() %>%
#     left_join( df.irf.sign ) %>%
#     mutate_at( c( var.names, fcast.names ), function(col) col * {.}$irf.sign )
#       # The re-signed IRFs
# 
#   # Now classify the shocks
#   df.irf.impact <- df.irf.rescaled %>%
#     filter(pd==1, type=='u', pc=='pt.est') %>%
#     mutate( shk.type=case_when( (infl<0 & y<0)~'Monetary policy',
#                                 (infl<0 & y>0)~'Supply',
#                                 (infl>0 & y>0)~'Demand',
#                                 (infl>0 & y>0)~'Phillips curve',
#                                 TRUE~'Unlabelled') ) %>%
#     select( dimension, shk.type, !!sign.var ) %>%
#     arrange(shk.type, !!sign.var ) %>%
#     group_by(shk.type) %>%
#     full_join(df.irf.rescaled %>%
#                 filter(pd==1, pc=='pt.est') %>% select( dimension, type) ) %>%
#     mutate( shk.name=paste0(shk.type, ' #', rep( ((n())/2):1, each=2),
#                             case_when((type=='epsilon'|shk.type.only)~'',
#                                       type=='u'~': \nContemp. shock',
#                                       type=='v'~': \nNews shock') )) %>%
#     ungroup() %>%
#     select( shk.name, type, dimension ) %>%
#     arrange( shk.name )
# 
#   df.out <- df.irf.rescaled %>%
#     left_join( df.irf.impact ) %>%
#     select(-irf.sign)
# 
#   return(df.out)
# }
# 
# make_irf_fiscal_labels <- function(df.irf.in, l.est, var.names, fcast.names, 
#                                    threshold, thresh.pds.x, fiscal.var='gov', 
#                                    horizon=4, sign.correction=1 ){
# # Labels the IRFs
# 
#   # Compute the signs of the interest rate IRFs and correct. Means can reapply
#   # the sign for the IRFs and not recompute all the bootstraps etc.
#   df.cum.by.shk <- df.irf.in %>% 
#     filter( pd <= horizon, type=='v' ) %>%
#     group_by( type, dimension ) %>%
#     summarise( gov.cum=sum(get(!!fiscal.var)), y.cum=sum(y) ) %>%
#     mutate( news.mult=y.cum/gov.cum ) %>%
#     ungroup() %>%
#     select(-type)
#       # The cumulative impulses and the multiplier
#   
#   df.var <- l.est$var.decomp$df %>%
#     filter(pd == horizon, outcome==fiscal.var, type=='epsilon') # %>%
#     # group_by(pd, dimension) %>%
#     # summarise( variance.ratio=sum(variance.ratio) ) %>%
#     # ungroup() %>%
#     # group_by(dimension) %>%
#     # summarise( var.avg=mean(variance.ratio), var.min=min(variance.ratio), 
#     #            var.min.thresh=nth(variance.ratio, 1+thresh.pds.x, order_by = variance.ratio ) )
#       # The share of the variance due to epsilon (adding u and v) by each dimension
#   
#   df.labels <- full_join( df.cum.by.shk, df.var, by='dimension') %>%
#     mutate( sign.crit = (sign(gov.cum)==sign.correction * sign(y.cum)),
#             var.crit = variance.ratio > threshold,
#             shk.name = ifelse( sign.crit&var.crit, 'Fiscal shock', 'Unlabelled shock') ) %>%
#     group_by( shk.name ) %>%
#     arrange( -variance.ratio ) %>%
#     mutate( shk.name=paste0( shk.name, ' #', 1:n(), 
#                              ifelse(sign.crit & var.crit, paste0(' (',round(100*variance.ratio),'%)'), '' ) ) ) %>%
#     ungroup() %>%
#     arrange(dimension)
#   
#   return(df.labels)
# }
# 
# make_epsilon_labels <- function(df.irf.in, l.est, var.names, fcast.names, 
#                                 threshold=0.1, mult.threshold=.5, fiscal.var='gov', horizon=4, 
#                                 sign.correction=1, mon.var='ir03' ){
# # Labels the shocks  using the epsilon IRF
#   
#   # Compute the signs of the interest rate IRFs and correct. Means can reapply
#   # the sign for the IRFs and not recompute all the bootstraps etc.
#   tax.in.data <- 'tax' %in% (df.irf.in %>% names)
#   
#   df.cum.by.shk <- df.irf.in %>% 
#     filter( pd <= horizon, type=='epsilon' ) %>%
#     group_by( dimension ) %>%
#     summarise( fiscal.cum=sum(get(!!fiscal.var)), y.cum=sum(y),
#                fiscal.1.4=mean(get(!!fiscal.var)[pd %in% 1:4]),
#                # fiscal.5.8=mean(get(!!fiscal.var)[pd %in% 5:8]),
#                # fiscal.9.12=mean(get(!!fiscal.var)[pd %in% 9:12]),
#                tax.1.4=if(tax.in.data) mean(tax[pd %in% 1:4]) else NA,
#                # tax.5.8=if(tax.in.data) mean(tax[pd %in% 5:8]) else NA,
#                # tax.9.12=if(tax.in.data) mean(tax[pd %in% 9:12]) else NA,
#                mon.1.4=mean(get(!!mon.var)[pd %in% 1:4]),
#                inf.1.8=mean(infl[pd %in% 1:8]),
#                y.1.4=mean(y[pd %in% 1:4]),
#                y.1.8=mean(y[pd %in% 1:8]) ) %>%
#     mutate( mult=y.cum/fiscal.cum ) %>%
#     ungroup() %>%
#     mutate( fiscal.sign.crit=(sign.correction*sign(fiscal.1.4)==sign(y.1.4)),
#             consolidation.crit=(sign.correction*sign(fiscal.1.4)==-sign(tax.1.4)),
#             fiscal.mult.crit = (mult>mult.threshold),
#               # &(sign.correction*sign(fiscal.5.8)==-sign(tax.5.8))&
#               # (sign.correction*sign(fiscal.9.12)==-sign(tax.9.12)),
#             demand.crit=(sign(y.1.8)==sign(inf.1.8))&(sign(mon.1.4)==sign(inf.1.8)),
#             supply.crit=(sign(y.1.8)==-sign(inf.1.8)),
#             mon.crit=(sign(y.1.4)==sign(inf.1.8))&(sign(mon.1.4)==-sign(inf.1.8)) )
#       # The cumulative impulses and the multiplier
#   
#   df.var <- l.est$var.decomp$df %>%
#     filter(pd == horizon, outcome %in% c(fiscal.var, mon.var), type=='epsilon') %>% 
#     select(-pd, -type) %>%
#     mutate( outcome=case_when(outcome==fiscal.var~'var.fiscal',
#                               outcome==mon.var~'var.monetary',
#                               TRUE~outcome ) ) %>%
#     spread( outcome, variance.ratio ) %>%
#     mutate( fiscal.var.crit = (var.fiscal>threshold) )
#       # Just for information.  No need for labelling.  
#   
#   df.labels <- full_join( df.cum.by.shk, df.var, by='dimension') %>%
#     mutate( shk.name = case_when( (fiscal.sign.crit & fiscal.var.crit & fiscal.mult.crit & !(consolidation.crit))~' Government spending shock',
#                                   (fiscal.sign.crit & fiscal.var.crit & consolidation.crit)~' Fiscal stimulus shock',
#                                   demand.crit~'Demand shock',
#                                   mon.crit~' Monetary policy shock',
#                                   supply.crit~'Supply shock',
#                                   TRUE~'Unlabeled shock') ) %>%
#     group_by( shk.name ) %>%
#     mutate( irf.sign = case_when( shk.name==' Government spending shock'~sign(fiscal.1.4),
#                                   shk.name==' Fiscal stimulus shock'~-sign(tax.1.4),
#                                   shk.name=='Demand shock'~sign(y.1.8),
#                                   shk.name=='Supply shock'~sign(y.1.8),
#                                   shk.name==' Monetary policy shock'~sign(y.1.8),
#                                   TRUE~1 ) ) %>%
#     mutate( shk.name=paste0( shk.name, if(n()>1) paste0(' #', 1:n() ) else NULL ) ) %>%
#     ungroup() %>%
#     arrange(dimension)
#   
#   return(df.labels)
# }
# 
# sd_trim <- function(x, trim=.2, ... ){
# # Computes the standard deviation after dropping outliers, and corrects for normality
#   x.trim <- x[ ( x>quantile(x,trim/2,...) & x<quantile(x,1-trim/2,...) ) ]
#     # The trimmed dataset
#   set.seed(42)
#   # x.sim <- rnorm(1e6)
#   # y.sim <- x.sim[ x.sim>quantile(x.sim,trim/2) & x.sim<quantile(x.sim,1-trim/2) ]
#   # correction <- 1/sd(y.sim)
#   correction <- if(trim==.1) 1.266571 else if(trim==.2) 1.51107
#     # Pre-computed using the exercise above
#   return( correction * sd(x.trim, ... ) )
# }
