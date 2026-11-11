##########################################################################################
#
# validate_via_sim_2025.R
#
# Code to validate the method for the "Identification by RE" project using
# a simulated data set. This works. It recovers the true coefficients from the
# simulated data. See "CONCLUSIONS" comments for more detail.
#
# Philip Barrett, Washington DC
# First version: 02dec2022
#
##########################################################################################


#### 0. Set up ####
rm(list=ls())
library(tidyverse)
library(data.table)
library(lubridate)
# library(zoo)
source('code/var_sim_fns.R')
source('code/make_rf.R')
source('code/bootstrap.R')
source('code/structural_identification.R')
source('code/plotting.R')
select <- dplyr::select
v.shk.types <- c('contemp', 'news', 'epsilon' ) 


#### 1. Controls ####
url.base <- 'https://raw.githubusercontent.com/jonathanjadams/news-shocks-from-forecasts/refs/heads/main/simulated-data'
run.name <- 'NK3_monetary_news_6shocks_baseline' 
# run.name <- 'NK3_monetary_news_6shocks_var_u_0.95' # 
# run.name <- 'NK3_monetary_news_9shocks' #  
n.horiz.dsge <- 2
url.prefix <- paste0(url.base, '/', run.name, '/')
m.A.file <- paste0( url.prefix, 'Amatrix.csv' )
m.C.file <- paste0( url.prefix, 'Cmatrix.csv' )
m.D.file <- paste0( url.prefix, 'parametervector.csv' )
m.var.decomp.file <- paste0( url.prefix, 'variance_decomposition.csv' )
sim.data.file <- paste0( url.prefix, 'simdata.csv' )
sim.shocks.file <- paste0( url.prefix, 'simshocks.csv' )

# m.D.u.2 <- c(.5, 2/3, 3/4) %>% diag
# m.D.v.2 <- diag( m.D.u.2 %>% nrow ) - m.D.u.2
data.names <- c( 'y', 'pi', 'i' ) # , 'gov' )
fcast.names <- paste0(data.names, '_f')
shk.names <- c( 'x', 'z', 'm' )
df.shk.cht.names <- data.frame( shk=c( shk.names, ' Total'), 
    shk.cht=c( 'Cost-push, x', 'Demand, z', 'Mon. Pol., h', ' Total' ) )
sign.var <- 'i'
n.burn <- 1000
n.lags <- 1
n.boot <- 100
n.vars <- data.names %>% length()
n.pds.irf <- 15
n.match.pds <- n.pds.irf
df.names <- data.frame( dimension=1:3, shk.name=shk.names )
l.irf.file <- lapply( data.names, function(x) paste0( url.prefix, 'irf_', x, '_allshocks.csv') ) %>%
  set_names(data.names)
v.long.shk.names <-  paste0( letters[ c(1:(n.horiz.dsge-1) + 21, 21) ], 
                             '_', rep(shk.names, each=n.horiz.dsge) )

#### 2. Get the simulation data ####

## 2.1 Read in the model coefficients ##
m.A <- read_csv(m.A.file, col_names = shk.names) %>%
  as.matrix %>%
  set_rownames(data.names)
m.C <- read_csv(m.C.file, col_names = shk.names) %>%
  as.matrix %>%
  set_rownames(data.names)
v.D.in <- read_csv(m.D.file, col_names = paste0( c('v_', 'u_'), rep(shk.names, each=2) ) )
m.D.u.2 <- (v.D.in %>% select(contains('u_')) %>% diag) ^2 %>%
  set_rownames( paste0('u_', shk.names) ) %>%
  set_colnames( paste0('u_', shk.names) )
m.D.v.2 <- (v.D.in %>% select(contains('v_')) %>% diag) ^2 %>%
  set_rownames( paste0('v_', shk.names) ) %>%
  set_colnames( paste0('v_', shk.names) )
m.var.decomp <- read_csv(m.var.decomp.file, 
                         col_names = v.long.shk.names ) %>%
  as.matrix %>%
  set_rownames(data.names)
n.vars <- data.names %>% length

## 2.2 Read in the simulation ##
df.sim.data <- read_csv(sim.data.file, col_names = c(fcast.names, data.names) ) %>%
  mutate( pd=1:n()) %>%
  select(pd, everything())
dt.shock.data.in <- 
  setDT( read_csv(sim.shocks.file, col_names = paste0( c('v_', 'u_'), rep(shk.names, each=2) ) ) )[, pd := .I]
long <- melt(dt.shock.data.in, id.vars = "pd", variable.name = "variable", value.name = "value")
long[, c("type", "shock") := tstrsplit(variable, "_", fixed = TRUE)]
wide <- dcast(long, pd + shock ~ type, value.var = "value")
wide[, epsilon := u + shift(v), by = shock]
  # Use data.table functions for such a large dataset
df.shock.data <- wide %>% 
  gather(type, value, -shock, -pd ) %>%
  mutate( variable=paste0(type,'_',shock)) %>%
  select(-type, -shock) %>%
  spread(variable, value) %>%
  as_tibble
df.sim <- cbind( df.sim.data, df.shock.data %>% select(-pd) ) %>%
  as_tibble

## 2.3 Read in the IRFs  ##
shk.idx <- seq(along=shk.names) %>%
  set_names(shk.names)
df.irfs.exact <- lapply( data.names, function(x) 
  read_csv( l.irf.file[[x]], col_names = paste0( c('v_', 'u_'), rep(shk.names, each=2) ) ) %>%
    mutate(pd=1:n()) %>% 
    gather( shk, value, -pd ) %>%
    separate(shk, c('type', 'shk'), sep='_') %>%
    mutate( dimension=shk.idx[shk] ) %>%
    select( pd, type, dimension, value) %>%
    group_by( type, dimension ) %>%
    mutate(fcast=lead(value)) %>%
    rename( !!x:=value, !!paste0(x, '_f'):=fcast ) ) %>%
  reduce( full_join, by=c('pd', 'type', 'dimension') )

# m.B <- read_csv(m.B.file, col_names = paste0(c(fcast.names,data.names),'.l1') ) %>%
#   as.matrix %>%
#   set_rownames(c(fcast.names,data.names))
# # l.B.j <- list( m.B )
# l.B.j <- lapply( 1:n.lags, function(i) m.B[ n.lags*(i-1) + 1:n.vars, 1:n.vars] )
# big.A <- rbind( cbind( l.B.j[[1]] %*% m.C + m.A, l.B.j[[1]] %*% m.A ),
#                 cbind(                      m.C,                m.A ) )
# m.big.Sigma <- big.A %*% diag( c( m.D.v.2 %>% diag, m.D.u.2 %>% diag ) ) %*% t(big.A)

# ## 2.2 Simulate the model ##
# x.0 <- matrix(0,nrow=1, ncol=n.vars)
# df.sim <- sim_var_fcast( x.0, l.B.j, m.A, m.C, m.D.u.2, m.D.v.2, 
#                          n.sim + n.burn ) %>%
#   tail(n.sim) %>%
#   set_names( c( 'pd', data.names, fcast.names, 
#                 paste( rep( c('v', 'u', 'epsilon' ), each=n.vars), 
#                        rep(shk.names,3), sep='_' ) ) )
# 



#### 3. Recover the model via regression ####

## 3.1 Directly from the DGP: Eq (6), including the errors ##
st.reg.rhs <- paste0( ' ~ 0 + ', paste0( 'lag(', data.names, ')' ) %>% paste( collapse=' + '), ' + ',
                      paste0( 'epsilon_', shk.names ) %>% paste( collapse=' + '), ' + ',
                      paste0( 'v_', shk.names ) %>% paste( collapse=' + ') )
l.lm.reg <- lapply( data.names, function(x) lm( paste0( x, st.reg.rhs ) %>% as.formula(), df.sim ) )
m.A.reg.1 <- sapply( l.lm.reg, function(x) (x %>% coef)[ x %>% coef %>% names %>% grepl('^epsilon_',.) ] ) %>% t()
m.B.reg.1 <- sapply( l.lm.reg, function(x) (x %>% coef)[ x %>% coef %>% names %>% grepl('^lag\\(',.) ] ) %>% t()
m.C.reg.1 <- sapply( l.lm.reg, function(x) (x %>% coef)[ x %>% coef %>% names %>% grepl('^v_',.) ] ) %>% t()
D.u.2.reg.1 <- df.sim %>% select( paste0('u_', shk.names) ) %>% as.matrix() %>% cov() %>% diag
D.v.2.reg.1 <- df.sim %>% select( paste0('v_', shk.names) ) %>% as.matrix() %>% cov() %>% diag
big.A.reg.1 <- rbind( cbind( m.B.reg.1 %*% m.C.reg.1 + m.A.reg.1, 
                             m.B.reg.1 %*% m.A.reg.1),
                      cbind( m.C.reg.1, m.A.reg.1 ) )
m.big.Sigma.reg.1 <- big.A.reg.1 %*% diag( c( D.v.2.reg.1, D.u.2.reg.1 ) ) %*% t( big.A.reg.1 )
r.sq.reg.1 <- sapply( l.lm.reg, function(x) x %>% summary() %>% .[['r.squared']] )
# The R2
A.err.1 <- (m.A.reg.1 - m.A) %>% abs %>% max
C.err.1 <- (m.C.reg.1 - m.C) %>% abs %>% max
  # CONCLUSION: The regression on the DGP with forecast errors recovers the true
  # values for A and B, and C. Thus, the data is being generated ok.
D.v.2.err.1 <- (m.D.v.2 %>% diag - D.v.2.reg.1) %>% abs %>% max
D.u.2.err.1 <- (m.D.u.2 %>% diag - D.u.2.reg.1) %>% abs %>% max
  # CONCLUSION: The variance estimate is much harder, even with 1 million
  # observations we have errors at the 3rd decimal place for D.v.2, D.u.2. I
  # guess that variances are just hard to estimate? This then feeds through to
  # larger errors on big.Sigma.
l.B.j <- list( m.B.reg.1 )

## 3.2 Check structural identification with the correct B, Sigma ##

# 3.2.1 Create the IRF matching #
df.pt.obj.irf <- irf_var_fcast( list( A=m.A, l.B.j=list(m.B.reg.1), C=m.C, D.u.2=m.D.u.2, D.v.2=m.D.v.2),
                             n.pds = n.pds.irf, var.names = data.names )
loss.target <- df.pt.obj.irf %>%
  arrange( dimension, type, pd ) %>%
  select(dimension, type, pd, data.names) %>%
  rename( data.names %>% set_names(paste0('x.', 1:n.vars)) ) %>%
  gather( variable, value, -dimension, -type, -pd ) %>%
  spread( dimension, value ) %>%
  select( -type, -pd, -variable ) %>%
  as.matrix

loss.fn <- function( cand.irf ){
  err <- ( cand.irf - loss.target )
  return( err^2 %>% colSums() )
}

# 3.2.2 Apply the identification restrictions: Near-exact B, Sigma #
# These diagnostics require positive news variance for every shock.
has.news.variance <- all(diag(m.D.v.2) > 0)
if(n.horiz.dsge==2 && has.news.variance){
  m.big.A <- rbind( cbind( m.B.reg.1 %*% m.C + m.A, 
                           m.B.reg.1 %*% m.A),
                        cbind( m.C, m.A ) )
  m.big.Sigma <- m.big.A %*% diag( c( diag(m.D.v.2), diag(m.D.u.2) ) ) %*% t( m.big.A )
  struct.reg.x <- structural_identification( list(m.B.reg.1), m.big.Sigma,
                                             ordering.loss.fn = loss.fn, 
                                             n.match.pds=n.match.pds )
  A.err.struct.x <- (struct.reg.x$A - m.A) %>% abs %>% max
  C.err.struct.x <- (struct.reg.x$C - m.C) %>% abs %>% max
  D.v.2.err.struct.x <- (struct.reg.x$D.v.2 - m.D.v.2) %>% abs %>% max
  D.u.2.err.struct.x <- (struct.reg.x$D.u.2 - m.D.u.2) %>% abs %>% max
    ## CONCLUSION: The identifying restrictions are *super* sensitive to
    ## inconsistent B, Sigma
  
  # 3.2.3 Apply the identification restrictions: Big regression B, Sigma #
  struct.reg <- structural_identification( list(m.B.reg.1), m.big.Sigma.reg.1,
                                           ordering.loss.fn = loss.fn, n.match.pds=n.match.pds )
  A.err.struct.1 <- (struct.reg$A - m.A) %>% abs %>% max
  C.err.struct.1 <- (struct.reg$C - m.C) %>% abs %>% max
  D.v.2.err.struct.1 <- (struct.reg$D.v.2 - m.D.v.2) %>% abs %>% max
  D.u.2.err.struct.1 <- (struct.reg$D.u.2 - m.D.u.2) %>% abs %>% max
    ## CONCLUSION: The structural identification method recovers the correct
    ## parameters when using the B, Sigma recovered from regression including the
    ## shocks.
  struct.reg.alt <- structural_identification( list(m.B.reg.1), m.big.Sigma.reg.1 + diag(6) * 1e-04,
                                           ordering.loss.fn = loss.fn, n.match.pds=n.match.pds )
    ## CONCLUSION: Small errors in estimating Sigma do not lead to big changes in
    ## the structural answers.  This looks good
}


## 3.3 Now excluding the errors: Eq (7) upper and lower ##
if(has.news.variance){
  st.reg.rhs.f.x <- paste0( ' ~ 0 + ', paste0( 'lag(',fcast.names, ')' ) %>% paste( collapse=' + ') )
  l.lm.reg.f.x <- lapply( fcast.names, function(x) lm( paste0( x, st.reg.rhs.f.x ) %>% as.formula(), df.sim ) )
  m.B.reg.f.x <- sapply( l.lm.reg.f.x, function(x) (x %>% coef)[ x %>% coef %>% names %>% grepl('^lag',.) ] ) %>%
    t() %>% set_rownames(fcast.names)
  err.hat.1.reg.x <- sapply( l.lm.reg.f.x, resid )
  err.hat.2.reg.x <- sapply( data.names, function(this.var)
    df.sim %>% mutate( err=eval(parse( text = paste0(this.var,
                                                     '-lag( ', this.var, '_f )' ) ) ) ) %>% pull(err) )[-1,]
  big.Sigma.reg.f.x <- cov( cbind( err.hat.1.reg.x, err.hat.2.reg.x ) )
  r.sq.f.x <- sapply( l.lm.reg.f.x, function(x) x %>% summary() %>% .[['r.squared']] )
      # NB: R squared is <1 because the errors are excluded from the regression
  B.err.f.x <- (m.B.reg.f.x - m.B.reg.1) %>% abs %>% max
  big.Sigma.err.f.x <- (big.Sigma.reg.f.x - m.big.Sigma.reg.1) %>% abs %>% max
      # CONCLUSION: We get some small error on the estimates for B and Sigma
  struct.f.x <- structural_identification( list(m.B.reg.f.x), big.Sigma.reg.f.x,
                                           ordering.loss.fn = loss.fn, n.match.pds=n.match.pds )
      # Apply the structural identification transformation
  A.err.f.x <- (struct.f.x$A - m.A) %>% abs %>% max
  B.err.f.x <- (m.B.reg.f.x - m.B.reg.1) %>% abs %>% max
  C.err.f.x <- (struct.f.x$C - m.C) %>% abs %>% max
      # CONCLUSION: Errors on the structural parameters are not nothing
}


#### 4. Apply our method ####

l.rf.ols.sim <- make.rf.var( df.sim %>% select( pd, data.names, fcast.names) %>% na.omit, data.names, fcast.names, pd.name='pd', 
                             lags = n.lags, restrict=TRUE, do.aic=FALSE )
l.rf.adj.sim <- list(rf=l.rf.ols.sim)
# l.rf.adj.sim <- rf_debias( l.rf.ols.sim, method='analytic' )
l.rf.adj.sim$rf$df <- l.rf.ols.sim$df %>%
  rename(date=pd)
# Adds the data as a data frame to the debiased RF 
l.struct.sim <- struct_id_rf( l.rf.adj.sim, sign.var, n.pds.irf-1, pt.est.only = TRUE )
l.order <- shk_sign_order( loss.fn, l.struct.sim$pt$A, l.struct.sim$pt$C, l.struct.sim$pt$D.v.2, 
                           l.struct.sim$pt$D.u.2, l.struct.sim$pt$l.B, n.match.pds )

message( 'Maximum eigenvalues' )
m.eig <- sapply( list( restricted=l.struct.sim$pt ), function(x) 
  sapply( list( x.form=x$big.B.x, f.form=x$big.B.f), function(y) ( y %>% eigen() )$values %>% abs %>% max) )
m.eig %>% print()
# Check the eigenvalues
B.err.est <- (l.struct.sim$pt$B - l.B.j[[1]]) %>% abs %>% max
big.Sigma.err.est <- (l.struct.sim$pt$Sigma - m.big.Sigma.reg.1) %>% abs %>% max
m.A.est <- (l.struct.sim$pt$A %*% l.order$sign)[,l.order$order]
m.C.est <- (l.struct.sim$pt$C %*% l.order$sign)[,l.order$order]
m.D.u.2.est <- l.struct.sim$pt$D.u.2[l.order$order,l.order$order]
m.D.v.2.est <- l.struct.sim$pt$D.v.2[l.order$order,l.order$order]
A.err.est <- (m.A.est - m.A) %>% abs %>% max
C.err.est <- (m.C.est - m.C) %>% abs %>% max
D.u.2.err.est <- (m.D.u.2.est - m.D.u.2) %>% abs %>% max
D.v.2.err.est <- (m.D.v.2.est - m.D.v.2) %>% abs %>% max


#### 5. Graphs #####

## 5.0 IRFs ##
l.est.reordered <- l.struct.sim$pt
l.est.reordered$A <- m.A.est
l.est.reordered$C <- m.C.est
l.est.reordered$D.u.2 <- m.D.u.2.est
l.est.reordered$D.v.2 <- m.D.v.2.est
df.pt.est.irf <- irf_var_fcast( l.est.reordered,
                             n.pds = n.pds.irf, var.names = data.names ) %>%
  left_join(df.names, by='dimension') %>% 
  mutate(irf.sign=1)

gg.irf.est <- irf_plot_pt(df.pt.est.irf, data.names, fcast.names)
irf_plot_pt_overlay( gg.irf.est,  df.irfs.exact %>% 
                       left_join(df.names, by='dimension') %>% mutate(irf.sign=1), 
                     data.names, fcast.names ) %>% print()

irf_plot_pt_overlay( gg.irf.est,  df.pt.obj.irf %>% 
                       left_join(df.names, by='dimension') %>% mutate(irf.sign=1), 
                     data.names, fcast.names ) %>% print()


if(FALSE){
  
  ## 5.1 Big Sigma ##
  plot( log(abs(m.big.Sigma), 10),  log( abs(big.Sigma.reg.f.x - m.big.Sigma ), 10 ), 
        main='Log absolute error: Big Sigma, regression by hand', xlab='True value', 
        ylab='Error', ylim=c( -6.5, 0.5 ))
  abline(h=0:(-3), lty=1:4)
  abline(a=0, b=1, lwd=0.1 )
  plot( log(abs(m.big.Sigma), 10),  log( abs(l.struct.sim$pt$Sigma - m.big.Sigma ), base = 10 ), 
        main='Log absolute error: Big Sigma, main function', xlab='True value', 
        ylab='Error', ylim=c( -6.5, 0.5 ))
  abline(h=0:(-3), lty=1:4)
  abline(a=0, b=1, lwd=0.1 )
  
  ## 5.2 A & C ##
  plot( log(abs(m.A), 10),  log( abs(m.A.est - m.A ), 10 ), 
        main='Log absolute error: A, main function', xlab='True value', 
        ylab='Error', ylim=c( -4.5, 0.5 )  )
  abline(h=0:(-3), lty=1:4)
  abline(a=0, b=1, lwd=0.1 )
  
  plot( log(abs(m.C), 10),  log( abs(m.C.est - m.C ), 10 ),
        main='Log absolute error: C, main function', xlab='True value',
        ylab='Error', ylim=c( -4.5, 0.5 )  )
  abline(h=0:(-3), lty=1:4)
  abline(a=0, b=1, lwd=0.1 )
  
  ## 5.3 B ##
  plot( log(abs(l.B.j[[1]]), 10),  log( abs(m.B.reg.f.x - l.B.j[[1]] ), 10 ),
        main='Log absolute error: B, regression by hand', xlab='True value',
        ylab='Error', ylim=c( -4.5, 0.5 )  )
  abline(h=0:(-3), lty=1:4)
  abline(a=0, b=1, lwd=0.1 )
  
  plot( log(abs(l.B.j[[1]]), 10),  log( abs(l.struct.sim$pt$B - l.B.j[[1]] ), 10 ),
        main='Log absolute error: B, main function', xlab='True value',
        ylab='Error', ylim=c( -4.5, 0.5 )  )
  abline(h=0:(-3), lty=1:4)
  abline(a=0, b=1, lwd=0.1 )
}

## 5.4 Variance decomposition ##
l.var.decomp <- l.est.reordered %>% var.decomp( n.pds=20 )
df.var.decomp <- lapply( v.shk.types, 
   function(x) lapply( 1:n.pds.irf, 
                       function(i) data.frame( l.var.decomp$var.decomp[[i]][[x]] ) %>% 
       rownames_to_column('variable') %>%
       gather( dimension, variance, -variable ) %>%
       mutate( pd =i, type=x, dimension=gsub('X', '', dimension) %>% as.numeric ) ) %>%
     reduce(.,rbind) ) %>%
  reduce(.,rbind) %>% as_tibble()

df.var.decomp.share <- df.var.decomp %>% 
  left_join( {.} %>% group_by(pd, variable) %>% 
               summarise( total=sum(variance[type=='epsilon']) ) ) %>%
  mutate( variance.pc = 100 * variance / total ) %>% 
  left_join(df.names, by='dimension')

m.var.decomp.est <- (df.var.decomp.share %>%
  filter( pd==max(pd) ) %>% 
  select( variable, shk.name, type, variance.pc) %>% 
  filter( type!='epsilon') %>% 
  mutate( type=case_when(type=='contemp'~'u',
                         type=='news'~'v')) %>%
  unite(temp, type, shk.name) %>% 
  spread(temp, variance.pc) %>% 
  column_to_rownames('variable') %>% 
  as.matrix())[data.names,paste0( c('v_', 'u_'), rep(shk.names, each=2) ) ]


df.var.decomp.compare <- 
  full_join( m.var.decomp.est %>% 
               as.data.frame %>% 
               rownames_to_column('variable') %>% 
               gather(shock, stat, -variable),
             m.var.decomp %>% 
               as.data.frame %>% 
               rownames_to_column('variable') %>% 
               gather(shock, dsge, -variable),
             by=c('variable', 'shock') ) %>%
  gather( source, value, -variable, -shock) %>%
  separate(shock, c('type','shk'), sep='_' ) %>%
  full_join( {.} %>% group_by( source, type, variable ) %>%
               summarise(value=sum(value)) %>%
               mutate(shk=' Total') ) %>%
  mutate(type.cht=case_when(type=='u'~'Surprise',
                            type=='v'~'News',
                            type=='w'~'News +2'),
         source.cht=case_when(source=='stat'~'Identified VAR',
                              source=='dsge'~'DSGE solution') ) %>%
  left_join(df.shk.cht.names)

gg.var.decomp <- ggplot( df.var.decomp.compare, aes(x=source.cht, y=value, fill=type.cht ) ) +
  geom_col() +
  theme_minimal() +
  ylab('Percent of total variation') +
  xlab('') +
  facet_grid( variable ~ shk.cht, scales='free_y') +
  scale_x_discrete(guide=guide_axis(angle=45)) +
  labs(fill='') +
  theme(legend.position = 'bottom') +
  scale_fill_grey()

gg.var.decomp %>% print()

dir.create('figures/simulations', recursive=TRUE, showWarnings=FALSE)
ggsave( paste0('figures/simulations/var_decomp_', gsub( '\\.', '_', run.name), '.pdf'), 
        gg.var.decomp, height=4, width=4 )


