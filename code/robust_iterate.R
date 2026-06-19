##########################################################################################
#
# robust_iterate.R
#
# Code to iterate over possible VAR specifications in the baseline
#
# Philip Barrett, Washington DC
# First version: 29jan2023
#
##########################################################################################


#### 0. Set up ####
rm(list=ls())
setwd('C:\\Users\\pbarrett\\OneDrive - International Monetary Fund (PRD)\\archive/identification from RE/')
  # Current working directory.  Change as needed.
library(tidyverse)
library(lubridate)
library(RColorBrewer)
library(stringr)
library(xtable)
library(data.table) # Useful for rbindlist
# library(zoo)
source('code/make_rf.R')
source('code/bootstrap.R')
source('code/var_sim_fns.R')
source('code/structural_identification.R')
source('code/plotting.R') 
source('code/fcast_cleaning.R')
source('code/policy_counterfactuals.R')
select <- dplyr::select
# Prevent conflict with MASS::select

#### 1. Controls ####

## 1.1 Data ##
data.in.file <- 'data/final_replication_data.csv'
ml.fcast.file <- 'data/ml_fcast.rdata'
out_prefix <- '/final_replication/'
  # Set the output file
mult.100.vars <- c('y', 'infl', 'gov', 'sur', 'fed', 'tax' )
data.labels <- c( ir03='3-month \ninterest rate', 
                  ir5='5-year \ninterest rate', 
                  un='Unemployment', infl='Inflation', 
                  y=' Output', ir5='Five-year interest rate', 
                  oil='Oil Price', usdcad='USD-CAD exchange rate',
                  sur='Surplus (percent of GDP)', tax=' Taxes',
                  hs='Housing starts', indpro='Industrial production',
                  fed='  Federal government \nconsumption',
                  gov='  Government \nspending' )
sign.var <- 'y'

## 1.2 Specifications ##
start.pd <- '1969-01-01' # '1983-01-12' # '1950-01-01'    # Set later to restrict the sample
n.lag.ml <- 8           # The number of lags available for the ML forecast construction
lag.max <- 6            # The maximum number of lags for automatic lag selection.
n.pds.irf <- 25         # IRF length
n.boot <- 1000          # Bootstrap size
l.var.names <- list( 
  baseline_ml_1_lag=list( vars=c('infl', 'ir03', 'y', 'gov', 'tax', 'hs' ), do.cumul=TRUE,
                   do.cf=TRUE, lags=1, ml.clean=TRUE, make.ml=TRUE,  force.f=FALSE,
                   label.def=TRUE ),
  baseline_ml_1_lag_spf=list( vars=c('infl', 'ir03', 'y', 'gov', 'tax', 'hs' ), do.cumul=FALSE,
                          do.cf=FALSE, lags=1, ml.clean=TRUE, make.ml=TRUE,  force.f=FALSE,
                          label.def=FALSE, use.spf.i=TRUE, init='1981-03-01', final='2022-03-01' ),
  baseline_ml_1_lag_short=list( vars=c('infl', 'ir03', 'y', 'gov', 'tax', 'hs' ), do.cumul=FALSE,
                          do.cf=FALSE, lags=1, ml.clean=TRUE, make.ml=FALSE,  force.f=FALSE,
                          label.def=FALSE, init='1988-03-01', final='2019-12-01' ),
  baseline_ml_2_lag=list( vars=c('infl', 'ir03', 'y', 'gov', 'tax', 'hs' ), do.cumul=FALSE,
                                do.cf=FALSE, lags=2, ml.clean=TRUE, make.ml=FALSE,  force.f=FALSE,
                                label.def=FALSE ),
  baseline_ml_4_lag=list( vars=c('infl', 'ir03', 'y', 'gov', 'tax', 'hs' ), do.cumul=FALSE,
                          do.cf=FALSE, lags=4, ml.clean=TRUE, make.ml=FALSE,  force.f=FALSE,
                          label.def=FALSE ),
  baseline_ml_aic=list( vars=c('infl', 'ir03', 'y', 'gov', 'tax', 'hs' ),
                   do.cf=TRUE, lags=NA, ml.clean=TRUE, make.ml=FALSE,  force.f=FALSE,
                   label.def=FALSE, do.cumul=FALSE),
  baseline_f_1_lag=list( vars=c('infl', 'ir03', 'y', 'gov', 'tax', 'hs' ),
                        do.cf=TRUE, lags=1, ml.clean=TRUE, make.ml=FALSE,  force.f=TRUE,
                        label.def=FALSE, do.cumul=FALSE),
  fed_ml_1_lag=list( vars=c('infl', 'ir03', 'y', 'fed', 'tax', 'hs' ), do.cumul=TRUE,
                        do.cf=FALSE, lags=1, ml.clean=TRUE, make.ml=FALSE,  force.f=FALSE,
                   label.def=FALSE ),
  unemp_ml_1_lag=list( vars=c('infl', 'ir03', 'y', 'gov', 'tax', 'un' ), do.cumul=TRUE,
                        do.cf=FALSE, lags=1, ml.clean=TRUE, make.ml=FALSE,  force.f=FALSE,
                     label.def=FALSE, match.drop=c('hs', 'un') ),
  fiscal_ml_1_lag=list( vars=c('infl', 'ir03', 'y', 'gov', 'tax', 'sur' ), do.cumul=TRUE,
                     do.cf=FALSE, lags=1, ml.clean=TRUE, make.ml=FALSE,  force.f=FALSE,
                     label.def=FALSE, match.drop=c('hs', 'sur') ),
  ir5_ml_1_lag=list( vars=c('infl', 'ir03', 'y', 'gov', 'ir1', 'ir5' ), do.cumul=FALSE,
                      do.cf=FALSE, lags=1, ml.clean=TRUE, make.ml=FALSE,  force.f=FALSE,
                      label.def=FALSE, match.drop=c('hs', 'tax', 'ir5', 'ir1') ),
  indpro_ml_1_lag=list( vars=c('infl', 'ir03', 'y', 'gov', 'tax', 'indpro' ), do.cumul=FALSE,
                      do.cf=TRUE, lags=1, ml.clean=TRUE, make.ml=FALSE,  force.f=FALSE,
                      label.def=FALSE, match.drop=c('hs', 'indpro') ),
  usdcadfx_ml_1_lag=list( vars=c('infl', 'ir03', 'y', 'gov', 'tax', 'usdcad' ), do.cumul=FALSE,
                        do.cf=TRUE, lags=1, ml.clean=TRUE, make.ml=FALSE,  force.f=FALSE,
                        label.def=FALSE, match.drop=c('hs', 'usdcad') ),
  oil_ml_1_lag=list( vars=c('infl', 'ir03', 'y', 'gov', 'tax', 'oil' ), do.cumul=FALSE,
                        do.cf=FALSE, lags=1, ml.clean=TRUE, make.ml=FALSE,  force.f=FALSE,
                        label.def=FALSE, match.drop=c('hs', 'oil') )
  ) 
    # Specification details for different runs

## 1.3 Plotting/charting ##
p.vals <- c( .025, .05, .1, .25 )
  # p-values to compute.  Will also compute 1-minus these
var.pc <- c(.1, .9)
  # CI for the variance decomp table
pc.plot <- c( .1, .25, .75, .9 )
  # Plotting percentiles for IRFs etc.
l.mults <- list( lewis=list( h=2, mu.g=.56, mu.t=.03, label=' Lewis 2021' ),
                 lewis=list( h=4, mu.g=.57, mu.t=-.09, label=' Lewis 2021' ),
                 lewis=list( h=8, mu.g=.57, mu.t=-.71, label=' Lewis 2021' ),
                 lewis=list( h=12, mu.g=.64, mu.t=-1.33, label=' Lewis 2021' ),
                 lewis=list( h=16, mu.g=.76, mu.t=-1.77, label=' Lewis 2021' ),
                 lewis=list( h=20, mu.g=.87, mu.t=-2.06, label=' Lewis 2021' ),
                 bp=list( h=2, mu.g=.61, mu.t=-.64, label='Blanchard and Perotti 2002' ),
                 bp=list( h=4, mu.g=.60, mu.t=-1.01, label='Blanchard and Perotti 2002' ),
                 bp=list( h=8, mu.g=.56, mu.t=-2.28, label='Blanchard and Perotti 2002' ),
                 bp=list( h=12, mu.g=.60, mu.t=-3.63, label='Blanchard and Perotti 2002' ),
                 bp=list( h=16, mu.g=.70, mu.t=-4.69, label='Blanchard and Perotti 2002' ),
                 bp=list( h=20, mu.g=.80, mu.t=-5.41, label='Blanchard and Perotti 2002' ),
                 ramey=list( h=20, mu.g=1.2, mu.t=-2.6, label=' Ramey 2011/Romer and Romer 2010' ),
                 caldara.penalty=list( h=2, mu.g=.05, mu.t=-1.05, label='Caldara and Kamps 2017, penalty function' ),
                 caldara.penalty=list( h=6, mu.g=.35, mu.t=-1.2, label='Caldara and Kamps 2017, penalty function' ),
                 caldara.penalty=list( h=12, mu.g=0.55, mu.t=-0.8, label='Caldara and Kamps 2017, penalty function' ),
                 caldara.penalty=list( h=20, mu.g=0.25, mu.t=-0.45, label='Caldara and Kamps 2017, penalty function' ),
                 pappa=list( h=6, mu.g=2.4, mu.t=-1.25, label='Ben Zeev and Pappa 2017/Romer and Romer 2010' ),
                 ricco=list( h=16, mu.g=1.5, mu.t=-2.6, label='Ricco 2015/Romer and Romer 2010' ) )
  # Comparison multipliers for the fiscal shock


make.dir.if.needed <- function(dir.path) if(!file.exists(dir.path)) dir.create(dir.path,recursive = TRUE)
for( this.vars in names(l.var.names) ) make.dir.if.needed(paste0('figures',out_prefix, this.vars))

############### MAIN LOOP #################
for( this.vars in names(l.var.names) ){
  
  message( '************** Starting specification "', this.vars, '" *************************')
  
  
  #### 2. Read and organize in the data ####
  data.names <- l.var.names[[this.vars]]$vars %>% sort
  ml.clean <- l.var.names[[this.vars]]$ml.clean
  make.ml <- l.var.names[[this.vars]]$make.ml
  lags <- l.var.names[[this.vars]]$lags
  do.cf <- l.var.names[[this.vars]]$do.cf
  force.f <- l.var.names[[this.vars]]$force.f
  label.def <- l.var.names[[this.vars]]$label.def
  match.drop <- l.var.names[[this.vars]]$match.drop
  do.cumul <- l.var.names[[this.vars]]$do.cumul
  use.spf.i <- if(is.null(l.var.names[[this.vars]]$use.spf.i)) FALSE else l.var.names[[this.vars]]$use.spf.i
  fcast.names <- paste0(data.names, '_f')
  init <- if(is.null(l.var.names[[this.vars]]$init)) start.pd else l.var.names[[this.vars]]$init
  final <- if(is.null(l.var.names[[this.vars]]$final)) Inf else l.var.names[[this.vars]]$final
  fig.dir <- paste0('figures',out_prefix, this.vars)
  fiscal.var <- case_when( 'gov' %in% data.names~'gov',
                           'fed' %in% data.names~'fed',
                           'sur' %in% data.names~'sur' )
  message( ' Main fiscal target variable: ', fiscal.var )
  fiscal.sign <- if(fiscal.var %in% c('sur', 'tax')) -1 else 1
  
  
  ### 2.0 Read the raw data ###
  df.in <- suppressMessages( data.in.file %>% read_csv )
  df <- df.in %>% 
    separate( quarter, c('year', 'quarter'), sep = 'q') %>%
    mutate( quarter=quarter %>% as.numeric,
            date=make_date( year, 3*quarter, 1 ) ) %>%
    select( date, everything(), -year, -quarter ) %>% 
    rename_with(~str_remove(., '_dt')) %>%
    # select( date, data.names, fcast.names ) %>%
    filter( date >= init, date <= final ) %>%
    mutate_if( is.numeric, na.approx, rule=2 ) %>%
        # NB: Means we are filling in the ends with the nearest observations
    mutate_if( is.numeric, function(x) x-mean(x,na.rm=T) ) %>%
    mutate_at(vars(contains(mult.100.vars)), function(x) x*100 )
  if(use.spf.i) df <- df %>% mutate( ir03_f = ir03_spf )
    # Switch in the SPF interest rate forecast if required

  if( ml.clean & make.ml ){
  
    ### 2.1 Clean the forecasts ###
    fcast.clean <- fcast_ml_all( df, data.names, n.lag=n.lag.ml )
        # Need to summarize how much info is coming from the forecasts here.
    df.f.clean <- fcast.clean$df %>%
      set_names( c('date', paste0( data.names, '_z' ) ) )
    df.clean <- df %>%
      select( date, data.names, fcast.names ) %>%
      left_join( df.f.clean ) %>%
      na.omit()
    
    glmnet.raw.rmse.ratio <- (sapply( fcast.clean$details, function(x) x$rmse )['glmnet',] / sapply( fcast.clean$details, function(x) x$rmse )['raw',]) %>% mean
      
    ### 2.2 Some plots ###
    fcast.importances <- lapply(1:length(data.names), 
                                function(i) (fcast.clean$details[[i]]$model %>% 
                                               varImp() )$importance %>% 
                                  rownames_to_column('variable') %>%
                                  set_names( c( 'variable', data.names[i]) ) ) %>% 
      reduce(full_join) %>%
      gather( outcome, importance, -variable ) %>%
      na.omit() %>%
      filter(importance >0 ) %>%
      separate( variable, c('variable', 'lag'), sep = '\\.l' ) %>%
      mutate( lag=replace_na(as.numeric(lag),0),
              fcast=(str_sub(variable,-2,-1)=='_f')) %>%
      select(outcome, lag, variable, fcast, importance) %>%
      arrange(outcome, lag, variable) %>%
      group_by(outcome) %>%
      mutate( importance.total=100*importance/sum(importance) )
    ggplot( fcast.importances, aes(x=lag, y=importance.total, fill=variable ))+
      geom_bar(stat='identity') +
      facet_wrap(~outcome) +
      theme_minimal() +
      theme(legend.position = 'bottom')
    ggsave( paste0( fig.dir, '/fcast_importance_by_variable.pdf'), height=8, width=12)
    ggplot( fcast.importances, aes(x=lag, y=importance.total, fill=fcast ))+
      geom_bar(stat='identity') +
      facet_wrap(~outcome) +
      theme_minimal() +
      theme(legend.position = 'bottom')
    ggsave(paste0( fig.dir,'/fcast_importance_by_lag.pdf'), height=8, width=12)
    
    ### 2.3 Summary table ###
    n.coef.poss <- sapply( fcast.clean$details, function(x) x$model$coefnames %>% length )
    # Number of possible coefficients
    df.fcast.sumy <- 
      full_join( fcast.importances %>% group_by(outcome, fcast) %>%
                   summarise( fcast.importance=sum(importance.total) ) %>%
                   filter(fcast) %>% select(-fcast),
                 fcast.importances %>% mutate(lag.0=(lag==0)) %>%
                   group_by(outcome, lag.0) %>%
                   summarise( lag.0.importance=sum(importance.total) ) %>%
                   filter(lag.0) %>% select(-lag.0) ) %>%
      cbind( alpha=sapply( fcast.clean$details, function(x) x$model$bestTune$alpha ),
             n.coef.poss=n.coef.poss ) %>%
      full_join(fcast.importances %>% group_by(outcome) %>% summarise( n.coeff=n() )) %>% 
      mutate( pc.non.zero=n.coeff/n.coef.poss*100 ) %>%
      select( outcome, alpha, fcast.importance, lag.0.importance, pc.non.zero ) %>%
      replace_na(list(fcast.importance=0))
    
    
  }else if(ml.clean & !make.ml){
    load(ml.fcast.file)
      # Loads the data frame with all variables, and with both explicit and generated forecasts
    if( force.f ){
      df.clean <- df %>%
        select( date, data.names, fcast.names ) %>%
        left_join( df.f.clean %>% select( 'date', paste0( data.names, '_z' )) ) %>%
        na.omit() %>%
        mutate_if( is.numeric, function(x) x - mean(x) )
          # Forces the inclusion of the "_f" forecasts from the data and adds
          # the ML "_z" forecasts
    }else{
      df.clean <- df %>%
        select( date, data.names ) %>%
        left_join( df.f.clean %>% select( 'date', paste0( data.names, '_z' )) %>%
                     set_names( c('date', fcast.names ) ) ) %>%
        na.omit() %>%
        mutate_if( is.numeric, function(x) x - mean(x) )
          # Replaces the "_f" forecasts entirely with the ML "_z" types
    }
  }else{
    df.clean <- df %>%
      select( date, data.names, fcast.names )
  }
  

  #### 3. Estimate coefficients and compute IRFs ####
  
  ## 3.1 Automatic lag selection
  if(is.na(lags)){
    message('Lag selection...')
    l.multi.est <- lapply( c(restricted=T,unrestricted=F), 
       function(restrict){ message(' ** Restriction: ', restrict) ; 
         lapply(1:lag.max, function(n){ 
           message( '   * lag = ', n)
           l.rf <- make.rf.var( df.clean, data.names, fcast.names, pd.name='date',
                                lags = n, restrict=restrict ) %>%
             rf_debias(n.boot=n.boot, pt.est.only = TRUE )
           aic.pt <- l.rf$rf$aic
           # aic.boot <- sapply( l.rf$boot, function(x) x$aic )
           # aic.pc <- aic.boot %>% 
           #   apply( ., 1, quantile, prob=c(p.vals, .5, 1-p.vals) )
           return( list(pt=aic.pt) ) #, boot=aic.boot, pc=aic.pc) )
           } ) } )

    df.aic <- data.frame( lag=1:lag.max, t(sapply(l.multi.est$restricted, function(x) x$pt)), restrict=TRUE ) %>%
      full_join( data.frame( lag=1:lag.max, t(sapply(l.multi.est$unrestricted, function(x) x$pt)), restrict=FALSE ) ) %>%
      gather(variable, value, -lag, -restrict )
    ggplot( df.aic, aes(x=lag, y=value, color=restrict)) +
      geom_line(lwd=.9) + 
      facet_wrap( ~variable, scales = 'free_y') +
      theme_minimal()
    ggsave(paste0( fig.dir,'/aic.pdf'), height=8, width=12)
    lags <- df.aic %>% 
      filter(restrict, variable=='f.form') %>% 
      summarize( lag=lag[which.min(value)]) %>% 
      unlist()
    message('...complete')
  }
  message( 'Lag length for main VAR = ', lags )
  
  ## 3.2 Estimation
  message('  *** Restricted estimates ***')
  l.rf.ols <- make.rf.var( df.clean, data.names, fcast.names, pd.name='date', 
                           lags = lags, restrict=TRUE )
  l.rf.adj <- rf_debias( l.rf.ols, n.boot=n.boot )
  l.rf.adj$rf$df <- l.rf.ols$df               # Adds the data as a data frame to the debiased RF 
  l.struct <- struct_id_rf( l.rf.adj, sign.var, n.pds.irf-1 )
  data.idx <- colnames(l.struct$pt$data ) %in% data.names
  var.aic.sim <- sapply( l.struct$boot, function(x) (x$data[,data.idx] %>% VAR(lag.max = 12, ic = 'AIC'))$p ) %>% table()
    # Computes the number of lags which would be estimated by just a VAR.
  
  message('\n  *** Unrestricted estimates ***')
  l.rf.adj.alt <- make.rf.var( df.clean, data.names, fcast.names, pd.name='date', 
                               lags = lags, restrict=FALSE ) %>% 
    rf_debias( pt.est.only = TRUE, n.boot=n.boot )
  
  message( 'Maximum eigenvalues' )
  m.eig <- sapply( list( restricted=l.rf.adj$rf, unrestricted=l.rf.adj.alt$rf), function(x) 
    sapply( list( x.form=x$big.B.x, f.form=x$big.B.f), function(y) ( y %>% eigen() )$values %>% abs %>% max) )
  m.eig %>% print()
    # Show the eigenvalues

  ## 3.3 Compute the news IRFs
  message('Constructing IRFs: point estimates...')
  df.str.irf.pt.init <- irf_var_fcast( l.struct$pt, var.names=data.names, n.pds = n.pds.irf ) %>%
    mutate_if(is.numeric, funs(ifelse(is.nan(.), 0, .))) %>%
    mutate( type=as.factor(type), pc='pt.est' ) %>%
    select( pd, type, dimension, everything() )
      # The point estimates
  message('...complete')
  
  ## 3.4 Define the shock labels
  if( label.def ){
    message('**** WARNING: Hard-coded shock labels in use: ****')
    df.labels <- data.frame( dimension=1:6,
                           shk.name=c( 'Supply shock', ' Fiscal stimulus shock', 
                                       'Unlabeled shock #1',' Monetary policy shock',
                                       'Unlabeled shock #2', 'Demand shock' ),
                           irf.sign=c( 1, -1, 1, 1, 1, 1 ) )
  }else{
    load(paste0('figures', out_prefix, 'label_def.rdata'))
      # Load the target data
    if( 'fed' %in% data.names ){
      df.target <- df.target %>% mutate( variable=gsub('gov', 'fed', variable ) )
    }
    m.target <- df.target %>% filter(!(variable %in% match.drop) ) %>% 
      arrange(type, pd, variable) %>%
      select( -type, -pd, -variable ) %>% as.matrix()
      # The target matrix
    df.candidate <- df.str.irf.pt.init %>% 
      select( type, pd, dimension, data.names ) %>% 
      gather( variable, value, -type, -pd, -dimension ) %>% 
      spread( dimension, value )
    m.candidate <- df.candidate %>% filter(!(variable %in% match.drop) ) %>% 
      arrange(type, pd, variable) %>%
      select( -type, -pd, -variable ) %>% as.matrix()
      # The target and candidate IRFs
    if( !identical( df.target %>% filter(!(variable %in% match.drop)) %>% select(type, pd, variable),
               df.candidate %>% filter(!(variable %in% match.drop)) %>% select(type, pd, variable) ) ){
      stop('ERROR: VARABLE MISMATCH FOR SHOCK LABELING')
    }
    m.cost.pos <- sapply( 1:ncol(m.candidate), function(i) (m.candidate[,i] - m.target)^2 %>% colSums() )
    m.cost.neg <- sapply( 1:ncol(m.candidate), function(i) (m.candidate[,i] + m.target)^2 %>% colSums() )
    m.cost <- pmin(m.cost.pos, m.cost.neg )
    m.sign <- (m.cost.pos==m.cost) - (m.cost.neg==m.cost)
      # Set up the cost function for the allocation problem
    opt.order <- solve_LSAP(m.cost) %>% c
    opt.sign <- sapply( seq(along=opt.order), function(i) m.sign[i, opt.order[i]] )
    df.labels <- data.frame( dimension=opt.order, shk.name=df.labels.def$shk.name, 
                             irf.sign=df.labels.def$irf.sign*opt.sign ) %>%
      arrange(dimension)
  }
  
  print(df.labels %>% select( dimension, shk.name, irf.sign))
  df.str.irf.pt <- df.str.irf.pt.init %>%
    left_join(df.labels, by='dimension' )
      # Label the point IRFs
  l.pt.v.scale <- l.struct$pt ; l.pt.v.scale$D.v.2 <- l.pt.v.scale$D.u.2
  df.str.irf.pt.v.scale <- irf_var_fcast( l.pt.v.scale, var.names=data.names, 
                                          n.pds = n.pds.irf ) %>%
    mutate_if(is.numeric, funs(ifelse(is.nan(.), 0, .))) %>%
    mutate( type=as.factor(type), pc='pt.est' ) %>%
    select( pd, type, dimension, everything() )%>%
    left_join(df.labels, by='dimension' )
      # The rescaled news-shock IRFs
  if( label.def ){
    df.labels.def <- df.labels
    df.target <- df.str.irf.pt %>% 
      select( type, pd, dimension, data.names ) %>% 
      gather( variable, value, -type, -pd, -dimension ) %>% 
      spread( dimension, value )
        # Make the target for the dimension matching
    save( df.target, df.labels.def, file=paste0('figures', out_prefix, 'label_def.rdata') )
      # Rename and save the labeled IRFs for comparison for the next variant
  } 
  
  ## 3.5 The bootstrapped IRFs
  message('Constructing IRFs: bootstrap...')
  df.str.irf.boot <- 
    lapply( 1:n.boot, function(i) irf_var_fcast( l.struct$boot[[i]], 
                                                 var.names=data.names, n.pds = n.pds.irf ) %>%
    mutate_if(is.numeric, funs(ifelse(is.nan(.), 0, .))) %>%
    mutate( type=as.factor(type), i.boot=i ) ) %>%
    data.table::rbindlist() %>%
    select( i.boot, pd, type, dimension, everything() ) %>%
    left_join(df.labels, by='dimension' )
      # The bootstrapped estimates
  df.str.irf.boot.v.scale <- 
    lapply( 1:n.boot, function(i){
      this.boot <- l.struct$boot[[i]]
      this.boot$D.v.2 <- this.boot$D.u.2
      return( irf_var_fcast( this.boot, 
                             var.names=data.names, n.pds = n.pds.irf ) %>%
                mutate_if(is.numeric, funs(ifelse(is.nan(.), 0, .))) %>%
                mutate( type=as.factor(type), i.boot=i ) )
      } ) %>%
    data.table::rbindlist() %>%
    select( i.boot, pd, type, dimension, everything() ) %>%
    left_join(df.labels, by='dimension' )
    # The bootstrapped estimates: with the news shocks rescaled to match the
    # surprises.
  message('...complete')
  message('Constructing IRFs: bootstrap percentiles...')
  df.str.irf.pc <- 
    lapply( c( p.vals, .5, rev(1-p.vals ) ), function(x)
      df.str.irf.boot %>% group_by( type, dimension, shk.name, pd ) %>%
        summarise_if( is.numeric, quantile, probs=x, na.rm=T ) %>%
        mutate( pc=paste0('pc.', x ) ) %>%
        select(-i.boot)
    ) %>%
    reduce(rbind) %>%
    select( type, dimension, pd, pc, everything() )
      # Bootstrap percentiles
  df.str.irf.pc.v.scale <- 
    lapply( c( p.vals, .5, rev(1-p.vals ) ), function(x)
      df.str.irf.boot.v.scale %>% group_by( type, dimension, shk.name, pd ) %>%
        summarise_if( is.numeric, quantile, probs=x, na.rm=T ) %>%
        mutate( pc=paste0('pc.', x ) ) %>%
        select(-i.boot)
    ) %>%
    reduce(rbind) %>%
    select( type, dimension, pd, pc, everything() )
  # Bootstrap percentiles
  
  message('...complete')
  
  ## 3.5 Create the epsilon shock IRFs ##
  ## NB: Differs from the news IRFs by timing and shock size
  message('Creating epsilon IRFs...')
  df.eps.irf.pt <- df.str.irf.pt %>% irf_var_epsilon
  df.eps.irf.boot <- suppressMessages( df.str.irf.boot %>% 
                                         group_by(i.boot) %>% 
                                         group_modify( ~irf_var_epsilon(.x) ) %>%
                                         ungroup() )
  df.eps.irf.pc <- lapply( c( p.vals, .5, rev(1-p.vals ) ), function(x)
    df.eps.irf.boot %>% group_by( type, dimension, shk.name, pd ) %>%
      summarise_if( is.numeric, quantile, probs=x, na.rm=T ) %>%
      mutate( pc=paste0('pc.', x ) ) %>%
      select(-i.boot) ) %>%
    reduce(rbind) %>%
    select( type, dimension, pd, pc, everything() )
  message('...complete')
  
  ## 3.7 Variance decomposition ##
  message('Computing variance decomposition...')
  v.shk.types <- c('contemp', 'news', 'epsilon' ) ; names(v.shk.types)=v.shk.types
  l.var.decomp <- list( pt=l.struct$pt %>% var.decomp( n.pds=n.pds.irf ),
                        boot=lapply(l.struct$boot, var.decomp, n.pds=n.pds.irf ) )
  message('...complete')
  message('Calculating variance share percentiles...')
  df.var.decomp <- lapply( v.shk.types, function(x) lapply( 1:n.pds.irf, 
             function(i) data.frame( l.var.decomp$pt$var.decomp[[i]][[x]] ) %>% 
               rownames_to_column('variable') %>%
               gather( dimension, variance, -variable ) %>%
               mutate( pd =i, type=x, dimension=gsub('X', '', dimension) %>% as.numeric ) ) %>%
               reduce(.,rbind) ) %>%
      reduce(.,rbind) %>% as_tibble()
    # The point error variance decomposition dataframe
  df.var.decomp.boot <- lapply( 1:n.boot, function(i.boot)
    lapply( v.shk.types, function(x) lapply( 1:n.pds.irf, 
              function(i) data.frame( l.var.decomp$boot[[i.boot]]$var.decomp[[i]][[x]] ) %>% 
                rownames_to_column('variable') %>%
                gather( dimension, variance, -variable ) %>%
                mutate( pd =i, type=x, dimension=gsub('X', '', dimension) %>% as.numeric ) ) %>%
                reduce(.,rbind) ) %>%
    reduce(.,rbind) %>% as_tibble() %>% mutate('i.boot'=i.boot) ) %>%
    reduce(.,rbind) %>% as_tibble() 
  # The bootstrap error variance decomposition dataframe
  df.var.decomp.boot.frac <- df.var.decomp.boot %>% 
    left_join( {.} %>% group_by(pd, variable, i.boot) %>% 
                 summarise( total=sum(variance[type=='epsilon']) ) ) %>%
    mutate( variance.pc = 100 * variance / total ) %>%
    left_join( df.labels %>% select(dimension, shk.name ) %>% 
                 mutate( shk.name=str_remove(shk.name, ' shock' ) ) )
  df.var.decomp.boot.pc <- 
    lapply( c( p.vals, .5, rev(1-p.vals ) ), function(x)
      df.var.decomp.boot.frac %>% 
        full_join( {.} %>% group_by(variable, pd, type, i.boot ) %>% 
                     summarise( variance.pc=sum(variance.pc) ) %>% 
                     mutate( shk.name='Total') ) %>%
        group_by( type, shk.name, variable, pd ) %>%
        summarise_if( is.numeric, quantile, probs=x, na.rm=T ) %>%
        mutate( pc=paste0('pc.', x ) )
    ) %>%
    reduce(rbind) %>%
    select( type, shk.name, pd, pc, everything() )
      # Bootstrap percentiles
  message('...complete')
  
    
  ## 3.8 Fiscal Multipliers ##
  l.df.cumul <- l.df.cumul.all <- l.df.cumul.pc <- 
    l.df.cumul.all.pc <- l.df.cumul.mults <- l.df.cumul.mults.all <-  NULL
  if(do.cumul){
    for( this.shk in c( 'v', 'u', 'epsilon' )){
      this.gov <- if('gov' %in% data.names ) 'gov' else 'fed'
      l.df.cumul[[this.shk]] <- df.eps.irf.pt %>%
        filter(type==this.shk, 
               shk.name %in% c(' Government spending shock', ' Fiscal stimulus shock') ) %>%
        select( pd, dimension, shk.name, y, gov=this.gov, tax, irf.sign ) %>%
        group_by( dimension, shk.name ) %>%
        mutate( y.cum=cumsum(y) * irf.sign,
                gov.cum=cumsum(gov) * irf.sign, 
                tax.cum=cumsum(tax)  * irf.sign ) %>%
        ungroup()
      l.df.cumul.pc[[this.shk]] <- df.eps.irf.boot %>%
        filter(type==this.shk, 
               shk.name %in% c(' Government spending shock', ' Fiscal stimulus shock') ) %>%
        select( i.boot, pd, dimension, shk.name, y, gov=this.gov, tax, irf.sign ) %>%
        group_by( i.boot, dimension, shk.name ) %>%
        mutate( y.cum=cumsum(y) * irf.sign,
                gov.cum=cumsum(gov) * irf.sign, 
                tax.cum=cumsum(tax) * irf.sign ) %>%
        ungroup() %>%
        group_by( dimension, shk.name, pd ) %>%
        summarise( pc=c( p.vals, .5, rev(1-p.vals ) ),
                   y.cum=quantile( y.cum, c( p.vals, .5, rev(1-p.vals ) ) ),
                   gov.cum=quantile( gov.cum, c( p.vals, .5, rev(1-p.vals ) ) ),
                   tax.cum=quantile( tax.cum, c( p.vals, .5, rev(1-p.vals ) ) ) )
      
      l.df.cumul.mults[[this.shk]] <- lapply( l.mults, function(x) l.df.cumul[[this.shk]] %>%
                                   filter( pd==x$h ) %>% 
                                   mutate( tax.mult=x$mu.t, 
                                           gov.mult=x$mu.g, 
                                           y.cum=gov.mult*gov.cum + tax.mult * tax.cum, label=x$label) %>%
                                  select( label, pd, tax.mult, tax.cum, gov.mult, gov.cum, y.cum) ) %>%
        reduce( rbind )
    }
    
    df.cumul.tab <- l.df.cumul.mults[['u']] %>%
      select( Source=label, h=pd, `$\\mu_G^h$`=gov.mult, `$\\mu_T^h$`=tax.mult ) %>%
      group_by(Source) %>%
      mutate( Source=c(Source[1], rep('', n()-1) ) )
    print( xtable(df.cumul.tab, caption='Tax and spending multipliers from the literature', 
                  digits = c(0, 0, 0, 2, 2), label = 't:lit_multipliers' ), 
           include.rownames=FALSE, sanitize.text.function=identity, 
           floating=FALSE, hline.after=c(-1,-1,0,nrow(df.cumul.tab) ),
           file=paste0( fig.dir, '/lit_multipliers.tex') )
    
    df.cumul.tab.wide <- 
      l.df.cumul.mults[['u']] %>% 
      select( Source=label, pd, tax.mult, gov.mult ) %>%
      gather( variable, value, -Source, -pd ) %>%
      mutate( variable=case_when( 
        variable=='gov.mult'~'$\\mu_G^h$',
        variable=='tax.mult'~'$\\mu_T^h$') ) %>%
      spread(pd, value ) %>%
      group_by(Source) %>%
      mutate( Source=c(Source[1], rep('', n()-1) ) )
    print( xtable(df.cumul.tab.wide, caption='Tax and spending multipliers from the literature', 
                  digits = c(0, 0, 0, rep(2,7)), label = 't:lit_multipliers' ), 
           include.rownames=FALSE, sanitize.text.function=identity, 
           floating=FALSE, hline.after=c(-1,-1,0,nrow(df.cumul.tab.wide) ),
           file=paste0( fig.dir, '/lit_multipliers_wide.tex') )
    
    
    for( this.shk in c( 'v', 'u', 'epsilon' )){
      l.df.cumul.all[[this.shk]] <- df.eps.irf.pt %>%
        filter(type==this.shk ) %>%
        select( pd, dimension, shk.name, y, gov=this.gov, tax, irf.sign ) %>%
        group_by( dimension, shk.name ) %>%
        mutate( y.cum=cumsum(y) * irf.sign,
                gov.cum=cumsum(gov) * irf.sign, 
                tax.cum=cumsum(tax)  * irf.sign ) %>%
        ungroup()
      l.df.cumul.all.pc[[this.shk]] <- df.eps.irf.boot %>%
        filter(type==this.shk ) %>%
        select( i.boot, pd, dimension, shk.name, y, gov=this.gov, tax, irf.sign ) %>%
        group_by( i.boot, dimension, shk.name ) %>%
        mutate( y.cum=cumsum(y) * irf.sign,
                gov.cum=cumsum(gov) * irf.sign, 
                tax.cum=cumsum(tax) * irf.sign ) %>%
        ungroup() %>%
        group_by( dimension, shk.name, pd ) %>%
        summarise( pc=c( p.vals, .5, rev(1-p.vals ) ),
                   y.cum=quantile( y.cum, c( p.vals, .5, rev(1-p.vals ) ) ),
                   gov.cum=quantile( gov.cum, c( p.vals, .5, rev(1-p.vals ) ) ),
                   tax.cum=quantile( tax.cum, c( p.vals, .5, rev(1-p.vals ) ) ) )
      
      l.df.cumul.mults.all[[this.shk]] <- lapply( l.mults, function(x) l.df.cumul.all[[this.shk]] %>%
                                                filter( pd==x$h ) %>% 
                                                mutate( tax.mult=x$mu.t, 
                                                        gov.mult=x$mu.g, 
                                                        y.cum=gov.mult*gov.cum + tax.mult * tax.cum, label=x$label) %>%
                                                select( label, shk.name, pd, tax.mult, tax.cum, gov.mult, gov.cum, y.cum) ) %>%
        reduce( rbind )
    }
  }
  
  #### 4. Policy counterfactuals ####
  if( do.cf ){
    gov.shk.idx <- df.labels %>% filter( shk.name ==' Government spending shock') %>% pull(dimension) %>% c
    consol.shk.idx <- df.labels %>% filter( shk.name ==' Fiscal stimulus shock') %>% pull(dimension) %>% c
    mon.shk.idx <- df.labels %>% filter( shk.name ==' Monetary policy shock') %>% pull(dimension) %>% c
    
    l.cf <- list()
    
    v.F.base <- rep(0, data.names %>% length() ) %>% set_names(data.names)
      # The baseline version of the objective function
    sd.wts <- 1 / ( df.clean %>% select(data.names) %>% apply(2, sd) %>% c )
      # Objective function weights for joint targeting
    
    if( length(consol.shk.idx)>0 ){
      
      l.cf$consol$y <- cf_policy_pc( df.str.irf.pt,
                                  l.struct, df.labels, data.names, 
                                  consol.shk.idx, p.vals, 
                                  v.F=replace( v.F.base, data.names=='y', 1 ) )
      l.cf$consol$gov <- cf_policy_pc( df.str.irf.pt,
                                     l.struct, df.labels, data.names, 
                                     consol.shk.idx, p.vals, 
                                     v.F=replace( v.F.base, data.names=='gov', 1 ) )
      l.cf$consol$tax <- cf_policy_pc( df.str.irf.pt,
                                     l.struct, df.labels, data.names, 
                                     consol.shk.idx, p.vals, 
                                     v.F=replace( v.F.base, data.names=='tax', 1 ) )
      l.cf$consol$infl <- cf_policy_pc( df.str.irf.pt,
                                     l.struct, df.labels, data.names, 
                                     consol.shk.idx, p.vals, 
                                     v.F=replace( v.F.base, data.names=='infl', 1 ) )
      l.cf$consol$infl_y <- cf_policy_pc( df.str.irf.pt,
                                        l.struct, df.labels, data.names, 
                                        consol.shk.idx, p.vals, 
                                        v.F=replace( v.F.base, data.names %in% c('y','infl'),
                                                     sd.wts[data.names %in% c('y','infl')] ) )
    }
    
    if( length(mon.shk.idx)>0 ){
      
      l.cf$mon$y <- cf_policy_pc( df.str.irf.pt,
                                   l.struct, df.labels, data.names, 
                                   mon.shk.idx, p.vals, 
                                   v.F=replace( v.F.base, data.names=='y', 1 ) )
      l.cf$mon$infl <- cf_policy_pc( df.str.irf.pt,
                                  l.struct, df.labels, data.names, 
                                  mon.shk.idx, p.vals, 
                                  v.F=replace( v.F.base, data.names=='infl', 1 ) )
      l.cf$mon$ir03 <- cf_policy_pc( df.str.irf.pt,
                                  l.struct, df.labels, data.names, 
                                  mon.shk.idx, p.vals, 
                                  v.F=replace( v.F.base, data.names=='ir03', 1 ) )
      l.cf$mon$infl_y <- cf_policy_pc( df.str.irf.pt,
                                  l.struct, df.labels, data.names, 
                                  mon.shk.idx, p.vals, 
                                  v.F=replace( v.F.base, data.names %in% c('y','infl'),
                                               sd.wts[data.names %in% c('y','infl')] ) )
      
      # l.cf$mon$infl <- cf_policy_pc( l.df.irfs, l.struct, df.labels, data.names, 
      #                                   'infl', mon.shk.idx, p.vals, phi.j=c(0,0,0,0,0) )
      # l.cf$mon$taylor <- cf_policy_pc( l.df.irfs, l.struct, df.labels, data.names, 
      #                                  'ir03', mon.shk.idx, p.vals, phi.j=c(0,0,-1.5,0,-.5) )
    }
    
    if( (length(consol.shk.idx)>0) & (length(mon.shk.idx)>0) ){
      
      l.cf$both$y <- cf_policy_pc( df.str.irf.pt,
                                   l.struct, df.labels, data.names, 
                                     c( consol.shk.idx, mon.shk.idx), p.vals, 
                                     v.F=replace( v.F.base, data.names=='y', 1 ) )
      l.cf$both$infl <- cf_policy_pc( df.str.irf.pt,
                                   l.struct, df.labels, data.names, 
                                   c( consol.shk.idx, mon.shk.idx), p.vals, 
                                   v.F=replace( v.F.base, data.names=='infl', 1 ) )
      l.cf$both$infl_y <- cf_policy_pc( df.str.irf.pt, 
                                   l.struct, df.labels, data.names, 
                                   c( consol.shk.idx, mon.shk.idx), p.vals, 
                                   v.F=replace( v.F.base, data.names %in% c('y','infl'),
                                                sd.wts[data.names %in% c('y','infl')] ) )
    }
    
  }
  
  
  #### 5. Make charts ####
  
  
  ### 5.1 The Main IRFs ###
  gg.irf.pt <- irf_plot_pt(df.str.irf.pt, data.names, fcast.names, data.labels = data.labels )
  ggsave(paste0( fig.dir, '/irfs_all_pt.pdf'), plot=gg.irf.pt, 
         height=8, width=12)
    
  df.str.all <- rbind( df.str.irf.pt, df.str.irf.pc )
  gg.irf.news <- irf_plot_all(df.str.all, data.names, fcast.names, 
                                   data.labels, pc=pc.plot )
  ggsave(paste0( fig.dir, '/irfs_all.pdf'), 
         plot=gg.irf.news, height=8, width=12)
  
  gg.irf.news.pc <- irf_plot_all(df.str.irf.pc, data.names, fcast.names, 
                                  data.labels, pc=pc.plot )
  ggsave(paste0( fig.dir, '/irfs_all_pc.pdf'), 
         plot=gg.irf.news, height=8, width=12)
  
  gg.irf.surp.pc <- irf_plot_all(df.str.irf.pc %>% filter(type=='u'), data.names, fcast.names, 
                              data.labels, pc=pc.plot )
  ggsave(paste0( fig.dir, '/irfs_surp_all_pc.pdf'), 
         plot=gg.irf.surp.pc, height=8, width=12)
  ggsave(paste0( fig.dir, '/irfs_surp_all_pc_tall.pdf'), 
         plot=gg.irf.surp.pc, height=11, width=9)
  
  
  ### 5.1b With scaled news shocks ###
  df.str.all.v.scale <- rbind( df.str.irf.pt.v.scale, df.str.irf.pc.v.scale )
  gg.irf.news.v.scale <- irf_plot_all(df.str.irf.pc.v.scale, data.names, fcast.names, 
                                      data.labels, pc=pc.plot )
  ggsave(paste0( fig.dir, '/irfs_all_v_scale.pdf'), 
         plot=gg.irf.news.v.scale, height=8, width=12)
  
  gg.irf.news.dem.sup.v.scale <- irf_plot_all(df.str.irf.pc.v.scale  %>% 
                     filter(shk.name %in% c('Demand shock', 'Supply shock')), 
                                   data.names, fcast.names, data.labels, pc=pc.plot )
  ggsave(paste0( fig.dir, '/irfs_coef_dem_sup_v_scale.pdf'), 
         plot=gg.irf.news.dem.sup.v.scale, height=10, width=10)
  
  gg.irf.news.dem.v.scale <- irf_plot_all(df.str.irf.pc.v.scale  %>% 
                                                filter(shk.name %in% c('Demand shock')), 
                                              data.names, fcast.names, data.labels, pc=pc.plot )
  ggsave(paste0( fig.dir, '/irfs_coef_dem_v_scale.pdf'), 
         plot=gg.irf.news.dem.v.scale, height=10, width=5)
  ggsave(paste0( fig.dir, '/irfs_coef_dem_v_scale_wide.pdf'), 
         plot=gg.irf.news.dem.v.scale, height=5, width=5)
  
  gg.irf.news.sup.v.scale <- irf_plot_all(df.str.irf.pc.v.scale  %>% 
                                                filter(shk.name %in% c('Supply shock')), 
                                              data.names, fcast.names, data.labels, pc=pc.plot )
  ggsave(paste0( fig.dir, '/irfs_coef_sup_v_scale.pdf'), 
         plot=gg.irf.news.sup.v.scale, height=10, width=5)
  ggsave(paste0( fig.dir, '/irfs_coef_sup_v_scale_wide.pdf'), 
         plot=gg.irf.news.sup.v.scale, height=5, width=5)
  
  gg.irf.news.pol <- irf_plot_all(df.str.irf.pc.v.scale  %>% 
                                        filter( grepl( 'Fiscal|Monetary', shk.name )), 
                                  data.names, fcast.names, data.labels, pc=pc.plot )
  ggsave(paste0( fig.dir, '/irfs_coef_pol_v_scale.pdf'), 
         plot=gg.irf.news.pol, height=10, width=10)
  ggsave(paste0( fig.dir, '/irfs_coef_pol_v_scale_wide.pdf'), 
         plot=gg.irf.news.pol, height=6, width=10)
  
  ### 5.2 The epsilon shocks ###
  df.epsilon.cht <- df.eps.irf.pt %>% 
    select( type, pd, shk.name, data.names, irf.sign ) %>%
    gather(variable, value, data.names) %>%
    mutate( type.cht=case_when(type=='u'~'Surprise',
                               type=='v'~'News',
                               type=='epsilon'~'Combined'),
            variable.cht=data.labels[variable] )
  gg.epsilon.pt <- 
    ggplot( df.epsilon.cht %>% filter( type!='epsilon' ), aes( x=pd-1, y=value*irf.sign ) ) +
      geom_path( data=df.epsilon.cht %>% filter( type=='epsilon' ), lwd=1.1 ) +
      facet_grid( variable.cht~shk.name, scales='free_y') +
      theme_minimal() +
      ylab('Percentage points') + xlab('Quarters') +
      labs(fill='') +
      scale_fill_grey() +
      scale_x_continuous(breaks=number_ticks(6)) +
      theme(legend.position = 'bottom') +
      geom_hline(aes(yintercept=0), lwd=.1)
  ggsave( paste0( fig.dir, '/irf_epsilon.pdf'), gg.epsilon.pt, height=8, width=12)
  gg.epsilon.pt.decomp <- gg.epsilon.pt +
      geom_col( aes(fill=type.cht), alpha=.8 ) +
      geom_path( data=df.epsilon.cht %>% filter( type=='epsilon' ), lwd=1.1 )
  ggsave( paste0( fig.dir, '/irf_epsilon_decomp.pdf'), gg.epsilon.pt.decomp, height=8, width=12 )
  ggsave( paste0( fig.dir, '/irf_epsilon_decomp_tall.pdf'), gg.epsilon.pt.decomp, height=11, width=9 )

  df.epsilon.cht.pc50 <- df.eps.irf.pc %>% filter( pc=='pc.0.5') %>%
    select( type, pd, shk.name, data.names, irf.sign ) %>%
    gather(variable, value, data.names) %>%
    mutate( type.cht=case_when(type=='u'~'Surprise',
                               type=='v'~'News',
                               type=='epsilon'~'Combined'),
            variable.cht=data.labels[variable] )
  gg.epsilon.pc50 <- 
    ggplot( df.epsilon.cht.pc50 %>% filter( type!='epsilon' ), aes( x=pd-1, y=value*irf.sign ) ) +
    geom_path( data=df.epsilon.cht.pc50 %>% filter( type=='epsilon' ), lwd=1.1 ) +
    facet_grid( variable.cht~shk.name, scales='free_y') +
    theme_minimal() +
    ylab('Percentage points') + xlab('Quarters') +
    labs(fill='') +
    scale_fill_grey() +
    scale_x_continuous(breaks=number_ticks(6)) +
    theme(legend.position = 'bottom') +
    geom_hline(aes(yintercept=0), lwd=.1)
  ggsave( paste0( fig.dir, '/irf_epsilon_pc50.pdf'), gg.epsilon.pc50, height=8, width=12)
  gg.epsilon.pt.decomp.pc50 <- gg.epsilon.pc50 +
    geom_col( aes(fill=type.cht), alpha=.8 ) +
    geom_path( data=df.epsilon.cht.pc50 %>% filter( type=='epsilon' ), lwd=1.1 )
  ggsave( paste0( fig.dir, '/irf_epsilon_decomp_pc50.pdf'), gg.epsilon.pt.decomp.pc50, height=8, width=12 )
  ggsave( paste0( fig.dir, '/irf_epsilon_decomp_pc50_tall.pdf'), gg.epsilon.pt.decomp.pc50, height=11, width=9 )
  
  st.pc.plot <- paste0( 'pc.', pc.plot )
  df.epsilon.pc <- df.eps.irf.pc %>% 
    filter( type=='epsilon', pc %in% st.pc.plot ) %>%
    select( type, pd, shk.name, pc, data.names, irf.sign ) %>%
    gather(variable, value, data.names) %>%
    mutate( value=irf.sign*value) %>%
    spread( pc, value) %>%
    mutate( type.cht=case_when(type=='u'~'Surprise',
                               type=='v'~'News',
                               type=='epsilon'~'Combined'),
            variable.cht=data.labels[variable] ) %>%
    left_join(df.epsilon.cht %>% filter(type=='epsilon'))
  gg.epsilon.pc <- gg.epsilon.pt +
    geom_path(data=df.epsilon.pc, aes_string(y=st.pc.plot[2]), lty=3, lwd=.5 ) +
    geom_path(data=df.epsilon.pc, aes_string(y=st.pc.plot[3]), lty=3, lwd=.5 ) +
    geom_path(data=df.epsilon.pc, aes_string(y=st.pc.plot[1]), lty=2, lwd=.6 ) +
    geom_path(data=df.epsilon.pc, aes_string(y=st.pc.plot[4]), lty=2, lwd=.6 ) 
  ggsave(paste0( fig.dir, '/irf_epsilon_pc.pdf'), gg.epsilon.pc, height=6, width=9)
  ggsave(paste0( fig.dir, '/irf_epsilon_pc_tall.pdf'), gg.epsilon.pc, height=11, width=9)
  
  ### 5.3 The cumulative fiscal multipliers ###
  if(do.cumul){
    for( this.shk in c( 'v', 'u', 'epsilon' )){
      ggplot( data=l.df.cumul.pc[[this.shk]] %>% filter(pc==.5), aes(x=pd, y=y.cum)) +
        geom_path( lwd=2 ) +
        geom_path( data=l.df.cumul.pc[[this.shk]] %>% filter( pc %in% pc.plot ),
                   aes(lty=as.factor(pc) ), lwd=.5 ) +
        scale_linetype_manual(values=c(2, 3, 3, 2) ) +
        geom_point( data=l.df.cumul.mults[[this.shk]], 
                    aes( group=label, pch=label, color=label ), size=5 ) +
        guides(linetype='none') +
        labs( pch=NULL, color=NULL ) +
        ylab('Percentage points') +
        xlab('Horizon, h') +
        theme_minimal() +
        theme(legend.position = 'bottom') +
        geom_hline(aes(yintercept=0), lwd=.5)
      ggsave(paste0( fig.dir, '/mult_verify_', this.shk, '.pdf'), height=6, width=9)
      ggsave(paste0( fig.dir, '/mult_verify_', this.shk, '_wide.pdf'), height=4, width=10)
      
      
      ggplot( data=l.df.cumul.all.pc[[this.shk]] %>% filter(pc==.5), aes(x=pd, y=y.cum)) +
        geom_path( lwd=2 ) +
        geom_path( data=l.df.cumul.all.pc[[this.shk]] %>% filter( pc %in% pc.plot ),
                   aes(lty=as.factor(pc) ), lwd=.5 ) +
        scale_linetype_manual(values=c(2, 3, 3, 2) ) +
        geom_point( data=l.df.cumul.mults.all[[this.shk]], 
                    aes( group=label, pch=label, color=label ), size=5 ) +
        guides(linetype='none') +
        facet_wrap(~shk.name) +
        labs( pch=NULL, color=NULL ) +
        ylab('Percentage points') +
        xlab('Horizon, h') +
        theme_minimal() +
        theme(legend.position = 'bottom') +
        geom_hline(aes(yintercept=0), lwd=.5)
      ggsave(paste0( fig.dir, '/mult_verify_all_', this.shk, '.pdf'), height=6, width=9)
      ggsave(paste0( fig.dir, '/mult_verify_all_', this.shk, '_wide.pdf'), height=4, width=9)
    }
  }
  
  ### 5.4 The counterfactual policies ###
  if( do.cf ){
    
    df.baseline.cht <- df.str.irf.pc %>% 
      filter(pc=='pc.0.5') %>%
      select( pd, type, dimension, pc, data.names, shk.name, irf.sign ) %>%
      gather(variable, pc.0.5, -pd, -type, -dimension, -pc, -shk.name, -irf.sign ) %>%
      mutate(type.cht=case_when(type=='u'~' Surprise',
                                type=='v'~'News',
                                type=='epsilon'~'Combined'),
             variable.cht=data.labels[variable])
    
    if( !is.null(l.cf$consol) ){
      df.baseline.consol <- df.baseline.cht %>% 
        filter(shk.name != ' Fiscal stimulus shock') %>%
        filter(shk.name %in% c('Supply shock', 'Demand shock'))
      for( this.target in names(l.cf$consol) ){
        gg.irf.cf.pc <- irf_plot_all(l.cf$consol[[this.target]]$df %>% 
                                       filter(shk.name %in% c('Supply shock', 'Demand shock')), data.names, fcast.names, 
                                       data.labels, pc=pc.plot, do.pt = TRUE ) + 
          # geom_point(data=df.baseline.consol, aes(x=pd, y=pc.0.5), col='blue', size=1.1, cex=13 )
          geom_line(data=df.baseline.consol, 
                    aes(x=pd-1, y=pc.0.5), col='blue', lwd=.9 )
        gg.irf.cf.pc$layers <- c( gg.irf.cf.pc$layers[[6]], gg.irf.cf.pc$layers[1:5] )
          # Change the order of the lines
        ggsave( paste0(fig.dir, '/cf_consol_', this.target, '.pdf'), gg.irf.cf.pc,
                height=8, width=8 )
        ggsave( paste0(fig.dir, '/cf_consol_', this.target, '_wide.pdf'), gg.irf.cf.pc,
                height=6, width=9 )
      }
    }
    
    if( !is.null(l.cf$mon) ){
      df.baseline.mon <- df.baseline.cht %>% 
        filter(shk.name != ' Monetary policy shock') %>%
        filter(shk.name %in% c('Supply shock', 'Demand shock') )
      for( this.target in names(l.cf$mon) ){
        gg.irf.cf.pc <- irf_plot_all(l.cf$mon[[this.target]]$df %>%
                                       filter(shk.name %in% c('Supply shock', 'Demand shock') ), data.names, fcast.names, 
                                     data.labels, pc=pc.plot, do.pt = TRUE ) + 
          geom_line(data=df.baseline.mon, aes(x=pd-1, y=pc.0.5), col='blue', lwd=.9 )
        gg.irf.cf.pc$layers <- c( gg.irf.cf.pc$layers[[6]], gg.irf.cf.pc$layers[1:5] )
          # Change the order of the lines
        ggsave( paste0(fig.dir, '/cf_mon_', this.target, '.pdf'), gg.irf.cf.pc,
                height=8, width=8 )
      }
    }
    
    if( !is.null(l.cf$both) ){
      df.baseline.both <- df.baseline.cht  %>%
        filter(shk.name %in% c('Supply shock', 'Demand shock'))
      for( this.target in names(l.cf$both) ){
        gg.irf.cf.pc <- irf_plot_all(l.cf$both[[this.target]]$df %>%
                                       filter(shk.name %in% c('Supply shock', 'Demand shock')), data.names, fcast.names, 
                                     data.labels, pc=pc.plot, do.pt = TRUE ) + 
          geom_line(data=df.baseline.both, aes(x=pd-1, y=pc.0.5), col='blue', lwd=.9 )
        gg.irf.cf.pc$layers <- c( gg.irf.cf.pc$layers[[6]], gg.irf.cf.pc$layers[1:5] )
        # Change the order of the lines
        ggsave( paste0(fig.dir, '/cf_both_', this.target, '.pdf'), gg.irf.cf.pc,
                height=8, width=8 )
      }
    }
    
    ## Make the inflation and dual mandate chart: Monetary policy ##
    df.cht.mon.cf <- full_join( df.str.irf.pc %>% mutate( label='Baseline' ),
                                l.cf$mon$infl$df %>% mutate( label='Inflation targeting' ) ) %>%
      full_join( l.cf$mon$y$df %>% mutate( label='Output stabilization' ) ) %>%
      full_join( l.cf$mon$infl_y$df %>% mutate( label='Dual mandate' ) ) %>%
      filter( shk.name %in% c('Supply shock', 'Demand shock'), pc %in% c('pc.0.5', 'pt') ) %>%
      ungroup() %>%
      select(-irf.sign, -shk.name) %>% 
      left_join( df.labels )
    irf_plot_compare_pt_multi( df.cht.mon.cf, 'label', data.names=c('infl', 'y', 'ir03'), 
                               data.labels = data.labels )
    ggsave(filename = paste0( fig.dir, '/mon_cf_compare.pdf'), height=8, width=8 )
    ggsave(filename = paste0( fig.dir, '/mon_cf_compare_wide.pdf'), height=5, width=8 )
    
    ## Make the inflation and dual mandate chart: Fiscal policy ##
    df.cht.fisc.cf <- full_join( df.str.irf.pc %>% mutate( label='Baseline' ),
                                l.cf$consol$infl$df %>% mutate( label='Inflation targeting' ) ) %>%
      full_join( l.cf$consol$y$df %>% mutate( label='Output stabilization' ) ) %>%
      full_join( l.cf$consol$infl_y$df %>% mutate( label='Dual mandate' ) ) %>%
      filter( shk.name %in% c('Supply shock', 'Demand shock'), pc %in% c('pc.0.5', 'pt') ) %>%
      ungroup() %>%
      select(-irf.sign, -shk.name) %>% 
      left_join( df.labels )
    irf_plot_compare_pt_multi( df.cht.fisc.cf, 'label', data.names=c('infl', 'y', 'gov', 'tax'), data.labels )
    ggsave(filename = paste0( fig.dir, '/consol_cf_compare.pdf'), height=8, width=8 )
    ggsave(filename = paste0( fig.dir, '/consol_cf_compare_wide.pdf'), height=5, width=8 )
    
    ## Make the passive policies chart
    df.cht.passive.cf <- full_join( df.str.irf.pc %>% mutate( label='Baseline' ),
                                    l.cf$mon$ir03$df %>% mutate( label='Interest rate peg' ) ) %>%
      full_join( l.cf$consol$gov$df %>% mutate( label='Fixed government spending' ) ) %>%
      full_join( l.cf$consol$tax$df %>% mutate( label='Fixed taxes' ) ) %>%
      filter( shk.name %in% c('Supply shock', 'Demand shock'), pc %in% c('pc.0.5', 'pt') ) %>%
      ungroup() %>%
      select(-irf.sign, -shk.name) %>% 
      left_join( df.labels )
    irf_plot_compare_pt_multi( df.cht.passive.cf, 'label', data.names=c('infl', 'y', 'ir03', 'gov', 'tax'), data.labels )
    ggsave(filename = paste0( fig.dir, '/passive_cf_compare.pdf'), height=10, width=8 )
    ggsave(filename = paste0( fig.dir, '/passive_cf_compare_wide.pdf'), height=5, width=8 )
    
    ## Make the variance reduction table
    v.target <- c( 'Inflation'='infl', 'Output'='y', 'Dual mandate'='infl_y' )
    df.irf.norm.ratio <- 
      lapply( v.target, function(y) sapply( l.cf, function(x) x[[y]]$irf.norm$ratio ) ) %>%
      reduce( cbind ) %>%
      as.data.frame
    df.irf.norm.ratio <- 
      cbind( Variable=data.labels[rownames(df.irf.norm.ratio)] %>%
               str_replace('Inf', '   Inf') %>% str_replace('Out', '  Out') , 
             df.irf.norm.ratio )
    df.irf.norm.ratio <- df.irf.norm.ratio[order(df.irf.norm.ratio$Variable),]
    
    col.rename.fn <- function(x){
      x %>% str_replace( 'consol.*', 'Fisc.' ) %>% 
        str_replace( 'mon.*', 'Mon.' ) %>% 
        str_replace( 'both.*', 'Joint' ) %>%
        str_replace( 'Variable', 'Policy used' ) %>%
        paste0( '\\textit{', ., '}')
    }
    print( xtable(df.irf.norm.ratio ),
           include.rownames=FALSE, sanitize.colnames.function=col.rename.fn,
           add.to.row=list(pos=list(-1),  
              command=paste0("\\toprule Target: & \\multicolumn{3}{c}{Inflation} &  \\multicolumn{3}{c}{Output} &  \\multicolumn{3}{c}{Dual Mandate} \\\\ \n",
                       " \\cmidrule(l){2-4} \\cmidrule(l){5-7} \\cmidrule(l){8-10} \n") ),
           floating=FALSE, hline.after=c(0,nrow(df.irf.norm.ratio) ), 
           booktabs=TRUE,
           file=paste0( fig.dir, '/cf_irf_norm.tex') )
    
    
  }
  
  
  ### 5.6 Plot the data ####
  ggplot( df.clean %>% gather(variable, value, -date ) %>%
            separate(variable, c('variable','fcast'), sep='_' ) %>%
            mutate(fcast=ifelse(fcast=='f', 'ML Forecast', ' Data') %>% replace_na(' Data')) %>%
            full_join( df %>% gather(variable, value, -date ) %>%
                         separate(variable, c('variable','fcast'), sep='_' ) %>% 
                         filter(!is.na(fcast)) %>%
                         mutate(fcast=ifelse(fcast=='f', 'Original Forecast', 'Data') )
            ), 
          aes(x=date, y=value, group=fcast, color=fcast) ) +
    geom_line(lwd=.9) +
    facet_wrap(~variable, scales='free_y') + 
    theme_minimal() +
    theme(legend.position = 'bottom') +
    labs(color='')
  ggsave( paste0(fig.dir, '/data_original.pdf'), height=8, width=12)
  
  ggplot( df.clean %>% gather(variable, value, -date ) %>%
            separate(variable, c('variable','fcast'), sep='_' ) %>%
            mutate(fcast=ifelse(fcast=='f', 'ML Forecast', ' Data') %>% replace_na('Data')), 
          aes(x=date, y=value, group=fcast, color=fcast) ) +
    geom_line(lwd=.9) +
    facet_wrap(~variable, scales='free_y') + 
    theme_minimal() +
    theme(legend.position = 'bottom') +
    labs(color='')
  ggsave( paste0(fig.dir, '/data_lasso.pdf'), height=8, width=12)
  
  ggplot( l.rf.ols$df %>%
            select(date, data.names, fcast.names) %>% 
            gather(variable, value, -date ) %>%
            separate(variable, c('variable','fcast'), sep='_' ) %>%
            mutate(fcast=ifelse(fcast=='f', 'Cleaned Forecast', ' Data')  %>% replace_na(' Data'),
                   variable.cht=data.labels[variable] ), 
          aes(x=date, y=value, group=fcast, color=fcast, lty=fcast) ) +
    geom_line(lwd=.9) +
    xlab('') + ylab('') +
    facet_wrap(~variable.cht, scales='free_y') + 
    theme_minimal() +
    theme(legend.position = 'bottom') +
    labs(lty='', color='')
  ggsave( paste0(fig.dir, '/data_cleaned.pdf'), height=8, width=12)
  ggsave( paste0(fig.dir, '/data_cleaned_wide.pdf'), height=6, width=12)
  
  ### 5.7 Plot the estimated shocks ####
  df.struct.resid <- l.struct$pt$df.resid %>%
    gather( variable, value, -date ) %>%
    mutate( type.cht=case_when(grepl('uhat', variable)~'Surprise',
                               grepl('vhat', variable)~'News',
                               grepl('epshat', variable)~'Combined'),
            dimension=str_extract( variable, '([0-9])*$') %>% as.numeric ) %>%
    left_join( df.labels )
  
  ggplot( df.struct.resid %>% filter(type.cht!='Combined'), 
          aes(x=date, y=value*irf.sign, color=type.cht, shape=type.cht, linetype=type.cht) ) +
    geom_line() +
    facet_wrap( ~shk.name ) +
    theme_minimal() + 
    theme(legend.position = 'bottom') +
    ylab('') + xlab('') +
    labs( color='Shock type', linetype='Shock type')
  ggsave(paste0( fig.dir, '/shocks.pdf'), 
         height=8, width=12)
  ggsave(paste0( fig.dir, '/shocks_tall.pdf'), 
         height=11, width=9)
  
  
  #### 6. Make Tables ####
  
  ## 6.1 Variance decomposition ##
  message('Making tables...')
  var.table <- df.var.decomp %>%
    filter( pd == max(pd) ) %>%
    left_join( {.} %>% group_by(variable) %>% summarise( total=sum(variance[type=='epsilon']) ) ) %>%
    mutate( variance.pc = 100 * variance / total ) %>%
    left_join( df.labels %>% select(dimension, shk.name ) %>% mutate( shk.name=str_remove(shk.name, ' shock' ) ) ) %>%
    full_join( {.} %>% group_by(variable, type ) %>% 
                 summarise( variance.pc=sum(variance.pc) ) %>% 
                 mutate( shk.name='Total') ) %>%
    mutate( type.name=case_when(type=='contemp'~' Surprise',
                                type=='news'~' News',
                                type=='epsilon'~'Total'),
            variable.name=data.labels[variable] %>% str_remove('\\n') %>%
              str_replace( 'Government', 'Gov.'),
            shk.name = shk.name %>%
              str_replace( 'Government', 'Gov.') %>% 
              str_replace( 'Monetary', 'Mon.') ) %>%
    select( -dimension, -variance, -total, -pd, -variable, -type ) %>%
    spread( shk.name, variance.pc ) %>%
    select( Variable=variable.name, Type=type.name, everything(), -Total, Total ) %>%
    arrange( Variable, Type) %>%
    mutate( Variable = c( Variable[1], ifelse(Variable[-1]==Variable[-n()], '', Variable[-1] ) ) )
      # The point estimate table
  
  var.table.pc50 <- df.var.decomp.boot.pc %>% filter(pc=='pc.0.5') %>%
    ungroup() %>%
    filter( pd == max(pd) ) %>%
    mutate( type.name=case_when(type=='contemp'~' Surprise',
                                type=='news'~' News',
                                type=='epsilon'~'Total'),
            variable.name=data.labels[variable] %>% str_remove('\\n') %>%
              str_replace( 'Government', 'Gov.'),
            shk.name = shk.name %>%
              str_replace( 'Government', 'Gov.') %>% 
              str_replace( 'Monetary', 'Mon.') ) %>%
    select( -dimension, -variance, -total, -pd, -variable, -type, -pc, -i.boot ) %>%
    spread( shk.name, variance.pc ) %>%
    select( Variable=variable.name, Type=type.name, everything(), -Total, Total ) %>%
    arrange( Variable, Type) %>%
    rbind( {.} %>% group_by(Type) %>% summarise_if(is.numeric, mean) %>%
             mutate( Variable='Unweighted average') %>% select(Variable, everything() ) ) %>%
    mutate( Variable = c( Variable[1], ifelse(Variable[-1]==Variable[-n()], '', Variable[-1] ) ) )
    # The median table
  
  st.var.pc <- paste0( 'pc.', var.pc )
  var.table.pc <- df.var.decomp.boot.pc %>% ungroup() %>% 
                 filter( pc %in% c( 'pc.0.5', st.var.pc), pd==max(pd) ) %>% 
                 select(-variance, -i.boot, -total) %>% 
                 spread( pc, variance.pc) %>% 
                 mutate( zz_range = paste0( '(', formatC(get(st.var.pc[1]), digits=1, format='f' ), 
                                            ', ', formatC(get(st.var.pc[2]), digits=1, format='f' ), ')' ) ) %>%
                 select( type, dimension, shk.name, variable, pc.0.5, zz_range ) %>%
    mutate( type.name=case_when(type=='contemp'~' Surprise',
                                type=='news'~' News',
                                type=='epsilon'~'Total'),
            variable.name=data.labels[variable] %>% str_remove('\\n') %>%
              str_replace( 'Government', 'Gov.'),
            shk.name = shk.name %>%
              str_replace( 'Government', 'Gov.') %>% 
              str_replace( 'Monetary', 'Mon.'),
            pc.0.5=formatC(pc.0.5, digits=1, format='f')) %>%
    select( -dimension, -variable, -type ) %>%
    gather( measure, value, -shk.name, -type.name, -variable.name ) %>%
    spread( shk.name, value ) %>%
    select( Variable=variable.name, Type=type.name, everything(), -Total, Total, -measure ) %>%
    arrange( Variable, Type) %>%
    mutate( Variable = c( Variable[1], ifelse(Variable[-1]==Variable[-n()], '', Variable[-1] ) ),
            Type = c( Type[1], ifelse(Type[-1]==Type[-n()], '', Type[-1] ) ) )
  # The CI table
  
  print( xtable(var.table, label = paste0('t:var_decomp_', this.vars ), digits = 1,
                caption=paste0( 'Forecast error variance decomposition, ', 
                                n.pds.irf, '-quarter-ahead horizon' ) ), 
         include.rownames=FALSE, hline.after=c(-1,-1,0,nrow(var.table) ), 
         # floating.environment = 'sidewaystable',
         file=paste0( fig.dir, '/var_decomp.tex'), floating=FALSE )
  
  
  print( xtable(var.table.pc50, label = 't:var_decomp_pc50', digits = 1,
                caption=paste0( 'Forecast error variance decomposition, ', 
                                n.pds.irf, '-quarter-ahead horizon' ) ), 
         include.rownames=FALSE, hline.after=c(-1,-1,0,nrow(var.table), nrow(var.table.pc50) ), 
         # floating.environment = 'sidewaystable',
         file=paste0( fig.dir, '/var_decomp_pc50.tex'), floating=FALSE )
  
  print( xtable(var.table.pc, label = paste0('t:var_decomp_pc_', this.vars ), digits = 1,
                caption=paste0( 'Forecast error variance decomposition, ', 
                                n.pds.irf, '-quarter-ahead horizon' ) ), 
         include.rownames=FALSE, hline.after=c(-1,-1,0,nrow(var.table.pc) ), 
         # floating.environment = 'sidewaystable',
         file=paste0( fig.dir, '/var_decomp_pc.tex'), floating=FALSE )
  message('...complete')
  # Save to Latex
  
  ## 6.2 The variance shares D ##
  l.D.u.2 <- lapply( 1:n.boot, function(i)
      data.frame( boot=i, dimension =1:l.rf.ols$n.var, D.u.2=l.struct$boot[[i]]$D.u.2 %>% diag )
    )
  df.D.2 <- l.D.u.2 %>%
    rbindlist() %>% mutate(D.v.2=1-D.u.2)
  df.D.2.smy <- lapply( c( .01, .05, .1, .5, .9, .95, .99), function(x) df.D.2 %>%
                          group_by(dimension) %>% 
                          summarise(D.u.2=quantile(D.u.2,x), D.v.2=quantile(D.v.2,x)) %>%
                          mutate(pc=x)) %>%
    rbindlist()
  tab.D.2.smy <- df.D.2.smy %>% 
    filter(pc %in% c(.1, .5, .9) ) %>%
    left_join(df.labels %>% select(-irf.sign), by='dimension') %>%
    gather( measure, value, -dimension, -pc, -shk.name ) %>%
    spread( pc, value ) %>%
    mutate( ` median`=sprintf("%.2f", `0.5`),
            ci.10.90=paste0( '(', sprintf("%.2f", `0.1`), ', ', sprintf("%.2f", `0.9`), ')' ) ) %>%
              select( shk.name, measure, ` median`, ci.10.90 ) %>%
    gather( variable, value, -shk.name, -measure ) %>%
    spread(shk.name, value) %>%
    select(-variable) %>%
    mutate( measure=case_when(measure=='D.u.2'~'$D_u^2$',
                              measure=='D.v.2'~'$D_v^2$' ),
            measure=c(measure[1], ifelse( measure[-n()]==measure[-1], '', measure[-1] ) )) %>%
    rename(' '=measure)
    
  xtable( tab.D.2.smy, align = c('l','l', rep('c', l.rf.ols$n.var)) ) %>%
    print(include.rownames=FALSE, sanitize.text.function=identity, 
          floating=FALSE, hline.after=c(-1,-1,0,nrow(tab.D.2.smy) ),
          file=paste0( fig.dir, '/D_shares.tex'))
  
  #### 7. Save the output ####
  if(do.cf){
    save( l.rf.ols, l.rf.adj.alt, l.struct, 
          df.eps.irf.pt, df.eps.irf.pc, df.var.decomp.boot.frac, 
          df.str.all, df.str.irf.pc, df.epsilon.cht, df.labels, df.clean, 
          df.str.irf.pt, df.str.irf.boot,
          l.df.cumul,l.df.cumul.pc,  l.df.cumul.mults, 
          df.str.irf.pc.v.scale, l.cf,
          file=paste0(fig.dir, '/est.rdata'))
  }else{
    save( l.rf.ols, l.rf.adj.alt, l.struct, 
          df.eps.irf.pt, df.eps.irf.pc, df.var.decomp.boot.frac, 
          df.str.all, df.str.irf.pc, df.epsilon.cht, df.labels, df.clean, 
          df.str.irf.pt, df.str.irf.boot,
          l.df.cumul,l.df.cumul.pc,  l.df.cumul.mults, 
          df.str.irf.pc.v.scale,
          file=paste0(fig.dir, '/est.rdata'))
  }

}

