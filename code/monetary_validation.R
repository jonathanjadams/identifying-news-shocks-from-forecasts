##########################################################################################
#
# monetary_validation.R
#
# Code to assemble the shocks for the monetary policy validation
#
# Philip Barrett, Washington DC
# First version: 29jun2022
#
##########################################################################################

#### 0. Set up ####
rm(list=ls())
library(tidyverse)
library(broom)
library(lubridate)
library(vars)
library(zoo)
library(readstata13)
library(lpirfs)
setwd('C:\\Users\\pbarrett\\OneDrive - International Monetary Fund (PRD)\\identification from RE/')
source("code/extract_varirf.R")
select <- dplyr::select

#### 1. Controls ####
date.start <- '1969-01-01'
# baseline.folder <- 'figures/28jul23/baseline_ml_1_lag/'
baseline.folder <- 'figures/final_replication/baseline_ml_1_lag/'

data.labels <- c( ir03='3-month \ninterest rate', 
                  ir5='5-year \ninterest rate', 
                  un='Unemployment', infl='Inflation', 
                  y=' Output', ir5='Five-year interest rate', 
                  sur='Surplus (percent of GDP)', tax=' Taxes',
                  hs='Housing starts', indpro='Industrial production',
                  fed='  Federal government \nconsumption',
                  gov='  Government \nconsumption' )
var.order <- c( 'mon.shk', 'ir03', 'hs', 'infl', 'gov', 'tax', 'y' )
st.shk.source <- 'data/monetary_compare/'
l.sources <- list( 
  romer_wieland=list( file='romer_romer_2004_updated/RR_monetary_shock_quarterly.dta',
                      shk_names=c( # resid='Romer-Romer 2004 (original)',
                                   # resid_romer='Romer-Romer 2004 (extended to 2007)',
                                   resid_full='Romer-Romer 2004 (2007 data)'),
                      read_with=read.dta13, freq='q' ),
  bauer_swanson=list( file='bauer_swanson/FOMC_Bauer_Swanson.csv',
                      shk_names=c( MPS='Bauer-Swanson HFI',
                                   MPS_ORTH='Bauer-Swanson HFI (orthogonalized)' ),
                      read_with=read_csv, freq='m' ),
  jarocinski_karadi=list( file='jarocinski_karadi/us_shocks.csv',
                      shk_names=c( ff4_hf='Jarocinski-Karadi, HFI from Fed Funds' #,
                                   ),
                      read_with=read_csv, freq='m' ),
  romer_romer_narrative=list( file='romer_romer/romer_romer_narrative.csv',
                          shk_names=c( RRDUMMY='Romer-Romer 2023' ),
                          read_with=read_csv, freq='q' ) )
n.lags <- 4 # 1
plot.pc <- c(.1)
st.plot.pc <- paste0( 'pc.', sapply( plot.pc, function(x) c(x,1-x) ) %>% c )
n.boot <- 10 ; n.pds <- 25 ; rescale.cum.pds <- 8
news.shk.type <- 'epsilon'  # Our shock type.  Can also be 'u' or 'v'

#### 2. Read in and consistently format the data ####
l.df.in <- NULL ; n.df.in <- l.sources %>% length()
for( i in 1:n.df.in ){
# Loop over the sources
  
  this.name <- l.sources %>% names() %>% .[i]
    # Short name of the data source
  this.df.in <- paste0( st.shk.source, l.sources[[i]]$file ) %>% 
    l.sources[[i]]$read_with(file = .)
      # Basic read in
  this.df <- this.df.in %>%
    { if( this.name == 'romer_romer_narrative' ) {.} %>% 
        separate(1, c('year', 'quarter'), sep=':' ) %>%
        mutate(month=3*as.numeric(quarter)) else . } %>%
    { if( all( c('year', 'month') %in% names(.) ) ) {.} %>% 
        mutate( date=make_date(year, month, 1 ) ) else . } %>%
      # Fix the date format
    { if( 'Date' %in% names(.) ) {.} %>% rename(date=Date) else . } %>%
    select( date, l.sources[[i]]$shk_names %>% names ) %>%
      # Isolate the shock variables
    mutate( date=(date %>% as.yearqtr %>% as.Date) %m+% months(2) ) %>%
    complete( date=seq.Date(min(date), max(date), by="quarter" ) ) %>%
    group_by(date) %>%
    summarise_all(sum) %>%
      # Compute quarterly averages
    gather( shk.short, mon.shk, -date ) %>%
    mutate( mon.shk = replace_na(mon.shk,0)) %>%
    mutate( shk.long = l.sources[[i]]$shk_names[shk.short],
            source = this.name ) %>%
      # Convert to long format
    select( source, date, shk.short, shk.long, mon.shk )
  
  l.df.in[[ this.name ]] <- this.df %>% ungroup()
}
df.shk <- l.df.in %>% reduce(rbind) %>% filter( date >= date.start )
  # The data frame of shocks

## 2.2 A summary table of the shocks ##
df.sumy <- df.shk %>%
  group_by(source, shk.short, shk.long) %>%
  summarise( Start=min(date), End=max(date), N=n() ) %>%
  mutate( `Orig. Freq.`=l.sources[[source]]$freq %>% toupper ) %>% 
  ungroup() %>%
  select( Shock=shk.long, `Orig. Freq.`, Start, End, N )

## 2.3 Merge in the regression variables ##
keep.vars <- c('df.clean', 'df.eps.irf.pc')
st.baseline.names <- load( paste0(baseline.folder, 'est.rdata') )
rm( list=st.baseline.names[!(st.baseline.names %in% keep.vars )] )
  # Load only the variables that we want to keep
df.all <- left_join( df.clean, df.shk, by='date' ) %>%
  arrange( shk.short, date) %>%
  select( source, shk.short, shk.long, date, mon.shk, everything() )


#### 3. Compare the sizes of the shocks

## 3.1 Raw data
ggplot( df.shk, aes(x=date, y=mon.shk, group=shk.long, color=shk.long, shape=shk.long ) ) +
  geom_line() +
  geom_point() +
  theme_minimal() + 
  xlab('') + ylab('') +
  theme(legend.position = 'bottom') +
  guides(color=guide_legend(nrow=3,byrow=TRUE)) +
  theme(legend.title=element_blank())


## 3.2 Project onto interest rates
df.lm <- 
  df.all %>% filter(!is.na(shk.short)) %>%
  nest_by(shk.short, shk.long) %>%
  mutate( mod=list(lm(ir03 ~ mon.shk, data=data ) )) 
    # Simple projection onto contemporaneous interest rates
df.lm.coefs <- df.lm %>% summarize( tidy( mod ) ) 
df.lm.sumy <- df.lm %>% summarize( glance( mod ) ) 
df.lm.fitted <- cbind( df.all %>% filter(!is.na(shk.short)) %>% arrange( shk.short ) %>% select(date),
                       df.lm %>% summarize( augment( mod ) ) )
  # The coefficients and the fitted values
ggplot( df.lm.fitted, aes(x=date, y=.fitted ) ) +
  geom_line(aes( color=shk.long, group=shk.long)) +
  geom_point(aes(color=shk.long, group=shk.long, shape=shk.long)) +
  geom_line( data=df.clean, aes(x=date, y=ir03), lwd=1.1, col='black' ) +
  theme_minimal() + 
  xlab('') + ylab('') +
  theme(legend.position = 'bottom') +
  guides(color=guide_legend(nrow=3,byrow=TRUE)) +
  theme(legend.title=element_blank())


#### 4. Run the dynamic estimates for each shock ####

## 4.1 Set up ##
df.names <- df.shk %>%
  group_by( source, shk.short, shk.long ) %>%
  summarise() %>% ungroup()
    # Set up the data frame of shock names
n.shks <- df.names %>% nrow() ; 
l.df.irf <- l.df.irf.var <- l.irf.lp <- l.df.irf.lp <- NULL

## 4.2 Main loop ##
for( i in 1:n.shks ){
  
  message('### Shock number ', i )
  
  this.shk.short <- df.names$shk.short[i]
  this.shk.long <- df.names$shk.long[i]
    # Shock names, short and long
  df.this.shk <- df.all %>% filter( shk.short == this.shk.short ) %>%
    select( var.order )
    # The data
  
  ## By VAR ##
  message('Calculating VAR')
  mon.var <- VAR(df.this.shk, type='none', p = n.lags )
  Bmat <- diag( ncol(df.this.shk) ) ; Bmat[lower.tri(Bmat, diag = TRUE)] <- NA
  mon.svar <- SVAR(mon.var, Bmat = Bmat )
  irf.mon <- irf(mon.svar, impulse='mon.shk', n.ahead = 24, 
                 cumulative=FALSE, ci=.9, runs=n.boot )
    # The structural IRF
  l.df.irf.var[[i]] <- extract_varirf(irf.mon) %>%
    rbind( {.} %>% head(1) %>% mutate_all( function(x) 0*x ) %>% mutate( period=-1), 
           . ) %>%
    select( !contains(c('upper', 'lower')) ) %>%
    gather( variable, irf.var, -period ) %>%
    separate( variable, into=c('temp1', 'temp2', 'variable'), sep='_' ) %>%
    select( -temp1, -temp2 ) %>%
    mutate( method=this.shk.long, period=period+1 ) %>%
    select( method, pd=period, everything() )
  
  ## By Local projection ##
  message('Calculating LP')
  l.irf.lp[[i]] <- lp_lin_iv( df.this.shk %>% select( -mon.shk ), #shock_type=1,
               lags_endog_lin=n.lags, num_cores=2,
               hor=n.pds, trend=0, confint=1.96,
               shock=df.this.shk %>% select(mon.shk) )
  l.df.irf.lp[[i]] <- cbind(0, l.irf.lp[[i]]$irf_lin_mean ) %>% 
    t %>% as.data.frame() %>%
    set_names( var.order[-1] ) %>%
    mutate( pd=0:n.pds, method=this.shk.long ) %>%
    select( method, pd, everything() ) %>%
    gather( variable, irf.lp, -method, -pd )
  
  l.df.irf[[i]] <- full_join( l.df.irf.var[[i]], l.df.irf.lp[[i]] )
  
}
df.irf <- l.df.irf %>% reduce(rbind) %>%
  mutate( variable.name=data.labels[variable] )
save(df.irf, date.start, baseline.folder, data.labels, 
     var.order, st.shk.source, l.sources, n.lags, plot.pc,
     st.plot.pc, n.boot, n.pds, rescale.cum.pds, rescale.cum.pds, 
     news.shk.type, df.all, 
     file=paste0(baseline.folder,'monetary_validation_',
                 n.lags, 'lags.rdata') )

## 4.3 Rescale ##
df.irf.rescale <- df.irf %>%
  left_join( df.irf %>% filter(variable=='ir03') %>%
               group_by(method) %>%
               summarise( ir.var.init=irf.var[pd==1], ir.lp.init=irf.lp[pd==1],
                          ir.var.cum=irf.var[pd %in% 1:rescale.cum.pds] %>% sum, 
                          ir.lp.cum=irf.lp[pd %in% 1:rescale.cum.pds] %>% sum ) ) %>%
  mutate( irf.var.init=irf.var / ir.var.init, irf.lp.init=irf.lp / ir.lp.init,
          irf.var.cum=irf.var / ir.var.cum, irf.lp.cum=irf.lp / ir.lp.cum )


#### 5. Our baseline estimates ####
df.baseline <- df.eps.irf.pc %>%
  filter( type==news.shk.type, shk.name==' Monetary policy shock', 
          pc %in% c( 'pc.0.5', st.plot.pc ) ) %>% 
  ungroup() %>%
  select( pd, pc, var.order[-1] ) %>%
  gather( variable, baseline, -pd, -pc ) %>%
  cbind( {.} %>% filter(variable=='ir03', pc=='pc.0.5') %>%
               summarise( ir.init=baseline[pd==1],
                          ir.cum=baseline[pd %in% 1:rescale.cum.pds] %>% sum ) ) %>%
  mutate( baseline.init=baseline / ir.init, baseline.cum=baseline / ir.cum ) %>%
  gather( case, value, -pd, -pc, -variable ) %>%
  unite( temp, case, pc, sep = '.' ) %>%
  spread( temp, value ) %>%
  mutate( variable.name=data.labels[variable] )


#### 6. Make charts ####

## 6.0 Code to highlight the median
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
gg_shapes <- function(n) {
  shapes <- c( 15:18, 3, 6, 12, 4 )
  shapes[1:n]
}

# Dynamically generate default color values, but have Baseline="black".
long.names <- df.irf.rescale$method %>% na.omit() %>% unique()
col.values <- gg_color_hue(length(long.names))
names(col.values) <- long.names
col.values <- c(col.values, c(Baseline="black"))
  # Set the color values with baseline in black
shape.values <- gg_shapes(length(long.names))
names(shape.values) <- long.names
shape.values <- c(shape.values, c(Baseline=NA))
  # Set the shape values

## 6.1 Naive charts ##
ggplot(df.irf %>% filter(!is.na(variable.name)), aes(x=pd ) ) +
  geom_hline(aes(yintercept=0)) +
  geom_path(aes(y=irf.var, group=method, color=method, shape=method)) + 
  geom_point(aes(y=irf.var, group=method, color=method, shape=method)) +
  theme_minimal() + 
  xlab('') + ylab('') +
  facet_wrap(~variable.name, scales='free_y') +
  theme( legend.position = 'bottom' )
ggsave(paste0(baseline.folder, 'monetary_compare_var_naive_', n.lags, '_lags.pdf' ),
       height=8, width=8 )
ggsave(paste0(baseline.folder, 'monetary_compare_var_naive_', n.lags, '_lags_wide.pdf' ),
       height=4, width=8 )

ggplot(df.irf %>% filter(!is.na(variable.name)), aes(x=pd) ) +
  geom_hline(aes(yintercept=0)) +
  geom_path(aes(y=irf.lp, group=method, color=method, shape=method)) + 
  geom_point(aes(y=irf.lp, group=method, color=method, shape=method)) +
  theme_minimal() + 
  xlab('') + ylab('') +
  facet_wrap(~variable.name, scales='free_y') +
  theme( legend.position = 'bottom' )
ggsave(paste0(baseline.folder, 'monetary_compare_lp_naive_', n.lags, '_lags.pdf' ),
       height=8, width=8 )
ggsave(paste0(baseline.folder, 'monetary_compare_lp_naive_', n.lags, '_lags_wide.pdf' ),
       height=4, width=8 )

## 6.2 Rescale on initial value ##
ggplot(df.irf.rescale %>% filter(!is.na(variable.name)), aes(x=pd ) ) +
  geom_hline(aes(yintercept=0)) +
  xlab('') + ylab('') +
  geom_path(aes(y=irf.var.init, group=method, color=method, shape=method)) + 
  geom_point(aes(y=irf.var.init, group=method, color=method, shape=method)) +
  geom_ribbon( data=df.baseline, aes( ymin=baseline.init.pc.0.1, ymax=baseline.init.pc.0.9 ), alpha=.2 ) +
  geom_path(data=df.baseline, aes(y=baseline.init.pc.0.5, group='Baseline', 
                                  color='Baseline'), lwd=1.1 ) +
  geom_point(data=df.baseline, aes(y=baseline.init.pc.0.5, group='Baseline', 
                                   color='Baseline', shape='Baseline'), lwd=1.1 ) +
  scale_colour_manual(values=col.values) +
  scale_shape_manual(values=shape.values) +
  labs(color=NULL, shape=NULL) +
  theme_minimal() + 
  facet_wrap(~variable.name, scales='free_y') +
  theme( legend.position = 'bottom' )
ggsave(paste0(baseline.folder, 'monetary_compare_var_init_', n.lags, '_lags.pdf' ),
       height=8, width=8 )
ggsave(paste0(baseline.folder, 'monetary_compare_var_init_', n.lags, '_lags_wide.pdf' ),
       height=4, width=8 )

ggplot(df.irf.rescale %>% filter(!is.na(variable.name)), aes(x=pd ) ) +
  geom_hline(aes(yintercept=0)) +
  geom_path(aes(y=irf.lp.init, group=method, color=method, shape=method)) + 
  geom_point(aes(y=irf.lp.init, group=method, color=method, shape=method)) +
  geom_ribbon( data=df.baseline, aes( ymin=baseline.init.pc.0.1, ymax=baseline.init.pc.0.9 ), alpha=.2 ) +
  geom_path(data=df.baseline, aes(y=baseline.init.pc.0.5, group='Baseline', 
                                  color='Baseline'), lwd=1.1 ) +
  geom_point(data=df.baseline, aes(y=baseline.init.pc.0.5, group='Baseline', 
                                   color='Baseline', shape='Baseline'), lwd=1.1 ) +
  scale_colour_manual(values=col.values) +
  scale_shape_manual(values=shape.values) +
  labs(color=NULL, shape=NULL) +
  theme_minimal() + 
  xlab('') + ylab('') +
  facet_wrap(~variable.name, scales='free_y') +
  theme( legend.position = 'bottom' )
ggsave(paste0(baseline.folder, 'monetary_compare_lp_init_', n.lags, '_lags.pdf' ),
       height=8, width=8 )
ggsave(paste0(baseline.folder, 'monetary_compare_lp_init_', n.lags, '_lags_wide.pdf' ),
       height=4, width=8 )

## 6.3 Rescale on cumulative impact value ##
ggplot(df.irf.rescale %>% filter(!is.na(variable.name)), aes(x=pd ) ) +
  geom_hline(aes(yintercept=0)) +
  geom_path(aes(y=irf.var.cum, group=method, color=method, shape=method)) + 
  geom_point(aes(y=irf.var.cum, group=method, color=method, shape=method)) +
  geom_ribbon( data=df.baseline, aes( ymin=baseline.cum.pc.0.1, ymax=baseline.cum.pc.0.9 ), alpha=.2 ) +
  guides(color=NULL, shape=NULL) +
  geom_path(data=df.baseline, aes(y=baseline.cum.pc.0.5, group='Baseline', 
                                  color='Baseline'), lwd=1.1 ) +
  geom_point(data=df.baseline, aes(y=baseline.cum.pc.0.5, group='Baseline', 
                                   color='Baseline', shape='Baseline'), lwd=1.1 ) +
  scale_colour_manual(values=col.values) +
  scale_shape_manual(values=shape.values) +
  labs(color=NULL, shape=NULL) +
  theme_minimal() + 
  xlab('') + ylab('') +
  facet_wrap(~variable.name, scales='free_y') +
  theme( legend.position = 'bottom' )
ggsave(paste0(baseline.folder, 'monetary_compare_var_cum_', n.lags, '_lags.pdf' ),
       height=8, width=8 )
ggsave(paste0(baseline.folder, 'monetary_compare_var_cum_', n.lags, '_lags_wide.pdf' ),
       height=4, width=8 )

ggplot(df.irf.rescale %>% filter(!is.na(variable.name)), aes(x=pd ) ) +
  geom_hline(aes(yintercept=0)) +
  geom_path(aes(y=irf.lp.cum, group=method, color=method, shape=method)) + 
  geom_point(aes(y=irf.lp.cum, group=method, color=method, shape=method)) +
  geom_ribbon( data=df.baseline, aes( ymin=baseline.cum.pc.0.1, ymax=baseline.cum.pc.0.9 ), alpha=.2 ) +
  geom_path(data=df.baseline, aes(y=baseline.cum.pc.0.5, group='Baseline', 
                                  color='Baseline'), lwd=1.1 ) +
  geom_point(data=df.baseline, aes(y=baseline.cum.pc.0.5, group='Baseline', 
                                   color='Baseline', shape='Baseline'), lwd=1.1 ) +
  scale_colour_manual(values=col.values) +
  scale_shape_manual(values=shape.values) +
  labs(color=NULL, shape=NULL) +
  theme_minimal() + 
  xlab('') + ylab('') +
  facet_wrap(~variable.name, scales='free_y') +
  theme( legend.position = 'bottom' )
ggsave(paste0(baseline.folder, 'monetary_compare_lp_cum_', n.lags, '_lags.pdf' ),
       height=8, width=8 )
ggsave(paste0(baseline.folder, 'monetary_compare_lp_cum_', n.lags, '_lags_wide.pdf' ),
       height=4, width=8 )



