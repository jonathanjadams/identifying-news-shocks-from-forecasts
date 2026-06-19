##########################################################################################
#
# var_decomp_compare.R
#
# Code to compare the variance decomposition across specifications
#
# Philip Barrett, Washington DC
# First version: 01aug2023
#
##########################################################################################


#### 0. Set up ####
rm(list=ls())
setwd('C:\\Users\\pbarrett\\OneDrive - International Monetary Fund (PRD)\\identification from RE/')
library(tidyverse)
source('code/structural_identification.R')

#### 1. Controls ####
# out_prefix <- '/28jul23/'
out_prefix <- '/final_replication/'
n.pds.irf <- 24
p.vals <- c( .025, .05, .1, .25 )

l.var.names <- list( 
  baseline_ml_1_lag=list( name='  Baseline'),
  baseline_ml_1_lag_short=list( name=' Short sample'),
  baseline_ml_2_lag=list( name=' 2 Lags'),
  baseline_ml_4_lag=list( name=' 4 lags'),
  baseline_ml_1_lag_spf=list(name=' SPF interest forecasts'),
  baseline_f_1_lag=list( name=' Literal forecasts'),
  fed_ml_1_lag=list( name='Inc. Federal spending'),
  unemp_ml_1_lag=list( name='Inc. Unemployment'),
  indpro_ml_1_lag=list( name='Inc. Industrial production'),
  oil_ml_1_lag=list( name='Inc. Oil Prices'),
  usdcadfx_ml_1_lag=list( name='Inc. Exch. rate')
) 

data.labels <- c( ir03='3-month \ninterest rate', 
                  ir5='5-year \ninterest rate', 
                  un='Unemployment', infl='Inflation', 
                  y=' Output', ir5='Five-year interest rate', 
                  sur='Surplus (percent of GDP)', tax=' Taxes',
                  oil='Oil Price', usdcad='USD-CAD exchange rate',
                  hs='Housing starts', indpro='Industrial production',
                  fed='  Federal government \nconsumption',
                  gov='  Government \nspending' )
chart.drop <- c('fed', 'un', 'indpro', 'oil', 'usdcad')

#### 2. Assemble the main dataset ####
l.df <- list()
i.counter <- 0


for( this.vars in names(l.var.names) ){
  
  i.counter <- i.counter + 1
  
  fig.dir <- paste0('figures',out_prefix, this.vars)
  load( paste0(fig.dir, '/est.rdata') )
  
  df.var.decomp.boot.frac.total <- df.var.decomp.boot.frac %>% 
    filter(pd==n.pds.irf+1) %>%
    filter(type!='epsilon') %>%
    full_join( {.} %>% group_by(variable, pd, type, i.boot ) %>% 
                 summarise( variance.pc=sum(variance.pc) ) %>% 
                 mutate( shk.name='Total') ) %>%
    filter(shk.name=='Total')
  
  df.var.decomp.pc.frac.total <- 
    lapply( c( p.vals, .5, rev(1-p.vals ) ), function(x) 
      df.var.decomp.boot.frac.total %>%
        group_by( type, shk.name, variable, pd ) %>%
        summarise_if( is.numeric, quantile, probs=x, na.rm=T ) %>%
        mutate( pc=paste0('pc.', x ) )
    ) %>%
    reduce(rbind) %>%
    select( type, shk.name, pd, pc, everything(), -dimension, -variance, -i.boot, -total ) %>%
    spread(pc,variance.pc) %>% mutate( model=l.var.names[[i.counter]]$name )
  
  l.df[[this.vars]] <- df.var.decomp.pc.frac.total
}

df.var.decomp.all <- l.df %>% reduce(rbind) %>%
  mutate( var.label=data.labels[variable],
          category=ifelse(type=='contemp', 'Surprise', 'News') )

pd <- position_dodge(.8)
ggplot(df.var.decomp.all %>% filter( !(variable %in% chart.drop)), aes(x=model, col=category ) ) +
  geom_point( aes(y=pc.0.5), position=pd, size=2 ) +
  geom_errorbar( aes(ymin=pc.0.1, ymax=pc.0.9), position=pd, width=.1, lwd=.9 ) +
  facet_wrap(~var.label, ncol=2 ) +
  ylim(c(0,100)) +
  ylab('Fraction of variance') +
  xlab('') +
  theme_minimal() +
  labs(color='') +
  theme(legend.position='bottom', axis.text.x=element_text(angle=90, vjust=0, hjust=1 ) )

ggsave( paste0('figures',out_prefix, 'var_decomp_compare.pdf'), height=10, width=8 )
ggsave( paste0('figures',out_prefix, 'var_decomp_compare_wide.pdf'), height=6, width=8 )

