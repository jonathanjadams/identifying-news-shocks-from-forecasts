##########################################################################################
#
# make_ml_fcast.R
#
# Code to make the set of all machine leraning forecasts
#
# Philip Barrett, Washington DC
# First version: 28feb2023
#
##########################################################################################


#### 0. Set up ####
rm(list=ls())
library(tidyverse)
source('code/fcast_cleaning.R')
select <- dplyr::select
# Prevent conflict with MASS::select

#### 1. Controls ####
data.in.file <- 'data/final_replication_data.csv'
mult.100.vars <- c('y', 'infl', 'gov', 'sur', 'fed', 'tax' )
ml.fcast.file <- 'data/ml_fcast.rdata'
n.lag.ml <- 8 
fig.dir <- 'figures/00_forecasts'
dir.create(fig.dir, recursive=TRUE, showWarnings=FALSE)

### 2.0 Read the raw data ###
df.in <- data.in.file %>% read_csv
df <- df.in %>% 
  separate( quarter, c('year', 'quarter'), sep = 'q') %>%
  mutate( quarter=quarter %>% as.numeric,
          date=make_date( year, 3*quarter, 1 ) ) %>%
  select( date, everything(), -year, -quarter ) %>% 
  rename_with(~str_remove(., '_dt')) %>%
  filter(!is.na(date)) %>%
  mutate_if( is.numeric, zoo::na.approx, rule=2 ) %>%
  mutate_if( is.numeric, function(x) x-mean(x,na.rm=T) ) %>%
  mutate_at(vars(contains(mult.100.vars)), function(x) x*100 )
data.names <- names(df)[ !(names(df) %>% str_detect('_f'))] %>%
  setdiff( 'date' )

### 2.1 Clean the forecasts ###
fcast.clean <- fcast_ml_all( df, data.names, n.lag=n.lag.ml )
  # Need to summarize how much info is coming from the forecasts here.
df.f.clean <- fcast.clean$df %>%
  set_names( c('date', paste0( data.names, '_z' ) ) )

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

save( df.f.clean, df.fcast.sumy, file=ml.fcast.file )
  # Save the output