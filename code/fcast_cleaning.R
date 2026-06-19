##########################################################################################
#
# fcast_cleaning.R
#
# Functions to clean the forecasts for the "Identification from Rational Expectations" project
#
# Philip Barrett, Washington DC
# First version: 21jul2022
#
##########################################################################################

# Return: 1) Cleaned RE forecast, 2) marginal impact of f
library(glmnet)
library(caret)
library(rpart)
library(earth)

add_lags <- function(df, var.names, pd.name='date', n.lag=1){
# Adds lags to a data frame
  df.select <- df %>%
    select(pd.name, var.names)
  df.out <- df.select %>%
    mutate( this.pd=lead(get(pd.name), n.lag) ) %>%
    filter( !is.na(pd) ) %>%
    select(this.pd, everything(), -pd.name) %>%
    set_names( c( pd.name, paste0(var.names, '.l', n.lag) ) )
  return( left_join( df, df.out) )
}

fcast_ml <- function( df, target, data.names, target.fcast=NA, 
                         pd.name='date', fixedWindow=FALSE,
                         initialWindow=20, horizon=20, n.lag=8 ){
# Cleans the forecast for a single variable using lasso
  
  # Set up
  tuneLength.num <- 5
    # Combinations of tuning parameters to try
  
  # browser()
  
  if(is.na(target.fcast)) target.fcast <- paste0( target, '_f' )
  df.lasso <- df %>% 
    # select( all_of(pd.name), data.names, all_of(target.fcast) ) %>%
    mutate( target.lead=lead(get(target)) ) %>%  
    na.omit() %>%
    rename( pd=pd.name )
  all.vars <- setdiff( df %>% names(), pd.name )
    # The names of all the variables
  for( i in 1:n.lag ) df.lasso <- df.lasso %>% 
    add_lags( ., all.vars, 'pd', i)
    # add_lags( ., c( data.names, target.fcast), 'pd', i)
  df.lasso <- df.lasso %>% na.omit()
    # The data frame for the ML forecasting
    # TO DO: ADD LAGS
  
  # Set seeds for reproducibility
  set.seed(42)
  
  # Controls for the cross-validation
  myTimeControl <- trainControl(method = "timeslice",
                                initialWindow = initialWindow,
                                horizon = horizon,
                                fixedWindow = fixedWindow,
                                allowParallel = FALSE)
    # Default is to have at least 5 years of data, and use a rolling windows
    # with out-of-sample forecast evaluation over the next 5 years
  
  # Set up the regularized linear regression
  glmnet.mod <- suppressWarnings(
    train(target.lead ~ . - pd,
                      data = df.lasso,
                      method = "glmnet",
                      family = "gaussian",
                      trControl = myTimeControl,
                      tuneLength=tuneLength.num,
                      preProc = c("center", "scale"),
                      metric='RMSE') )
    # NB: lambda is the scalar on the penalty function and alpha is the
    # weighting on lasso vs. ridge regression (alpha=1 implies fully lasso)
  
  
  # Set up the standard linear regression
  lm.mod <- suppressWarnings(
    train(target.lead ~ . - pd,
                  data = df.lasso,
                  method = "lm",
                  trControl = myTimeControl,
                  tuneLength=tuneLength.num,
                  preProc = c("center", "scale"),
                  metric='RMSE') )
  
  # Set up the adaptive regression splines
  # earth.mod <- train(target.lead ~ . - pd,
  #                    data = df.lasso,
  #                    method = "earth",
  #                    trControl = myTimeControl,
  #                    tuneLength=tuneLength.num,
  #                    preProc = c("center", "scale"),
  #                    metric='RMSE')
  
  # Set up the regression tree
  rpart.mod <- suppressWarnings(
    train(target.lead ~ . - pd,
                     data = df.lasso,
                     method = "rpart",
                     trControl = myTimeControl,
                     tuneLength=tuneLength.num,
                     preProc = c("center", "scale"),
                     metric='RMSE') )
  
  # Just the raw forecast 
  raw.mod.rmse <- (df.lasso[[target.fcast]] - df.lasso$target.lead)^2 %>% mean %>% sqrt()
  
  # Summarize the different models
  l.mod <- list(glmnet = glmnet.mod,
                lm = lm.mod,
                # earth=earth.mod,
                rpart=rpart.mod)
  resamps <- resamples( l.mod )
  resamp.sumy <- resamps %>% summary()
  i.best.mod <- resamp.sumy$statistics$RMSE[,'Mean'] %>% which.min()
  best.mod <- l.mod[[i.best.mod]]
  best.pred <- predict( best.mod )
  
  v.rmse <- c( glmnet=glmnet.mod$results$RMSE %>% min,
               raw=raw.mod.rmse,
               lm=lm.mod$results$RMSE,
               rpart=rpart.mod$results$RMSE %>% min )
  
  df.out <- df.lasso %>% 
    select(pd) %>%
    cbind( best.pred ) %>%
    set_names( c( pd.name, paste0(target, '_re') ) ) %>%
    as.tibble()
      # The output predictions
  
  return( list( df=df.out, model=best.mod,  rmse=v.rmse,
                model.type=names(l.mod)[i.best.mod], resamp=resamps ) )
}

fcast_ml_all <- function( df, data.names, ... ){
# Applies the ML 
  
  l.fcast <- lapply( data.names, function(x) fcast_ml( df, x, data.names, ... ) )
  df.out <- lapply( l.fcast, function(x) x$df ) %>%
    reduce(full_join)
  model.types <- sapply( l.fcast, function(x) x$model.type )
  return( list( df=df.out, model.types=model.types, details=l.fcast ) )
}