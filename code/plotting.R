##########################################################################################
#
# plotting.R
#
# Plotting functions for the "Identification from Rational Expectations" project
#
# Philip Barrett, Washington DC
# First version: 20jul2022
#
##########################################################################################

library(ggh4x)

number_ticks <- function(n) {function(limits) pretty(limits, n)}

irf_plot_all <- function( df.in, data.names, fcast.names, data.labels, 
                          shk.names=NA, sd.mult=NA, trim=FALSE, 
                          lin.approx=0, pc=c(.05, .2, .8, .95 ), do.pt=FALSE ){
# Makes a ggplot object plotting  the IRFs
  
# lin.approx is a value in [0,0.5] which determines how much of the tails of the
# distribution to throw away when computing the CE.  For lin.approx==0, none is
# thrown away, and in the limit where lin.approx==0.5, the whole tails are. 
  
  # browser()
  if(is.na(sd.mult)){
    # browser()
    v.bounds <- paste0( 'pc.', c(pc, .5) ) %>% 
      set_names( c( paste0( 'upper', c('', '.1') ),
                    paste0( 'lower', c( '.1', '') ), 'median' ) )
    # lin.ratio <- if(lin.approx==0) 1 else qnorm(1-(1-ce)/2)  / qnorm(1-lin.approx/2)
      # Scalar for the approximate CIs
    # v.rename <- c( 'Point estimate'='pt.est', 'Median'='pc.0.5', v.bounds )
    if(do.pt) v.bounds['median'] <- 'pt'
      # Make the point estimate the "median" for plotting
    df.cht <- df.in %>% 
      select( pd, type, dimension, pc, data.names, fcast.names, shk.name, irf.sign ) %>%
      gather(variable, value, -pd, -type, -dimension, -pc, -shk.name, -irf.sign ) %>%
      separate(variable, c('variable','fcast'), sep='_' ) %>%
      mutate(fcast=(fcast=='f') %>% replace_na(FALSE),
             # shk=if(is.na(shk.names)) paste0(type,'.',dimension) else get(shk.names),
             type.cht=case_when(type=='u'~' Surprise',
                                type=='v'~'News',
                                type=='epsilon'~'Combined'),
             variable.cht=data.labels[variable]) %>%
      spread( pc, value ) %>%
      rename( all_of(v.bounds) )
  }else{
    # browser()
    df.cht <- df.in %>% 
      select( pd, type, dimension, pc, data.names, fcast.names, shk.name, irf.sign ) %>%
      gather(variable, value, -pd, -type, -dimension, -pc, -shk.name, -irf.sign ) %>%
      separate(variable, c('variable','fcast'), sep='_' ) %>%
      mutate(fcast=(fcast=='f') %>% replace_na(FALSE),
             type.cht=case_when(type=='u'~' Surprise',
                                type=='v'~'News',
                                type=='epsilon'~'Combined'),
             variable.cht=data.labels[variable] ) %>%
      spread( pc, value ) %>%
      mutate( upper=irf.sign*(irf.sign*pt.est+sd.mult*(if(trim) sd.trim else sd)), 
              lower=irf.sign*(irf.sign*pt.est-sd.mult*(if(trim) sd.trim else sd)) )
  }
  
  gg.out <- ggplot( df.cht %>% filter(!fcast), aes(x=pd-1) ) +
                # geom_ribbon( aes(ymin=lower, ymax=upper), fill = "grey70") +
                ylab('Percent') +
                xlab('Quarters') +
                facet_nested( variable.cht ~ shk.name+type.cht, scales='free_y') +
                theme_minimal() +
                scale_x_continuous(breaks=number_ticks(6)) +
                theme(legend.position = 'bottom') +
                geom_hline(aes(yintercept=0), lwd=.1)
  
  
  if( is.na(sd.mult)){
    gg.out <- gg.out +
      # geom_line(aes(y=irf.sign*upper.1), lty=3, lwd=.6) +
      # geom_line(aes(y=irf.sign*lower.1), lty=3, lwd=.6) +
      geom_line(aes(y=irf.sign*upper), lty=2, lwd=.5) +
      geom_line(aes(y=irf.sign*lower), lty=2, lwd=.5) +
      geom_line(aes(y=irf.sign*median), lty=1, lwd=.5)
  }else{
    gg.out <- gg.out +
      geom_line(lwd=.5, aes(y=irf.sign*pc.0.5) ) + #, color='grey20') +
      geom_line(aes(y=irf.sign*upper), lty=2, lwd=.6) +
      geom_line(aes(y=irf.sign*lower), lty=2, lwd=.6)
  }
  
  if('pt.est' %in% (df.cht %>% names())){
    gg.out <- gg.out + geom_line(lwd=.9, aes(y=irf.sign*pt.est)) +
      geom_point( aes(y=irf.sign*pt.est), shape=16)
  }else{
    gg.out <- gg.out +
      geom_line(aes(y=irf.sign*median), lwd=.9 )
  }
  
  return(gg.out)
}


irf_plot_pt <- function( df.in, data.names, fcast.names, data.labels=NA ){
# Makes a ggplot object plotting  the IRFs

  if(all(is.na(data.labels))) data.labels <- data.names %>% set_names(data.names)
  
  df.cht <- df.in %>% 
    select( pd, type, dimension, data.names, fcast.names, shk.name, irf.sign ) %>%
    gather(variable, value, -pd, -type, -dimension, -shk.name, -irf.sign ) %>%
    separate(variable, c('variable','fcast'), sep='_' ) %>%
    mutate(fcast=(fcast=='f') %>% replace_na(FALSE),
           type.cht=case_when(type=='u'~'Surprise',
                              type=='v'~'News',
                              type=='epsilon'~'Combined'),
           variable.cht=data.labels[variable] )
  
  gg.out <- ggplot( df.cht %>% filter(!fcast), aes(x=pd-1) ) +
    geom_line(lwd=.9, aes(y=value*irf.sign)) +
    # geom_point( aes(y=value), shape=1) +
    theme_minimal() +
    ylab('Percent') +
    xlab('Quarters') +
    facet_nested( variable.cht ~ shk.name+type.cht, scales='free_y') +
    scale_x_continuous(breaks=number_ticks(6)) +
    theme(legend.position = 'bottom') +
    geom_hline(aes(yintercept=0), lwd=.1)
  
  return(gg.out)
}

irf_plot_pt_overlay <- function( gg.in, df.in, data.names, fcast.names, data.labels=NA ){
# Adds the point IRFs to a pre-existing plot
  
  if(all(is.na(data.labels))) data.labels <- data.names %>% set_names(data.names)
  df.cht <- df.in %>% 
    select( pd, type, dimension, data.names, fcast.names, shk.name, irf.sign ) %>%
    gather(variable, value, -pd, -type, -dimension, -shk.name, -irf.sign ) %>%
    separate(variable, c('variable','fcast'), sep='_' ) %>%
    mutate(fcast=(fcast=='f') %>% replace_na(FALSE),
           type.cht=case_when(type=='u'~'Surprise',
                              type=='v'~'News',
                              type=='epsilon'~'Combined'),
           variable.cht=data.labels[variable] ) %>%
    filter(!fcast)
  
  gg.out <- gg.in +
    geom_line(data=df.cht, aes(x=pd-1, y=value*irf.sign), lwd=1, lty=2, color='blue') #+
    # geom_point(data=df.cht, aes(x=pd, y=value*irf.sign), pch=1 )
  
  return(gg.out)
}

irf_plot_compare_pt <- function( df.1, df.2, data.names, fcast.names, 
                                 irf.names=c('IRF 1', 'IRF 2'), data.labels=NA ){
  # Makes a ggplot object plotting  the IRFs
  
  if(is.na(data.labels)) data.labels <- data.names %>% set_names(data.names)
  
  df.cht <- full_join( df.1 %>% mutate( irf=irf.names[1]), df.2 %>% mutate( irf=irf.names[2]) ) %>% 
    select( pd, type, dimension, data.names, fcast.names, shk.name, irf.sign, irf ) %>%
    gather(variable, value, -pd, -type, -dimension, -shk.name, -irf.sign, -irf ) %>%
    separate(variable, c('variable','fcast'), sep='_' ) %>%
    mutate(fcast=(fcast=='f') %>% replace_na(FALSE),
           type.cht=case_when(type=='u'~'Surprise',
                              type=='v'~'News',
                              type=='epsilon'~'Combined'),
           variable.cht=data.labels[variable] )
  
  gg.out <- ggplot( df.cht %>% filter(!fcast), aes(x=pd-1, lty=irf) ) +
    geom_line(lwd=.9, aes(y=value*irf.sign)) +
    # geom_point( aes(y=value), shape=1) +
    theme_minimal() +
    ylab('Percent') +
    xlab('Quarters') +
    labs(lty='') +
    facet_nested( variable.cht ~ shk.name+type.cht, scales='free_y') +
    scale_x_continuous(breaks=number_ticks(6)) +
    theme(legend.position = 'bottom') +
    geom_hline(aes(yintercept=0), lwd=.1)
  
  return(gg.out)
}

irf_plot_compare_pt_multi <- function( df,  group.name, data.names,
                                 data.labels=NA ){
  # Makes a ggplot object plotting  the IRFs
  
  if(any(is.na(data.labels))) data.labels <- data.names %>% set_names(data.names)
  
  df.cht <- df %>% 
    select( pd, type, dimension, data.names, shk.name, irf.sign, irf=group.name ) %>%
    gather(variable, value, -pd, -type, -dimension, -shk.name, -irf.sign, -irf ) %>%
    # separate(variable, c('variable','fcast'), sep='_' ) %>%
    mutate(type.cht=case_when(type=='u'~' Surprise',
                              type=='v'~'News',
                              type=='epsilon'~'Combined'),
           variable.cht=data.labels[variable] )
  
  gg.out <- ggplot( df.cht, aes(x=pd-1, y=value*irf.sign, group=irf, color=irf, shape=irf) ) +
    geom_line(lwd=.9) +
    geom_point( size=1.5 ) +
    theme_minimal() +
    ylab('Percent') +
    xlab('Quarters') +
    labs(color='', shape='') +
    facet_nested( variable.cht ~ shk.name+type.cht, scales='free_y') +
    scale_x_continuous(breaks=number_ticks(6)) +
    theme(legend.position = 'bottom') +
    geom_hline(aes(yintercept=0), lwd=.1)
  
  return(gg.out)
}
