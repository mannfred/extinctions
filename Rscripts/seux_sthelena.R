# to install SEUX:
# library(devtools)
# devtools::install_github("nadiahpk/seux", subdir="seux")

library(here)
library(seux)
library(tidyverse)

helena <- 
  read.csv(here('Data/sthelena_dark_extinction.csv'))

detections <- 
  data.frame(frstDetn = helena$frstDetn,
             lastDetn = helena$lastDetn) %>% 
  drop_na()

modelinputs <- 
  seux::get_model_inputs( 
    detections$frstDetn, detections$lastDetn)

CIs_estimates <- 
  seux::get_CI_estimate(modelinputs$S, modelinputs$E)

old_estimates <- 
  get_old_estimate(modelinputs$S, modelinputs$E )

df <- cbind(modelinputs, CIs_estimates, old_estimates)

plot_output <- function(df) {
  
  p <- ggplot(df) + 
    geom_line(aes(year, S,      color="S", linetype="data")) + 
    geom_line(aes(year, E,      color="E", linetype="data")) + 
    geom_line(aes(year, U_mean, color="U", linetype="hyper")) + 
    geom_line(aes(year, U_old,  color="U", linetype="old")) + 
    geom_line(aes(year, X_mean, color="X", linetype="hyper")) + 
    geom_line(aes(year, X_old,  color="X", linetype="old")) +
    scale_color_manual(name="Species class", 
                       values=c(
                         "S"     = "darkgreen", 
                         "E"     = "red",
                         "U"     = "orange",
                         "X"     = "blue"
                       ),
                       labels=c(
                         "S"     = expression(S[t]),
                         "E"     = expression(E[t]),
                         "U"     = expression(U[t]),
                         "X"     = expression(X[t])
                       )
    ) +
    scale_linetype_manual(name="Method of inferrence", 
                          values=c(
                            "data"  = "solid", 
                            "hyper" = "longdash",
                            "old"   = "dotted"
                          ),
                          labels=c(
                            "data"  = "from data",
                            "hyper" = "hypergeometric",
                            "old"   = "Chisholm et al. (2016)"
                          )
    ) +
    geom_ribbon(aes(x=year, ymin=X_lo, ymax=X_hi), fill="blue", alpha=0.3) +
    geom_ribbon(aes(x=year, ymin=U_lo, ymax=U_hi), fill="orange", alpha=0.3) +
    xlab('years') +
    ylab('number of species') +
    theme_bw()
  
  return(p)
  
}

plot_output(df)
