# to install SEUX:
# library(devtools)
# devtools::install_github("nadiahpk/seux", subdir="seux")

library(here)
library(seux)
library(tidyverse)

# QUESTIONS -------------------
# why does adding in earlier discovery dates (eg 1500) not 
# change anything in 'modelinputs'?


# data import ------------------
helena <- 
  read.csv(here('Data/sthelena_dark_extinction.csv'))

detections <- 
  data.frame(frstDetn = helena$frstDetn,
             lastDetn = helena$lastDetn) 
  #drop_na()



# experiments -----------------

# try making some spp go extinct earlier
detections[10:20, 2] <- c(1890, 1890, 1890, 1890, 1890, 1950, 1950, 1950, 1950, 1999, 1999)

# try making some spp be discovered later
detections[40:48, 1] <- c(1870, 1880, 1890, 1900, 1910, 1920, 1930, 1940, 1950)

# try forcing model to reach back to 1500
# by adding in an extinction event bw 1500 and 1501
# warning! cannot have gaps in the data (eg between 1500 and the next discovery)
detections[1,] <- c(1450, 1500)
detections[2,] <- c(1500, 2020)

# modelling -------------------

# d = discoveries
modelinputs <- 
  seux::get_model_inputs(
    detections$frstDetn, 
    detections$lastDetn,
    collapse_timesteps = TRUE)

CIs_estimates <- 
  seux::get_CI_estimate(modelinputs$S, modelinputs$E)

old_estimates <- 
  get_old_estimate(modelinputs$S, modelinputs$E )

df <- cbind(modelinputs, CIs_estimates, old_estimates)

plot_output <- function(df) {
  
  p <- ggplot(df) + 
    geom_line(aes(year, S,      color="S", linetype="data"), size = 2) + 
    geom_line(aes(year, E,      color="E", linetype="data"), size = 2) + 
    geom_line(aes(year, U_mean, color="U", linetype="hyper"), size = 2) + 
    geom_line(aes(year, U_old,  color="U", linetype="old"), size = 2) + 
    geom_line(aes(year, X_mean, color="X", linetype="hyper"), size = 2 ) + 
    geom_line(aes(year, X_old,  color="X", linetype="old"), size = 2 ) +
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
    theme_bw() +
    xlim(1500, 2020)
  
  return(p)
  
}

plot_output(df)
