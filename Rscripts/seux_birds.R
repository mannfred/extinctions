library(here)
library(seux)
library(tidyverse)


# data import ------------------
birds <- 
  read.csv(here('Data/bird_extinctions_cronk.csv'))

detections <- 
  data.frame(frstDetn = birds$"DESCRIPTION.DATE" ,
             lastDetn = birds$"Extinction.date") %>% 
  drop_na()

# try dropping pre-DNA descriptions (~1995)
detections <-
  detections %>% 
  filter(frstDetn <= 1995) %>% 
  mutate(
    lastDetn = case_when(
      lastDetn > 1995 ~ 1995,
      TRUE ~ as.numeric(lastDetn))) # this line leaves the column alone if it doesn't meet the above condition


# modelling -------------------

# d = discoveries
modelinputs <- 
  seux::get_model_inputs(
    detections$frstDetn, 
    detections$lastDetn,
    collapse_timesteps = TRUE)

CIs_estimates <- 
  seux::get_CI_estimate(modelinputs$S, modelinputs$E)
# CIs_estimates <- readRDS('birdanalysis.rds')

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
    theme_bw() #+
    #xlim(1500, 2020)
  
  return(p)
  
}

plot_output(df)

 write.csv(df, file = 'bird_DEs_oct23.csv')
