library(here)
library(tidyverse)


# data import
birds2 <- 
  read.csv(here('Data/figure2_data.csv')) %>% 
  slice(7:12) # capture data from 1800 - present


# ----------------------------
# linear model with no DEs

# extinctions = mu*year + b
# extinctions = 0.462*year - 792
model1 <- 
  lm(birds2$extinct ~ birds2$date)

mu1 <- as.numeric(model1$coefficients[2])
b1  <- as.numeric(model1$coefficients[1]) 


# ---------------------------
# add DEs from SEUX analysis
# and spread out evenly from 1800-2012

# add DE rate to mu1
mu2 <- mu1 + (55/(2012-1800))

# solve for intercept b2: 
# we know that by 2012 there should
# be 55 more extinctions than the OG dataset
# b2 = -1262.557
(max(birds2$extinct) + 55) - (mu2*2012)

f_seuxDE <- function(x) (mu2*x) - 1262.557

# create dataset from new model
birds_seuxDEs <-
  birds2 %>% 
  mutate(extinct = f_seuxDE(birds2$date))


# ---------------------------
# extrapolate so that there are 
# zero extinctions in year 1500

# b3 = -(0.721*1500) = -1081.5
f_pretax <- function(x) (mu2*x) - 1081.5

# zero extinctions at year 1500
f_pretax(1500)


# -------------------------
# create new data points to compare the SEUX and pretax models

# for pretaxonomic DEs
# calculate extinctions for dates of interest
birds_pretax <-
  birds2 %>% 
  mutate(extinct = f_pretax(birds2$date))

# -------------------------
# plotting

plot_data <-
  rbind(birds2, birds_seuxDEs, birds_pretax) %>% 
  mutate(model = c(rep('no DEs', 6), rep('SEUX DEs', 6), rep('SEUX DEs + Cronk Method', 6)))

ggplot(plot_data) +
  aes(x = date, y = extinct, colour = model) +
  geom_point(size = 3) + 
  expand_limits(x = 1500) +
  theme_bw() + 
  # geom_smooth(method = 'lm') +
  stat_smooth(method = 'lm', fullrange = TRUE)
# -----------------------------------------
# from SEUX we know there were 55 dark extinctions
# we need to add them while keeping extinctions at 1500 = 0
# so let's add them by increasing extinction rate

# cronk model
# extinctions = 0.462(2012) - 693 

# add 55 DEs while keeping ext = 0 at year 1500
# 0 = mu(1500) - 693 - 55
# mu = 0.4986667

# where x = year
f <- function(x) 0.4986667*x - 748

# calculate extinctions for dates of interest
tot_ext <- sapply(birds3$date, f)

# add new extinction data to dataset
birds4 <-
  birds3 %>% 
  mutate(extinct = tot_ext)

plot_data <-
  rbind(birds2, birds3, birds4) %>% 
  mutate(model = c(rep('no DEs', 6), rep('Cronk Method DEs', 6), rep('Cronk Method + SEUX DEs', 6)))

ggplot(plot_data) +
  aes(x = date, y = extinct, colour = model) +
  geom_point(size = 3) + 
  expand_limits(x = 1500) +
  theme_bw() + 
  # geom_smooth(method = 'lm') +
  stat_smooth(method = 'lm', fullrange = TRUE)





# ------------------------------------------------------
# let's calculate the slope ("extinction rate") manually
# using Quentin's bird data 

birds2 <- 
  read.csv(here('Data/figure2_data.csv')) %>% 
  slice(7:12) # capture data from 1800 - present

x <- birds2$date
y <- birds2$extinct
xhat <- mean(x)
yhat <- mean(y)

library(magrittr)

# slope (extinction rate) is 0.462
er <- 
  sum(
    (x - xhat) %>% 
      multiply_by(y - yhat)) %>% 
  
  divide_by(
    
    sum(
      (x - xhat) %>% 
        multiply_by(x - xhat))
  )


