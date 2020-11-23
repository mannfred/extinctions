library(here)
library(tidyverse)

# extinctions per date

birds1 <-
  read.csv(here('Data/bird_extinctions_figure2.csv')) %>% 
  select(c(2, 4))

# dates of interest
years <- c(1790, 1852, 1900, 1950, 2000, 2019)

# count extinctions for dates of interest
reg_ext <- vector("numeric", 6L)

for (i in 1:length(years)) {
  reg_ext[i] <- 
    birds1 %>% 
    filter(Extinction.date < years[i]+1) %>% 
    tally() %>% 
    as.numeric()
}


# merge into data.frame
birds2 <- 
  data.frame("date" = years, "extinct" = reg_ext)

# ----------------------------
# linear model with no DEs

# extinctions = mu*year + b
# extinctions = 0.462*year - 792
model1 <- 
  lm(birds2$extinct ~ birds2$date)

mu1 <- as.numeric(model1$coefficients[2])
b1  <- as.numeric(model1$coefficients[1]) 


# ---------------------------
# add DEs from SEUX analysis by year
# data below obtained from SEUX output
# 1790 = 13.175 
# 1852 = 44.8126 - 13.175 = 31.6376
# 1900 = 55.5284 - 44.8126 = 10.7158
# 1950 = 56.3457 - 55.5284 = 0.8173
# 2000 = 56.359 - 56.3457 = 0.0133
# 2019 = 56.3615 - 56.359 = 0.0025

DEs <- c(13.175, 31.6376, 10.7158, 0.8173, 0.0133, 0.0025)
DEs_cumlative <- c(13.175, 44.8126, 55.5284, 56.3457, 56.359, 56.3615)

# add DEs to regular extinctions at stated years
tot_ext <- reg_ext + DEs_cumlative


# ----------------------------------
# fit linear model to dataset that includes DEs

birds_DEs <- 
  cbind(years, tot_ext) %>% 
  as_tibble() %>% 
  rename(date = years) %>% 
  rename(extinct = tot_ext) 


model2 <- 
  lm(birds_DEs$extinct ~ birds_DEs$date)


mu2 <- as.numeric(model2$coefficients[2])
b2  <- as.numeric(model2$coefficients[1]) 

# ---------------------------
# extrapolate so that there are 
# zero extinctions in year 1500

# b3 = -(mu2*1500) = -936.4972
f_pretax <- function(x) (mu2*x) - 910.373

# zero extinctions at year 1500
f_pretax(1500)


# -------------------------
# create new data points to compare the SEUX and pretax models

# for pretaxonomic DEs
# calculate extinctions for dates of interest
birds_pretax <-
  birds_DEs %>% 
  mutate(extinct = f_pretax(birds_DEs$date))




# -------------------------
# plotting

plot_data <-
  rbind(birds2, birds_DEs, birds_pretax) %>% 
  mutate(model = c(rep('observed extinctions', 6), rep('observed + dark extinctions', 6), rep('pre-taxonomic extinctions', 6)))

ggplot(plot_data) +
  aes(x = date, y = extinct, colour = model) +
  geom_point(size = 3) + 
  expand_limits(x = 1500) +
  theme_bw() + 
  # geom_smooth(method = 'lm') +
  stat_smooth(method = 'lm', fullrange = TRUE) +
  ylab('cumulative extinctions')


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


