library(here)
library(tidyverse)

birds2 <- 
  read.csv(here('Data/figure2_data.csv')) %>% 
  slice(7:12) # capture data from 1800 - present

plot(birds2)

# first model: extinctions = 0.462(year) - 792.454
# where B_1 = -792.454
model1 <- 
  lm(birds2$extinct ~ birds2$date)


# extrapolate year when extinctions = 0
# year = 1715
792.454/0.462

# now force extinctions = 0 for year 1500
# 0 = 0.462(1500) - B_2
# B_2 = -693
# B_1 - B_2 = dB
# dB = 792.454 - 693 = 99.454
0.462*1500 

# with the adjusted model, how many extinctions at year 2012?
# extinctions = 0.462(2012) - 693
# total extinctions = 236
# observed extinctions = 134
# dark extinctions = 236 - 134 = 102
0.462*2012 - 693

# now modify data to give 0 extinct at year 1500
birds3 <-
  birds2 %>% 
  mutate(extinct = extinct + 102)

# plotting
plot_data <-
  full_join(birds2, birds3) %>% 
  mutate(type = c(rep('no DEs', 6), rep('DEs added', 6)))


mycolours <- c(rep("#CC79A7", 6), rep("#56B4E9", 6))


ggplot(plot_data) +
  aes(x = date, y = extinct, colour = type) +
  geom_point(size = 3) + 
  expand_limits(x = 1500) +
  theme_bw() + 
  # geom_smooth(method = 'lm') +
  stat_smooth(method = 'lm', fullrange = TRUE)

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


