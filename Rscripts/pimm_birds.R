library(here)
library(tidyverse)

# data import ------------------
birds <- 
  read.csv(here('Data/bird_extinctions.csv'))

head(birds)


