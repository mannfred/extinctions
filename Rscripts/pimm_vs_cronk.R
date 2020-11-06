library(here)
library(magrittr)
library(tidyverse)


# ------------------------------------------------------
# let's calculate the slope ("extinction rate") manually
# using nadiah kristensen's example data from
# https://nadiah.org/2017/05/10/comparing-e-msy-and-chisholm-method/ 



nadia_data <- 
  read.csv(here("Data/nadia_data.csv"))

year <- nadia_data$year
extinct <- nadia_data$extinct


# cronk method --------------------
# regression equations from:
# https://stats.libretexts.org/Bookshelves/Introductory_Statistics/Book%3A_Introductory_Statistics_(Shafer_and_Zhang)/10%3A_Correlation_and_Regression/10.04%3A_The_Least_Squares_Regression_Line

model2 <- lm(extinct ~ year)

# summary(model2)
# plot(extinct ~ year)
# abline(model2)

x <- year
y <- extinct
xhat <- mean(year)
yhat <- mean(extinct)


# slope (extinction rate) is 31.714
er <- 
  sum(
    (x - xhat) %>% 
      multiply_by(y - yhat)) %>% 
  
  divide_by(
    
    sum(
      (x - xhat) %>% 
        multiply_by(x - xhat))
  )

# intercept is -15.1428
b_0 <- yhat - er*xhat

# pimm method -------------------

sum(diff(extinct)) / sum(nadia_data$extant[1:6])
