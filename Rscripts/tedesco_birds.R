library(here)
library(stats4)
library(tidyverse)

# get_Tedesco_est() function is written in "tedesco_reanalysis.R"

# data import
birds <- 
  read.csv(here('Data/bird_extinctions_cronk.csv'))


detections <- 
  data.frame(description_year = birds$"DESCRIPTION.DATE" ,
             extinction_year = birds$"Extinction.date") %>% 
   drop_na() %>%
  mutate(extinction_year = na_if(extinction_year, 2019))

# build empty dataframe
all_results = data.frame(matrix(nrow = 1, ncol = 16))
names(all_results) = c("dataset", 
                       "Tedesco_mu", "Tedesco_mu_lo", "Tedesco_mu_hi",
                       "Tedesco_nu", "Tedesco_nu_lo", "Tedesco_nu_hi",
                       "Tedesco_N0", "Tedesco_N0_lo", "Tedesco_N0_hi",
                       "Tedesco_UE", "Tedesco_UE_lo", "Tedesco_UE_hi",
                       "Tedesco_pE", "Tedesco_pE_lo", "Tedesco_pE_hi")

all_results[,1] = "World Birds: Cronk Dataset"

# Tedesco uses date of first extinction (sometimes subtracting 10)
# our first extinction is 1784
start_years <- 
  detections %>%  
  drop_na() %>% 
  select(2) %>% 
  min() %>% 
  as.numeric()

end_years   <- 2018

# tedesco start/end years
start_years <- 1842
end_years <- 2010

# hi
my_df <- detections

# last_ext takes the most recent extinction and subtracts the start year 1985-1979 = 6. 
# this way we start at year 0 (instead of 1979), and the last extinction is at year 6
# instead of 1985
last_ext = max(my_df$extinction_year, na.rm = TRUE) - start_years

# for all NAs (where extinctions haven't yet been observed), input "2007+1" which is "2008"
# effectively setting extinction year to 2008
my_df$extinction_year[is.na(my_df$extinction_year)] = end_years+1

# shifting extinction years to accomodate for setting start year as 0
my_df$extinction_year = my_df$extinction_year - start_years

# subtracting 1979 from description year
# if result is negative, then replace with year = 0
my_df$description_year = sapply(my_df$description_year - start_years, function(x) max(0,x))


# if a species was described after it went extinct, 
# then subtract 1 year from the extinction year until fixed
# otherwise, leave unaltered
my_df$description_year = ifelse(my_df$description_year > my_df$extinction_year, 
                                my_df$extinction_year-1, my_df$description_year)

# works once abs() added to get_Tedesco_est()
t_result = get_Tedesco_est(my_df, last_ext)

# Bootstrapping to get Tedesco model CI
set.seed(100)
temp_UE = temp_N0 = temp_pE = temp_mu = temp_nu = numeric(1000)
for(jj in 1:1000){
  new_my_df = my_df[sample(nrow(my_df), nrow(my_df), replace = TRUE),]
  new_t_est = get_Tedesco_est(new_my_df[1:2], last_ext) # changed from [2:3]
  temp_UE[jj] = tail(new_t_est$UE,1)
  temp_N0[jj] = new_t_est$N0
  temp_pE[jj] = new_t_est$pE
  temp_mu[jj] = new_t_est$mu
  temp_nu[jj] = new_t_est$nu
}

all_results[1,-1] = c(round(c(t_result$mu, quantile(temp_mu, c(0.025, 0.975))),7),
                      round(c(t_result$nu, quantile(temp_nu, c(0.025, 0.975))),5),
                      round(c(t_result$N0, quantile(temp_N0, c(0.025, 0.975))),0),
                      round(c(tail(t_result$UE,1), quantile(temp_UE, c(0.025,0.975))),1),
                      round(c(t_result$pE, quantile(temp_pE, c(0.025, 0.975))),3))

