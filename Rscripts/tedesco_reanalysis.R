library(here)
library(stats4)
library(tidyverse)


# -------------------------------------------
## big function for running the Tedesco model
get_Tedesco_est = function(my_df, last_ext_year = NA){
  
  names(my_df) = c("det_year", "ext_year")
  
  ## Step 1: Getting extinction rate, mu
  max_T = max(my_df$ext_year)
  
  # creates an empty numeric object with length 169
  x_t = numeric(max_T)
  
  # find total extants for every year
  for ( i in 0:(max_T-1) ){
    x_t[i+1] = sum(my_df$det_year == i)
  }
  
  # number of years that extinctions were below their max #
  NE = sum(my_df$ext_year != max_T)
  

  if(is.na(last_ext_year)){last_ext_year = max_T}
  
  eq13 = function(mu){
    
    # Probability of non-persistence over all years
    i = seq(last_ext_year, 1, -1)
    p = 1-(1-mu)^i
    
    # Species are expected to be observed extinct if they do not persist until the tt-th year
    EE = 0
    for(j in 1:last_ext_year){
      EE = EE + x_t[j]*p[j]
    }
    
    diff_NE = abs(EE-NE)
    return(diff_NE)
  }
  
  fit = optimise(eq13, lower = 0, upper = 0.1, tol = 1e-25) 
  est_mu = fit$minimum
  
  ## Step 2: Finding N0 and nu
  library(stats4)
  eq6 = function(par){
    
    spp = x_t[-1]
    years = 1:(max_T-1)
    
    est_nu = par[1]
    Nm = par[2]
    ED = ((1-est_nu-est_mu)^(years-1))*est_nu*Nm
    
    # Negative log likelihood function
    nll = -sum(spp*(log(abs(ED)))-ED)
    return(nll)
  }
  
  est = nlminb(start = c(0.001,500), objective = eq6, 
               lower = c(1e-10, 1), upper = c(1-1e-10, Inf))
  
  ## Step 3: Entering fitted parameters back into equations
  time = seq(1,(max_T-1),1)
  est_nu = est$par[1]
  N0 = est$par[2]
  est_D = cumsum(c(0,(1-est_nu-est_mu)^(time-1)*est_nu*N0))
  est_E = cumsum(c(0,(1-est_nu-est_mu)^(time-1)*est_mu*N0))
  
  return(list(mu = est_mu,
              nu = est_nu,
              N0  = N0 + x_t[1],
              dets = est_D,
              UE = est_E,
              UT = N0+x_t[1] - est_D - est_E,
              pE = (tail(est_E,1)+NE)/(x_t[1]+N0)))
}



# -------------------------------------------
# reproduce Tedesco results:

file_names = Sys.glob(here("Data/lum_etal_2020_data/S4*.csv"))

# Start and end years as in Tedesco et al. (2014)
start_years = c(1979, 1805, 1840, 1825, 1855, 1842, 1894, 1758)
end_years = c(2007, 2010, 2008, 2010, 2009, 2010, 2010, 2005)

# build empty dataframe
all_results = data.frame(matrix(nrow = 8, ncol = 16))
names(all_results) = c("dataset", 
                       "Tedesco_mu", "Tedesco_mu_lo", "Tedesco_mu_hi",
                       "Tedesco_nu", "Tedesco_nu_lo", "Tedesco_nu_hi",
                       "Tedesco_N0", "Tedesco_N0_lo", "Tedesco_N0_hi",
                       "Tedesco_UE", "Tedesco_UE_lo", "Tedesco_UE_hi",
                       "Tedesco_pE", "Tedesco_pE_lo", "Tedesco_pE_hi")

# use file names to fill in "group_names"
# set substr to delete the first 82 chars and last 4 chars 
group_names = substr(file_names, 82, nchar(file_names)-4)
group_names = gsub("_", " ", group_names)
all_results[,1] = group_names

set.seed(100)
for (i in seq_along(file_names)){
  
  my_df = read.csv(file_names[i])
  
  # Rescaling description/extinction years for use with functions
  last_ext = max(my_df$extinction_year, na.rm = TRUE) - start_years[i]
  my_df$extinction_year[is.na(my_df$extinction_year)] = end_years[i]+1
  my_df$description_year = sapply(my_df$description_year - start_years[i], function(x) max(0,x))
  my_df$extinction_year = my_df$extinction_year - start_years[i]
  
  # For compatability with SEUX model
  my_df2 = my_df
  if(group_names[i] == "Australian mammals"){ 
    my_df2 = subset(my_df2, species != "Bettongia.pusilla")}
  if(group_names[i] == "World mammals"){ 
    my_df2 = subset(my_df2, species != "Bos.primigenius" & species != "Hippotragus.leucophaeus")}
  
  my_df2$description_year = ifelse(my_df2$description_year>my_df2$extinction_year, 
                                   my_df2$extinction_year-1, my_df2$description_year)
  
  # Running the model
  t_result = get_Tedesco_est(my_df[2:3], last_ext)
  
  # Bootstrapping to get Tedesco model CI
  temp_UE = temp_N0 = temp_pE = temp_mu = temp_nu = numeric(1000)
  for(jj in 1:1000){
    new_my_df = my_df[sample(nrow(my_df), nrow(my_df), replace = TRUE),]
    new_t_est = get_Tedesco_est(new_my_df[2:3], last_ext)
    temp_UE[jj] = tail(new_t_est$UE,1)
    temp_N0[jj] = new_t_est$N0
    temp_pE[jj] = new_t_est$pE
    temp_mu[jj] = new_t_est$mu
    temp_nu[jj] = new_t_est$nu
  }
  
  all_results[i,-1] = c(round(c(t_result$mu, quantile(temp_mu, c(0.025, 0.975))),7),
                        round(c(t_result$nu, quantile(temp_nu, c(0.025, 0.975))),5),
                        round(c(t_result$N0, quantile(temp_N0, c(0.025, 0.975))),0),
                        round(c(tail(t_result$UE,1), quantile(temp_UE, c(0.025,0.975))),1),
                        round(c(t_result$pE, quantile(temp_pE, c(0.025, 0.975))),3))
                        
}

# all_results agrees with Table S3 in Lum et al 2020


# -----------------------------------------------
# try to reproduce just one dataset at a time (amphibians)

my_df <- read.csv(here('Data/lum_etal_2020_data/S4_World_birds.csv'))

# Appendix S3 indicates data analysed for 1842-2010
# but data goes back to 1758
start_years <- 1842
end_years <- 2010

# build empty dataframe
all_results = data.frame(matrix(nrow = 1, ncol = 16))
names(all_results) = c("dataset", 
                       "Tedesco_mu", "Tedesco_mu_lo", "Tedesco_mu_hi",
                       "Tedesco_nu", "Tedesco_nu_lo", "Tedesco_nu_hi",
                       "Tedesco_N0", "Tedesco_N0_lo", "Tedesco_N0_hi",
                       "Tedesco_UE", "Tedesco_UE_lo", "Tedesco_UE_hi",
                       "Tedesco_pE", "Tedesco_pE_lo", "Tedesco_pE_hi")

all_results[,1] = "Australian Birds"

# Rescaling description/extinction years for use with functions
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
t_result = get_Tedesco_est(my_df[2:3], last_ext)

# works!
# Bootstrapping to get Tedesco model CI
set.seed(100)
temp_UE = temp_N0 = temp_pE = temp_mu = temp_nu = numeric(1000)
for(jj in 1:1000){
  new_my_df = my_df[sample(nrow(my_df), nrow(my_df), replace = TRUE),]
  new_t_est = get_Tedesco_est(new_my_df[2:3], last_ext)
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
