# from:
# https://conbio.onlinelibrary.wiley.com/doi/epdf/10.1111/cobi.13562 
# for Tedesco model, see line 347


## S5 - functions for running simulations and analysis

## Functions ============================================================================

## Appendix Equation (1) : Probability of extinction at time t given starting abundance n_0
P_T = function(lambda, p_sigma, n_0, t_star){
  exp_term = exp(lambda*p_sigma*t_star)
  output1 = ((1-exp_term)/(1-p_sigma-exp_term))^n_0
  return(output1)
}

## Appendix Equation (7) : Generating function for species abundance
n_gen_func = function(x, lambda, omega, p_sigma, t, n_0){
  lambdastar = lambda*(1-p_sigma)
  num = omega*(x-1)*exp(lambdastar*t)+(omega-lambdastar*x)*exp(omega*t)
  denom = lambdastar*(x-1)*exp(lambdastar*t)+(omega-lambdastar*x)*exp(omega*t)
  gen_func = ((num/denom))^n_0
  return(gen_func)
}

## Appendix Equation (11) : Generating function for number of detections
dets_gen_func = function(y, lambda, omega, xi, p_sigma, t, n_0){
  lambdastar = lambda*(1-p_sigma)
  sqrt_factor = sqrt(as.complex(-(lambdastar^2)+2*lambdastar*(omega-xi*(1-y))-((omega+xi*(1-y))^2)))
  gen_func_a = (1/(2*lambdastar))*(lambdastar+omega+xi*(1-y)+
                                     sqrt_factor*tan((1/2)*t*sqrt_factor+atan((lambdastar-omega-xi*(1-y))/sqrt_factor)))
  gen_func = gen_func_a^n_0
  return(gen_func)
}

## Bornemann (2011) algorithm to extract mu probabilities
freqdbn = function(n, lambda, omega, p_sigma, time, n_0){
  r = 1
  fac = exp(-n*log(r))
  cauchy = function(t){
    return(fac*(exp(-n*t)*(n_gen_func(r*exp(t), lambda, omega, p_sigma, time, n_0))))
  }
  m = max(n+1,8)
  tol = 10^(-10)
  im = complex(real=0,imaginary=1)
  s = cauchy(2*im*pi*(1:m)/m)
  val1 = mean(s)
  err1 = NaN
  n_step = 0
  while(m<=2^20){
    # Algorithm is forced to run for a minimum of n>=n_steps
    m = 2*m
    s = c(s,cauchy(2*im*pi*seq(1,m,2)/m))
    val = mean(s)
    kappa = mean(abs(s))/abs(val)
    err0 = abs(val-val1)/abs(val)
    err = ((err0/err1)^2)*err0
    if(!is.na(err1) & n_step >= 4){
      if(err<=kappa*tol|!is.finite(kappa)){
        m = 2^20+1
      }
    }
    val1 = val
    err1 = err0
    n_step = n_step + 1
    if(n>n_0){n_step = 4}
  }
  if(!is.na(val1)){
    return(abs(Re(val1)))
  }
}

## Bornemann (2011) algorithm to extract nu probabilities
freqdbn_2 <- function(k, lambda, omega, xi, p_sigma, time, n_0){
  r = 1
  fac = exp(-k*log(r))
  cauchy = function(t){
    return(fac*(exp(-k*t)*(dets_gen_func(r*exp(t), lambda, omega, xi, p_sigma, time, n_0))))
  }
  m = max(k+1,8)
  tol = 10^(-10)
  im = complex(real=0,imaginary=1)
  s = cauchy(2*im*pi*(1:m)/m)
  val1 = mean(s)
  err1 = NaN
  n_step = 0
  while(m<=2^20){
    # Algorithm is forced to run for a minimum of n>=n_steps
    m = 2*m
    s = c(s,cauchy(2*im*pi*seq(1,m,2)/m))
    val = mean(s)
    kappa = mean(abs(s))/abs(val)
    err0 = abs(val-val1)/abs(val)
    err = ((err0/err1)^2)*err0
    if(!is.na(err1) & n_step >= 4){
      if(err<=kappa*tol|!is.finite(kappa)){
        m = 2^20+1
      }
    }
    val1 = val
    err1 = err0
    n_step = n_step + 1
    if(k>n_0){n_step = 4}
  }
  if(!is.na(val1)){
    return(abs(Re(val1)))
  }
}

## Probability Q
prob_Q <- function(n, lambda, omega, p_sigma, n_0, t){
  
  # manually set probabilities at t=0 to avoid unintended behaviour
  if(t==0){
    if(n==n_0){
      Q = 1
    } else {
      Q = 0
    }
  } else {
    Q = try(freqdbn(n, lambda, omega, p_sigma, t, n_0),silent=TRUE)
  }
  return(Q)
}

## Probability R
prob_R <- function(k, lambda, omega, xi, p_sigma, n_0, t){
  R = try(freqdbn_2(k, lambda, omega, xi, p_sigma, t, n_0),silent=TRUE)
  return(R)
}

## Appendix Equation (2) : Probability of extinction at time t conditional on survival up to t
get_simulation_mus = function(lambda, p_sigma, n_0, t, DeltaT){
  mu = (P_T(lambda, p_sigma, n_0, t*DeltaT) - P_T(lambda, p_sigma, n_0, (t-1)*DeltaT))/
    (1 - P_T(lambda, p_sigma, n_0, (t-1)*DeltaT))
  return(mu)
}

## Appendix Equation (3) : Probability of detection at time t conditional on survival up to t
get_simulation_nus = function(lambda, omega, xi, p_sigma, n_0, t, DeltaT, k = 0){
  num = 0
  sum_num1 = prob_Q(0, lambda, omega, p_sigma, n_0, (t-1)*DeltaT)
  n = 1
  while(sum_num1<0.9999){
    xi_n = xi/(n^k) # Abundance-specific per-capita detection rate
    num1 = prob_Q(n, lambda, omega, p_sigma, n_0, (t-1)*DeltaT)
    num2 = 1-prob_R(0, lambda, omega, xi_n, p_sigma, n, DeltaT)
    sum_num1 = sum_num1 + num1
    num = num + num1*num2
    n = n+1
  }
  nu = num/(1-prob_Q(0, lambda, omega, p_sigma, n_0, (t-1)*DeltaT))
  return(nu) 
}

## Simulate species abundances
get_spp_abundances = function(J, p_sigma){
  
  # Expected number of species with abundance n (log-series SAD)
  S_n = function(J, p_sigma, n){
    J*p_sigma/(1-p_sigma)/n*(1-p_sigma)^n
  }
  flag = 0
  n = 0
  while(flag == 0){
    n = n+1
    if (round(S_n(J,p_sigma,n)) == 0){flag = 1}
  }
  
  S_n_values = sapply(1:(n-1), function(x) round(S_n(J,p_sigma,x)))
  abundances = unlist(sapply(1:length(S_n_values), function(x) rep(x, S_n_values[x])))
  return(abundances)
} 

## Simulate detection/extinction events
run_sim = function(all_det_probs, all_ext_probs){
  
  det_events = apply(all_det_probs, c(1,2), function(x) rbinom(1,1,x))
  ext_events = apply(all_ext_probs, c(1,2), function(x) rbinom(1,1,x))
  
  det_year = numeric(nrow(all_det_probs))
  ext_year = numeric(nrow(all_ext_probs))
  for ( spp_index in 1:nrow(all_det_probs)){
    det_year[spp_index] = if ( sum(det_events[spp_index,] == 1) != 0){
      min(which(det_events[spp_index,]==1))
    } else {  
      sim_TT+1
    }
    
    ext_year[spp_index] = if ( sum(ext_events[spp_index,] == 1) != 0){
      min(which(ext_events[spp_index,]==1))
    } else { 
      sim_TT+1 
    }
  }
  spp_years = data.frame(det_year = det_year, ext_year = ext_year)
  return(spp_years)
}

## To run the SEUX model
get_SEUX_est = function(my_df, CI = FALSE){
  
  ## Adapted from Chisholm et al. (2016)
  # Compute p
  compute_p_t <- function(S_t,E_t){
    stopifnot(length(S_t)==length(E_t))
    # Method assumes that initially zero species have been detected
    stopifnot(S_t[1]==0 & E_t[1]==0)
    TT = length(S_t)-1
    mu_t_est = diff(E_t)/S_t[-(TT+1)]
    # By assumption, the extinction rate in the initial year is zero
    mu_t_est[1] = 0
    p_t = c(0,1-cumprod(1-mu_t_est))
    return(p_t)
  }
  # Compute Xt
  compute_X_t <- function(S_t,E_t){
    p_t = compute_p_t(S_t,E_t)
    TT = length(S_t)-1
    X_t_est = S_t[TT+1]*p_t/(1-p_t[TT+1])-E_t
    return(X_t_est)
  }
  # Compute Xt*
  compute_X_t_star <- function(S_t,E_t,mu_t_star){
    stopifnot(length(S_t)==length(E_t) & length(mu_t_star)==length(S_t)-1)
    
    TT0 = length(S_t)-1
    Y0 = (S_t[-1]+E_t[-1]-S_t[-(TT+1)]-E_t[-(TT+1)])/(E_t[TT+1]+S_t[TT+1])
    
    iii = which(S_t+E_t==max(S_t+E_t))
    # This assumes that the number of undetected species reached zero after the
    # last detection. This is necessary, because it's impossible to estimate
    # the rate of detection if no species are being detected (and so we
    # effectively assume the rate is zero after this time).
    TT = iii[1]-1
    
    Y = Y0[1:TT]
    mu_t_star = c(mu_t_star,0)
    
    temp = numeric(TT)
    
    for ( j in 0:(TT-2) )
    {
      temp[j+1] = prod(1-mu_t_star[((j+1):(TT-1))+1])
    }
    # Special case for j = TT-1 (in the product above, the start index would be
    # greater than the end index for this case, which means the product should be
    # equal to 1)
    temp[(TT-1)+1] = 1
    
    nu_est = numeric(TT+1)
    for ( i in 0:(TT-1) )
    {
      nu_est[i+1] = Y[i+1]*prod(1-mu_t_star[(i:(TT-1))+1])/sum((Y*temp)[(i:(TT-1))+1])
    }
    
    temp2 = numeric(TT-1+1)
    for ( j in 0:(TT-1) )
    {
      if ( j==0 )
        temp2[j+1] = mu_t_star[j+1]
      else
        temp2[j+1] = mu_t_star[j+1]*prod((1-mu_t_star[1:TT]-nu_est[1:TT])[(0:(j-1))+1])
    }
    
    X_t_est = (E_t[TT+1]+S_t[TT+1])*c(0,cumsum(temp2))/(1-sum(temp2))
    if ( TT!=TT0 ) X_t_est = c(X_t_est,rep(tail(X_t_est,1),TT0-TT))
    return(X_t_est)
  }
  
  # Preparing time series data
  names(my_df) = c("det_years", "ext_years")
  stopifnot(all(my_df$det_years <= my_df$ext_years))
  year_min = min(my_df$det_years)
  year_max = max(my_df$ext_years)
  years = year_min:year_max
  
  TT = length(years)
  S_t = numeric(TT+1)
  E_t = numeric(TT+1)
  
  # Compute the cumulative number of detected extant and detected extinct species in each year
  for ( t in 1:TT ){
    year = years[t]
    S_t[t+1] = sum(my_df$det_years<=year & my_df$ext_years>=year)
    E_t[t+1] = sum(my_df$ext_years<year)
  }
  
  # Estimate the number of undetected extinctions
  X_t_est = compute_X_t(S_t,E_t)
  
  # Estimate number of undetected extant species
  NN = S_t[TT+1]+E_t[TT+1]+X_t_est[TT+1]
  U_t_est = NN-S_t-E_t-X_t_est
  U_t_est_rounded = round(U_t_est) 
  
  result = list(UE = X_t_est, 
                S_t = S_t, 
                E_t = E_t,
                pE = (X_t_est+E_t)/(S_t+E_t+X_t_est))
  
  if (CI == TRUE){
    # Bootstrap to get CI on the number of undetected extinct species
    n_rep = 1000
    all_X = matrix(NA,n_rep,length(X_t_est))
    
    mu_t_est = diff(E_t)/S_t[-(TT+1)]
    mu_t_est[1] = 0
    
    for ( j in 1:n_rep ){
      # This randomisation procedure generates new mu values from the binomial
      # distribution, based on the observed mu and S, and estimated U.
      # A resampled set of mu values is discarded if one of the mu values is 1, since this causes 
      # the corresponding estimate of X to be undefined.
      success_flag = 0
      while(success_flag == 0){
        mu_resample = rbinom(TT,S_t[-(TT+1)],mu_t_est)/S_t[-(TT+1)]
        mu_resample[1] = 0
        mu_star_resample = rbinom(TT,U_t_est_rounded[-(TT+1)],mu_resample)/U_t_est_rounded[-(TT+1)]
        U_t_est_rounded_indzero = which(U_t_est_rounded[-(TT+1)]==0)
        mu_star_resample[U_t_est_rounded_indzero] = 0
        if(max(mu_star_resample)<1){
          XX = compute_X_t_star(S_t,E_t,mu_star_resample)
          success_flag = 1
        }
      }
      all_X[j,] = XX
    }
    X_lo = apply(all_X,2,function(x) quantile(x,0.025))
    X_hi = apply(all_X,2,function(x) quantile(x,0.975))
    
    # pE_lo = round((E_t+X_lo)/(S_t+E_t+X_lo), 3)
    # pE_hi = round((E_t+X_hi)/(S_t+E_t+X_hi), 3)
    
    pE_lo = (E_t+X_lo)/(S_t+E_t+X_lo)
    pE_hi = (E_t+X_hi)/(S_t+E_t+X_hi)
    
    CIs = list(UE_lo = X_lo,
               UE_hi = X_hi,
               pE_lo = pE_lo,
               pE_hi = pE_hi,
               all_X = all_X)
    
    result = c(result, CIs)
  }
  return(result)
}

## To run the Tedesco model
get_Tedesco_est = function(my_df, last_ext_year = NA){
  
  names(my_df) = c("det_year", "ext_year")
  
  ## Step 1: Getting extinction rate, mu
  max_T = max(my_df$ext_year)
  
  x_t = numeric(max_T)
  for ( i in 0:(max_T-1) ){
    x_t[i+1] = sum(my_df$det_year == i)
  }
  NE = sum(my_df$ext_year != max_T)
  
  if(is.na(last_ext_year)){last_ext_year = max_T}
  
  eq13 = function(mu){
    
    # Probability of non-persistence over all years
    i = seq(last_ext_year, 1, -1)
    p = 1-(1-mu)^i
    
    # Species are expected to be observed extinct if they do not persist till the tt-th year
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
    nll = -sum(spp*(log(ED))-ED)
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


## Simulations ==========================================================

# Parameters
lambda = omega = 1             # Birth and death rates (for the neutral BD process)
J = 1e6                        # Community size  (for log-series, analagous to J in a neutral model)
p_sigma = 5e-5                 # Speciation rate (for log-series, analagous to nu in a neutral model)
# p_sigma = 1e-6               # Speciation rate for simulations with higher abundances
sim_TT = 50                    # Total time-steps in the simulation
DeltaT = 1                     # Length of each time-step
xi = 0.1                       # Per-capita detection rate
k_values = seq(0, 1, 0.05)     # k determines the rate at which the per-capita detection rate changes with abundance

spp_abundances = get_spp_abundances(J, p_sigma)
cat("Total number of species: ", length(spp_abundances), "\n",
    "Highest abundance: ", max(spp_abundances), "\n", sep = "")

k = k_values[1] # For a single k_value

# for (k in k_values){   # Uncomment to iterate over all k_values

# Getting probabilities (long run-time, but only needs to be run once for a given parameter set)
det_probs = ext_probs = matrix(nrow = length(unique(spp_abundances)), ncol = sim_TT+1)
count = 0
for ( spp in unique(spp_abundances) ){
  count = count+1
  ext_probs[count,] = c(spp, sapply(1:sim_TT, function(tt){
    get_simulation_mus(lambda, p_sigma, spp, tt, DeltaT)}))
  det_probs[count,] = c(spp, sapply(1:sim_TT, function(tt){
    get_simulation_nus(lambda, omega, xi, p_sigma, spp, tt, DeltaT, k = k)}))
  if(count%%20 == 0) { print(count) }
}

# Assigning probabilities
all_ext_probs = all_det_probs = matrix(nrow = length((spp_abundances)), ncol = sim_TT)
for (abun_no in unique(spp_abundances)){
  to_update = which(spp_abundances==abun_no)
  for ( ii in to_update ){
    all_ext_probs[ii,] = ext_probs[which(ext_probs[,1] == abun_no),-1]
    all_det_probs[ii,] = det_probs[which(det_probs[,1] == abun_no),-1]
  }
}

## Running simulation =======================================================
n_sims = 1000         # Number of simulations

true_UE = SEUX_est = Tedesco_est = numeric(n_sims)
true_pE = SEUX_pE = Tedesco_pE = numeric(n_sims)
U_T = E_T = Tedesco_N = numeric(n_sims)
for(sim in 1:n_sims){
  
  # Generate data set
  sim_results = run_sim(all_det_probs[,1:sim_TT], all_ext_probs[,1:sim_TT])
  observed_data = sim_results[sim_results$det_year<=sim_results$ext_year &
                                sim_results$det_year != sim_TT+1,]
  true_UE[sim] = sum(sim_results$det_year>sim_results$ext_year & sim_results$ext_year != sim_TT+1)
  true_pE[sim] = sum(sim_results$ext_year != sim_TT+1)/nrow(sim_results)
  
  # Run models
  SEUX_model = get_SEUX_est(observed_data)
  Tedesco_model = get_Tedesco_est(observed_data)
  
  SEUX_est[sim] = tail(SEUX_model$UE,1)
  Tedesco_est[sim] = tail(Tedesco_model$UE,1)
  SEUX_pE[sim] = tail(SEUX_model$pE,1)
  Tedesco_pE[sim] = Tedesco_model$pE
  
  Tedesco_N[sim] = Tedesco_model$N0
  U_T[sim] = sum(sim_results$det_year == sim_TT+1 & sim_results$ext_year == sim_TT+1)
  E_T[sim] = sum(observed_data$ext_year != sim_TT+1 & observed_data$det_year<=observed_data$ext_year)
 
  if(sim%%100==0){print(sim)} 
}

cat("Average undetected extinctions", "\n",
    " Actual: ", mean(true_UE), "\n",
    "   SEUX: ", round(mean(SEUX_est),2), 
    " (", round(mean((SEUX_est-true_UE)/true_UE),3)*100, "%)","\n",
    "Tedesco: ", round(mean(Tedesco_est),2),
    " (", round(mean((Tedesco_est-true_UE)/true_UE),3)*100, "%)", "\n", sep = "")

cat("Average proportional extinctions", "\n",
    " Actual: ", mean(true_pE), "\n",
    "   SEUX: ", round(mean(SEUX_pE),2), 
    " (", round(mean((SEUX_pE-true_pE)/true_pE),3)*100, "%)","\n",
    "Tedesco: ", round(mean(Tedesco_pE),2),
    " (", round(mean((Tedesco_pE-true_pE)/true_pE),3)*100, "%)", "\n", sep = "")

res_df = data.frame(true_UE, SEUX_est, Tedesco_est,
                    true_pE, SEUX_pE, Tedesco_pE,
                    Tedesco_N, U_T, E_T)

#}  # Uncomment to iterate over all k_values


## Real-world data ==================================================
file_names = Sys.glob("./S4*.csv") # with datasets in working directory

# Start and end years as in Tedesco et al. (2014)
start_years = c(1979, 1805, 1840, 1825, 1855, 1842, 1894, 1758)
end_years = c(2007, 2010, 2008, 2010, 2009, 2010, 2010, 2005)

all_results = data.frame(matrix(nrow = 8, ncol = 25))
names(all_results) = c("dataset", 
                       "Tedesco_mu", "Tedesco_mu_lo", "Tedesco_mu_hi",
                       "Tedesco_nu", "Tedesco_nu_lo", "Tedesco_nu_hi",
                       "Tedesco_N0", "Tedesco_N0_lo", "Tedesco_N0_hi",
                       "Tedesco_UE", "Tedesco_UE_lo", "Tedesco_UE_hi",
                       "Tedesco_pE", "Tedesco_pE_lo", "Tedesco_pE_hi",
                       "SEUX_UE", "SEUX_UE_lo", "SEUX_UE_hi", 
                       "SEUX_pE", "SEUX_pE_lo", "SEUX_pE_hi",
                       "Hybrid_UE", "Hybrid_UE_lo", "Hybrid_UE_hi")

group_names = substr(file_names, 6, nchar(file_names)-4)
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
  
  # Running the models
  s_result = get_SEUX_est(my_df2[2:3], CI = TRUE)
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
  
  # Getting results from hybrid model
  hybrid_UE = tail(s_result$UE+((t_result$N0-nrow(my_df))*s_result$pE)/(1-s_result$pE),1)
  hybrid_UE_lo = tail(s_result$UE_lo+((t_result$N0-nrow(my_df))*s_result$pE_lo)/(1-s_result$pE_lo),1)
  hybrid_UE_hi = tail(s_result$UE_hi+((t_result$N0-nrow(my_df))*s_result$pE_hi)/(1-s_result$pE_hi),1)
  
  all_results[i,-1] = c(round(c(t_result$mu, quantile(temp_mu, c(0.025, 0.975))),7),
                        round(c(t_result$nu, quantile(temp_nu, c(0.025, 0.975))),5),
                        round(c(t_result$N0, quantile(temp_N0, c(0.025, 0.975))),0),
                        round(c(tail(t_result$UE,1), quantile(temp_UE, c(0.025,0.975))),1),
                        round(c(t_result$pE, quantile(temp_pE, c(0.025, 0.975))),3),
                        round(c(tail(s_result$UE,1), tail(s_result$UE_lo,1), tail(s_result$UE_hi,1)),1),
                        round(c(tail(s_result$pE,1), tail(s_result$pE_lo,1), tail(s_result$pE_hi,1)),3),
                        round(c(hybrid_UE, hybrid_UE_lo, hybrid_UE_hi),1))

  cat("dataset", i, "of", length(file_names), "\n")
}

all_results

## Logistic regression ==============================
lr_results = data.frame(dataset = group_names,
                        delta_log_odds = numeric(8),
                        SE_log_odds = numeric(8),
                        deviance_val = numeric(8),
                        p_val = numeric(8))
for (i in 1:8){
  my_df = read.csv(file_names[i])
  
  extinct = ifelse(is.na(my_df$extinction_year), 0, 1)
  temp_data = data.frame(d_time = my_df$description_year - min(my_df$description_year),
                         extinct = extinct,
                         group = group_names[i])
  
  # Logistic regression model
  temp_lm = glm(temp_data$extinct ~ temp_data$d_time, family = "binomial")
  temp_lrt = anova(temp_lm, test = "LRT")
  
  lr_results$delta_log_odds[i] = round(summary(temp_lm)$coefficients[2,1], 5)
  lr_results$SE_log_odds[i] = round(summary(temp_lm)$coefficients[2,2], 5)
  lr_results$deviance_val[i] = round(temp_lrt$Deviance[2], 3)
  lr_results$p_val[i] = round(temp_lrt$`Pr(>Chi)`[2], 3)
  
}

lr_results
