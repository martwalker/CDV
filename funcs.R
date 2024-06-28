#####################################
## functions
####################################
funcs <- list(
  
  #######################
  # stochastic SEIRV model simulation 
  #######################   
  SEIRV = function(states, par) {
    
    #######################
    # initialise and populate compartments
    #######################
    S <- E <- I <- R <- V <- N <- Ex <- vector("integer", length=par$n_patches)
    
    for (i in 1:par$n_patches) {
      S[i] <- (states$N0[i] - states$E0[i])*(1-par$vacc[i])
      E[i] <- states$E0[i]
      I[i] <- 0
      R[i] <- 0
      V[i] <- (states$N0[i] - states$E0[i])*par$vacc[i]
      N[i] <- S[i] + E[i] + I[i] + R[i] + V[i]
      Ex[i] <- N[i]==0
    }
    
    # recovery rate
    gamma <- 1/par$dur_infectious
    
    # proggression rate 
    sigma <- 1/par$dur_latent
    
    # mortality rate
    mu <- par$mort*gamma/(1-par$mort)
    
    #  contact matrix
    beta <- par$R0mat*(gamma+mu)
    
    # rate of losing immunity
    delta <- 1/par$dur_immun
    
    # initialise storage matrices 
    S_series <- E_series <- I_series <-
      R_series <- V_series <- N_series <- Ex_series <- matrix(0, nrow = par$T, ncol = par$n_patches)    
    
    #######################
    # main time loop
    #######################     
    for (t in 1:par$T) {
      S_series[t, ] <- S
      E_series[t, ] <- E
      I_series[t, ] <- I
      R_series[t, ] <- R
      V_series[t, ] <- V
      N_series[t, ] <- S + E + I + R + V
      Ex_series[t, ] <- Ex
      
      new_infection <- new_endogenous_infection <-  new_exogenous_infection <- 
        new_progression <- new_loss <- new_recovered <-
          new_dead <- new_vaccinated <- new_lost_immunity <- rep(0, par$n_patches)
      
      for (patch in 1:par$n_patches) {
        
        #######################
        ## Event probabilities
        #######################
        
        if (!is.na(N[patch]) && N[patch] > 0) {
          
          p_infection <- 1 - exp(-(
            (beta[patch,patch] * I[patch] / N[patch]) +                # endogenous infection
              (sum(beta[patch,-patch]*I[-patch] / N[-patch], na.rm=TRUE))  # exogenous infection
          ))
          
          # endogenous infection event
          p_endogenous_infection <- 1 - exp(-(
            (beta[patch,patch] * I[patch] / N[patch])
          ))
          
          # exogenous infection event
          p_exogenous_infection <- 1 - exp(-(
            (sum(beta[patch,-patch]*I[-patch] / N[-patch])) 
          ))
          
          # progression from exposed to infectious
          p_progression <- 1 - exp(-(
            sigma
          ))
          
          # any loss event (recovery or death)
          p_loss <- 1 - exp(-(
            gamma +                                            # recovery 
              mu                                               # death
          ))
          
          # recovery event
          p_recovery <- 1 - exp(-(
            gamma
          ))
          
          # death event
          p_death <- 1 - exp(-(
            mu
          ))
          
          # losing immunity event
          p_lost_immunity <- 1 - exp(-(
            delta
          ))
          
        } else {
          
          p_infection <- p_progression <- p_endogenous_infection <- p_exogenous_infection <-  p_loss <- p_recovery <-  p_death <- p_vaccination <- p_lost_immunity <- 0
          
        }
        
        #######################
        ## numbers of events
        #######################
        
        
        if (p_infection>0) {
          # new infection (endogenous or exogenous)
          new_infection[patch] <- rbinom(1, S[patch], p_infection)
          # endogenous infections
          new_endogenous_infection[patch] <- rbinom(1, new_infection[patch], p_endogenous_infection/p_infection)
          # exogenous infections
          new_exogenous_infection[patch] <- new_infection[patch] - new_endogenous_infection[patch]
        } else {
          new_infection[patch] <- new_endogenous_infection[patch] <- new_exogenous_infection[patch] <- 0
        }
        
        if (p_progression) {
          new_progression[patch] <- rbinom(1, E[patch], p_progression)
        }
        
        if (p_loss>0) {
          # new loss (recovery or death)
          new_loss[patch] <- rbinom(1, I[patch], p_loss)
          # new recovery
          new_recovered[patch] <- rbinom(1, new_loss[patch], p_recovery/p_loss)
          # new death 
          new_dead[patch] <- new_loss[patch] - new_recovered[patch]
        } else {
          new_loss[patch] <- new_recovered[patch] <- new_dead[patch] <- 0
        }
        
        if (p_lost_immunity > 0) {
          new_lost_immunity[patch] <- rbinom(1,V[patch], p_lost_immunity)
        } else {
          new_lost_immunity[patch] <- 0
        }
        
      }
      
      #######################
      ## updates states
      #######################
      
      S <- S - new_infection - new_vaccinated + new_lost_immunity
      E <- E + new_infection - new_progression
      I <- I + new_progression - new_recovered - new_dead
      R <- R + new_recovered
      V <- V + new_vaccinated - new_lost_immunity
      N <- S + E + I + R + V
      Ex <- N==0
    }
    
    #######################
    ## output time series
    #######################
    
    
    list(S=S_series,E=E_series, I=I_series, R=R_series, V=V_series, N=N_series, Ex=Ex_series)
    
  }, 
  
  runSEIRV = function(repeats, states, par) {
    
    multi <- replicate(repeats, funcs$SEIRV(states = states, par = par))
    
    S <- as.data.frame(cbind(do.call(rbind, multi["S", 1:repeats]), time = seq(1, par$T)))
    E <- as.data.frame(cbind(do.call(rbind, multi["E", 1:repeats]), time = seq(1, par$T)))
    I <- as.data.frame(cbind(do.call(rbind, multi["I", 1:repeats]), time = seq(1, par$T)))
    R <- as.data.frame(cbind(do.call(rbind, multi["R", 1:repeats]), time = seq(1, par$T)))
    V <- as.data.frame(cbind(do.call(rbind, multi["V", 1:repeats]), time = seq(1, par$T)))
    N <- as.data.frame(cbind(do.call(rbind, multi["N", 1:repeats]), time = seq(1, par$T)))
    Ex <- as.data.frame(cbind(do.call(rbind, multi["Ex", 1:repeats]), time = seq(1, par$T)))
    
    
    S_df <- as.data.frame(S) %>% pivot_longer(-time, names_to = "patch", values_to = "S")
    E_df <- as.data.frame(E) %>% pivot_longer(-time, names_to = "patch", values_to = "E")
    I_df <- as.data.frame(I) %>% pivot_longer(-time, names_to = "patch", values_to = "I")
    R_df <- as.data.frame(R) %>% pivot_longer(-time, names_to = "patch", values_to = "R")
    V_df <- as.data.frame(V) %>% pivot_longer(-time, names_to = "patch", values_to = "V")
    N_df <- as.data.frame(N) %>% pivot_longer(-time, names_to = "patch", values_to = "N")
    Ex_df <- as.data.frame(Ex) %>% pivot_longer(-time, names_to = "patch", values_to = "Ex")
    
    
    S_out <- group_by(S_df, patch, time) %>% summarise(mean = mean(S), lwr = quantile(S, probs = 0.025, na.rm=TRUE), upr = quantile(S, probs = 0.975, na.rm=TRUE))
    E_out <- group_by(E_df, patch, time) %>% summarise(mean = mean(E), lwr = quantile(E, probs = 0.025, na.rm=TRUE), upr = quantile(E, probs = 0.975, na.rm=TRUE))
    I_out <- group_by(I_df, patch, time) %>% summarise(mean = mean(I), lwr = quantile(I, probs = 0.025, na.rm=TRUE), upr = quantile(I, probs = 0.975, na.rm=TRUE))
    R_out <- group_by(R_df, patch, time) %>% summarise(mean = mean(R), lwr = quantile(R, probs = 0.025, na.rm=TRUE), upr = quantile(R, probs = 0.975, na.rm=TRUE))
    V_out <- group_by(V_df, patch, time) %>% summarise(mean = mean(V), lwr = quantile(V, probs = 0.025, na.rm=TRUE), upr = quantile(V, probs = 0.975, na.rm=TRUE))
    N_out <- group_by(N_df, patch, time) %>% summarise(mean = mean(N), lwr = quantile(N, probs = 0.025, na.rm=TRUE), upr = quantile(N, probs = 0.975, na.rm=TRUE))
    Ex_out <- group_by(Ex_df, patch, time) %>% summarise(mean = mean(Ex), lwr = quantile(Ex, probs = 0.025, na.rm=TRUE), upr = quantile(Ex, probs = 0.975, na.rm=TRUE))
    
    
    return(list(S = S_out, E = E_out, I = I_out,R = R_out, V = V_out, N = N_out, Ex=Ex_out))
  }
)
