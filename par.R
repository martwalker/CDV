#####################################
## parameters
####################################
par <- list(n_patches = 7,    
            T = 100, 
            dur_infectious = 13, ## duration of infectiousness (days)
            dur_latent = 7,      ## duration of latent period (days)
            dur_immun = 547,     ## duration of immunity (days)
            mort = c(0.68, 0.68, 0.68, 0.68, 0.85, 0.68, 0.68),         ## proportion dying  
            beta_intra = 0.4,       ## intra-patch R0
            vacc = c(0, 1, 0, 0, 0, 0, 0))   ## indicator for whether a pack is vaccinated

par$beta <- matrix(c(par$beta_intra, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                      0.05, par$beta_intra, 0.05, 0.05, 0.05, 0.05, 0.05, 
                      0.05, 0.05, par$beta_intra, 0.05, 0.05, 0.05, 0.05,
                      0.05, 0.05, 0.05, par$beta_intra, 0.05, 0.05, 0.05, 
                      0.05, 0.05, 0.05, 0.05, par$beta_intra, 0.05, 0.05,
                      0.05, 0.05, 0.05, 0.05, 0.05, par$beta_intra, 0.05,
                      0.05, 0.05, 0.05, 0.05, 0.05, 0.05, par$beta_intra),
                    nrow = par$n_patches, 
                    ncol = par$n_patches, byrow = T)


#####################################
## states
####################################
states = list(N0 = c(2, 6, 4, 1, 10, 5, 3))