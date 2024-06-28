#####################################
## parameters
####################################
par <- list(n_patches=7,    
            T=100, 
            dur_infectious=13, ## duration of infectiousness (days)
            dur_latent=7,      ## duration of latent period (days)
            mort=0.68,         ## proportion dying  
            R0intra=5)         ## intra-patch R0

par$R0mat <- matrix(c( par$R0intra, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                       0.5, par$R0intra, 0.5, 0.5, 0.5, 0.5, 0.5, 
                       0.5, 0.5, par$R0intra, 0.5, 0.5, 0.5, 0.5,
                       0.5, 0.5, 0.5, par$R0intra, 0.5, 0.5, 0.5, 
                       0.5, 0.5, 0.5, 0.5, par$R0intra, 0.5, 0.5,
                       0.5, 0.5, 0.5, 0.5, 0.5, par$R0intra, 0.5,
                       0.5, 0.5, 0.5, 0.5, 0.5, 0.5, par$R0intra),
                    nrow=par$n_patches, 
                    ncol=par$n_patches,byrow=T)

######################################
## data from Gorden et al. 2015
####################################
df <- read.csv("data.csv")

#####################################
## mortality rates
####################################
n0 <- df$total[df$time==0]
n1 <- df$total[df$time==183]
n2 <- df$total[df$time==366]

#####################################
## states
####################################
states = list(N0 = df$total[df$time==0], 
              E0 = c(0, 0, 0, 0, 1,  0, 0))
