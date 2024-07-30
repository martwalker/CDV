
##################################################
## Stochastic SEIRV metapopulation model simulation
#################################################

## Clear workspace
rm(list=ls())

## Load required packages
source("packages.R")

## Load default model parameters
source("par.R")

## Load functions for running model
source("funcs.R")

# Plotting function
plot_results <- function(out) {
  plot_susceptible <- ggplot(out$S, aes(x = time, y = mean, color = patch)) +
    geom_line() +
    labs(y = "Susceptible")
  
  plot_exposed <- ggplot(out$E, aes(x = time, y = mean, color = patch)) +
    geom_line() +
    labs(y = "Exposed")
  
  plot_infected <- ggplot(out$I, aes(x = time, y = mean, color = patch)) +
    geom_line() +
    labs(y = "Infected")
  
  plot_recovered <- ggplot(out$R, aes(x = time, y = mean, color = patch)) +
    geom_line() +
    labs(y = "Recovered")
  
  plot_population <- ggplot(out$N, aes(x = time, y = mean, color = patch)) +
    geom_line() +
    labs(y = "Population")
  
  plot_extinction <- ggplot(out$Ex, aes(x = time, y = mean, color = patch)) +
    geom_line() +
    labs(y = "Probability extinction")
  # Arrange plots in a grid
  grid.arrange(
    plot_susceptible, plot_exposed, plot_infected,
    plot_recovered, plot_population, plot_extinction,
    ncol = 2
  )
}

# Example usage:
out <- funcs$runSEIRV(repeats = par$repeats, states = states, par = par)

plot_results(out)

## Plot number patches going extict
#ggplot(out$NEx, aes(x=time, y=mean)) + 
 # geom_line() +
  #labs(y = "Number extinct")

out$Ex[,"mean"][out$Ex[,"time"]==par$T]

out$NEx[,"mean"][out$NEx[,"time"]==par$T]


##################################################
## run vaccination strategies
#################################################

## list packs by vulnerability
vulner <- order(out$Ex[,"mean"][out$Ex[,"time"]==par$T], decreasing=T)

## vaccinate most vulnerable packs
par$vacc[vulner[1]] <- 1
par$vacc[-vulner[1]] <- 0
out_v1 <- funcs$runSEIRV(repeats = par$repeats, states = states, par = par)
## vaccinate most and second most vulnerable packs
par$vacc[vulner[c(1,2)]] <- 1
par$vacc[-vulner[c(1,2)]] <- 0
out_v2 <- funcs$runSEIRV(repeats = par$repeats, states = states, par = par)
## vaccinate most, second and third most vulnerable packs
par$vacc[vulner[c(1,2,3)]] <- 1
par$vacc[-vulner[c(1,2,3)]] <- 0
out_v3 <- funcs$runSEIRV(repeats = par$repeats, states = states, par = par)
## vaccinate most, second, third and fourth most vulnerable packs
par$vacc[vulner[c(1,2,3,4)]] <- 1
par$vacc[-vulner[c(1,2,3,4)]] <- 0
out_v4 <- funcs$runSEIRV(repeats = par$repeats, states = states, par = par)
## vaccinate most, second, third, fourth and fifth most vulnerable packs
par$vacc[vulner[c(1,2,3,4,5)]] <- 1
par$vacc[-vulner[c(1,2,3,4,5)]] <- 0
out_v5 <- funcs$runSEIRV(repeats = par$repeats, states = states, par = par)
## vaccinate most, second, third, fourth, fifth and sixth most vulnerable packs
par$vacc[vulner[c(1,2,3,4,5,6)]] <- 1
par$vacc[-vulner[c(1,2,3,4,5,6)]] <- 0
out_v6 <- funcs$runSEIRV(repeats = par$repeats, states = states, par = par)
## vaccinate most, second, third, fourth, fifth, sixth and seventh most vulnerable packs
par$vacc[vulner[c(1,2,3,4,5,6,7)]] <- 1
par$vacc[-vulner[c(1,2,3,4,5,6,7)]] <- 0
out_v7 <- funcs$runSEIRV(repeats = par$repeats, states = states, par = par)

## bind results together
vdat <- as.data.frame( cbind(c(out$NEx[,"mean"][out$NEx[,"time"]==par$T],
          out_v1$NEx[,"mean"][out_v1$NEx[,"time"]==par$T],
           out_v2$NEx[,"mean"][out_v2$NEx[,"time"]==par$T],
            out_v3$NEx[,"mean"][out_v3$NEx[,"time"]==par$T],
              out_v4$NEx[,"mean"][out_v4$NEx[,"time"]==par$T],
                out_v5$NEx[,"mean"][out_v5$NEx[,"time"]==par$T],
                  out_v6$NEx[,"mean"][out_v6$NEx[,"time"]==par$T],
                    out_v7$NEx[,"mean"][out_v7$NEx[,"time"]==par$T]),
          c(out$NEx[,"lwr"][out$NEx[,"time"]==par$T],
              out_v1$NEx[,"lwr"][out_v1$NEx[,"time"]==par$T],
                out_v2$NEx[,"lwr"][out_v2$NEx[,"time"]==par$T],
                  out_v3$NEx[,"lwr"][out_v3$NEx[,"time"]==par$T],
                    out_v4$NEx[,"lwr"][out_v4$NEx[,"time"]==par$T],
                      out_v5$NEx[,"lwr"][out_v5$NEx[,"time"]==par$T],
                        out_v6$NEx[,"lwr"][out_v6$NEx[,"time"]==par$T],
                          out_v7$NEx[,"lwr"][out_v7$NEx[,"time"]==par$T]),
          c(out$NEx[,"upr"][out$NEx[,"time"]==par$T],
              out_v1$NEx[,"upr"][out_v1$NEx[,"time"]==par$T],
                out_v2$NEx[,"upr"][out_v2$NEx[,"time"]==par$T],
                  out_v3$NEx[,"upr"][out_v3$NEx[,"time"]==par$T],
                    out_v4$NEx[,"upr"][out_v4$NEx[,"time"]==par$T],
                      out_v5$NEx[,"upr"][out_v5$NEx[,"time"]==par$T],
                        out_v6$NEx[,"upr"][out_v6$NEx[,"time"]==par$T],
                          out_v7$NEx[,"upr"][out_v7$NEx[,"time"]==par$T]),
      c(0,1,2,3,4,5,6,7)))
colnames(vdat) <- c("mean", "lwr", "upr", "NEx")

ggplot(data=vdat, aes(y=mean,x=NEx)) +
    geom_point() +
      geom_line() +
        scale_y_continuous(name="mean number of packs extinct", limits=c(0, max(vdat$mean))) +
          scale_x_continuous(name="number of packs vaccinated")

