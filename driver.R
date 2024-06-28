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
out <- funcs$runSEIRV(repeats = 100, states = states, par = par)

plot_results(out)

