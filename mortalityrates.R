######################################
## data from Gorden et al. 2015
####################################
df <- read.csv("data.csv")

#####################################
## mortality rates (subadults and adults)
####################################
n0 <- df$total[df$time==0]
n1 <- df$total[df$time==183]
n2 <- df$total[df$time==366]

1-n2/n0
