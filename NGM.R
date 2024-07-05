# Define parameters
n_patches <- 7
dur_infectious <- 13  # Duration of infectiousness (days)
dur_latent <- 7       # Duration of latent period (days)

# Empirical transition rates (example values)
gamma <- 1 / dur_infectious  # Recovery rate
sigma <- 1 / dur_latent      # Progression rate

# Contact rates (example values)
beta_intra <- 0.3  # Within-pack transmission rate
beta_inter <- 0.05 # Between-pack transmission rate

# Create the contact matrix (for beta values)
contact_matrix <- matrix(beta_inter, nrow = n_patches, ncol = n_patches)
diag(contact_matrix) <- beta_intra  # Setting the intra-pack transmission rate

# Create the F matrix (new infections produced)
F <- matrix(0, nrow = n_patches, ncol = n_patches)
for (i in 1:n_patches) {
  for (j in 1:n_patches) {
    F[i, j] <- contact_matrix[i, j] * sigma  # Use sigma as the progression rate for new infections
  }
}

# Create the V matrix (transition between states)
V <- diag(rep(gamma + sigma, n_patches))

# Calculate the next generation matrix K
K <- F %*% solve(V)

# Calculate the spectral radius of K (R0)
R0 <- max(Re(eigen(K)$values))  # Take the real part in case of numerical inaccuracies

# Print the basic reproduction number
print(R0)
