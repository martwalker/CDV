# Install and load the DiagrammeR package if not already installed
install.packages("DiagrammeR")
library(DiagrammeR)

# Create and render the diagram
grViz("
digraph SEIRV {
    node [shape=box style=rounded]
    
    S [label='Susceptible (S)']
    E [label='Exposed (E)']
    I [label='Infectious (I)']
    R [label='Recovered (R)']
    V [label='Vaccinated (V)']
    D [label='Death (D)']
    
    S -> E [label='Infection']
    E -> I [label='Progression']
    I -> R [label='Recovery']
    I -> D [label='Death']
    S -> V [label='Vaccination']
    V -> S [label='Loss of Immunity']
}
")

# Render the diagram
render_graph(grViz("
digraph SEIRV {
    node [shape=box style=rounded]
    
    S [label='Susceptible (S)']
    E [label='Exposed (E)']
    I [label='Infectious (I)']
    R [label='Recovered (R)']
    V [label='Vaccinated (V)']
    D [label='Death (D)']
    
    S -> E [label='Infection']
    E -> I [label='Progression']
    I -> R [label='Recovery']
    I -> D [label='Death']
    S -> V [label='Vaccination']
    V -> S [label='Loss of Immunity']
}
"))
