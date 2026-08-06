# =============================================================================
# Title:       Euler / Venn Diagram
# Description: Draws a proportional Euler diagram for three sets (EZE25, ATO5,
#              EZE25_ATO5) using the eulerr package. Regions unique to or
#              shared with EZE25_ATO5 are filled with the brand colour
#              #00A2BC; all other regions are left white.
# Inputs:      None (values are hard-coded).
# Outputs:     Euler diagram rendered to the active graphics device.
# Author:      Federica D'Annunzio
# =============================================================================

library(eulerr)

# Define set sizes and pairwise / triple overlap sizes
euler_plot <- euler(c(
  "A"     = 2,    "B"     = 2,    "C"     = 2,
  "A&B"   = 0.5,  "A&C"   = 0.5,  "B&C"   = 0.5,
  "A&B&C" = 0.25
))

# Plot the Euler diagram with custom fills and set labels
plot(
  euler_plot,
  quantities = FALSE,
  fills = c(
    "A"     = "white",   "B"     = "white",   "C"     = "#00A2BC",
    "A&B"   = "white",   "A&C"   = "#00A2BC", "B&C"   = "#00A2BC",
    "A&B&C" = "#00A2BC"
  ),
  labels = c("EZE25", "ATO5", "EZE25_ATO5")
)
