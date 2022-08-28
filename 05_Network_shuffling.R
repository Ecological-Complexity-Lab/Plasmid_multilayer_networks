####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 5: Randomize the observed network by shuffling the the layers in the extended edge-list
# 
# The following outputs are used for downstream analysis: 
# shuff.lf.list.Rda
# 
# Script tested for R version 4.1.1
####################################################################################################################

####################################################################################################################
# SCRIPT SET-UP
####################################################################################################################

# Set working directory to wherever your files are located

# Load necessary packages:
library(tidyverse)

# Starting files:
# Edgelist with edgeweights created in script 03_Network_setup.R
load("net.dat2k.ew.Rda")
####################################################################################################################


####################################################################################################################
# Section 1: SHUFFLE THE NETWORK
####################################################################################################################
# Select the columns from the edgelist that will be shuffled:
pre.shuff <- net.dat2k.ew %>%
  ungroup() %>%
  select(layer_from, node_from, layer_to, node_to, weight)

# Let's shuffle
set.seed(123)

# Function to shuffle the layer_from
shuff.lf.fun <- function(x) {
  transform(pre.shuff, layer_from = sample(layer_from, replace = F) )
}

# Choose how many times to shuffle the data using the function created above.
# Creates a list with each shuffling
shuff.lf.list <- (lapply(seq(1000),shuff.lf.fun))

# Give each shuffled list a name corresponding to the number of its shuffling
names(shuff.lf.list) <- seq(1000)

save(shuff.lf.list, file="shuff.lf.list.Rda")
####################################################################################################################