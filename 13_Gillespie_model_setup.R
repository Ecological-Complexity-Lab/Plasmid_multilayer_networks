####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 13: Set up Gillespie stochastic simulations code for gene dispersion in a network 
# 
# 
# Script tested for R version 4.1.1
####################################################################################################################

####################################################################################################################
# SCRIPT SET-UP
####################################################################################################################
# Set working directory to wherever your files are located

# Load the necessary packages:
library(tidyverse)
library(igraph)

# Starting files:
# Edgelist with edgeweights created in script 03_Network_setup.R
load("net.dat2k.ew.Rda")

# Degree of physical nodes, created in script 04_Basic_network_statistics.R
load("deg.str.2k.all.phys.Rda")

# Infomap data created in script 08_Infomap_analysis_full_network.R
load("plas_mods.df.Rda")

# Statistics on infomap modules created in script 08_Infomap_analysis_full_network.R
load("plas_mods.stats.Rda")
####################################################################################################################


####################################################################################################################
# Section 1: CREATE GILLESPIE FUNCTION
####################################################################################################################
# Gillespie function: 
gillespie <- function(boolean, # Your state nodes with or without gene
                       sim.matrix, # Matrix defining the flow
                       pars, # Parameters of the model
                       time_v){ #And a vector of times when things happen
  # Rates
  p <- pars$success_probability # For gene transfer
  cont <- pars$contact_rate # Rate at which plasmids contact (NOTE that I see this as a total rate, not a per capita one)
  loss <- pars$loss_rate # Gene loss per capita rate (NOTE that this changes when we gain or loose genes)

  # States: note that the state of the system is encoded in the boolean but it
  # will be handy to track the number of genes in the (meta)population.
  G <- sum(boolean)
  
  n_plasmids <- length(boolean) # Needed for later
  
    # Propensities
  r1 <- function(){cont}  # Two plasmids collide
  r2 <- function(){loss * G} # Losing a gene

  # Vector of propensities
  rates <- c(r1(), r2())
  
  # Names of events
  events <- c("f1", "f2")
  # Events: each one actualizes the corresponding states and propensities.
  f1 <- function(){ #Two plasmids collide
    x1 <- sample.int(n_plasmids, 1)
    x2 <- sample.int(n_plasmids, 1)
    
    if(boolean[x1] != boolean[x2]){
      gain <- runif(1) < (p * sim.matrix[x1, x2])
      if (gain) boolean[x1] <<- T; boolean[x2] <<- T; G <<- G + 1; rates[2] <<- r2() 
    }
  }
  f2 <- function(){ #A gene is loss at random
    index <- which(boolean)
    lost <- sample.int(G, 1)
    boolean[index[lost]] <<- F
    G <<- G - 1
    rates[2] <<- r2()
    rm(index); rm(lost)
    }

  
  # Initialize variables
  t <- time_v[1]
  out <- matrix(nrow = n_plasmids, ncol = length(time_v))
  out[, 1] <- boolean
  
  # Note: a little fraction of time is out of sync. Correct when everything else is settled.
  # Size of error depends on contact rate (higher contact rate = lower error)
  
  # Start simulations
  for(i in 2:length(time_v)){
    while (t < time_v[i]){
      
      # Sampling
      dt <- stats::rexp(1, sum(rates))
      
      z <- sample.int(length(rates), 1, prob = rates)
      do.call(events[z], args = list())
      
      t <- t + dt
    }
    out[, i] <- boolean #saving the states
  }
  
  out
}
####################################################################################################################


####################################################################################################################
# Section 2: CREATE A MATRIX FROM THE EDGELIST FOR THE GILLESPIE MODEL
####################################################################################################################

pre.mat1 <- net.dat2k.ew %>%
  rowwise() %>%
  mutate(fr_node = paste("N",node_from,sep = "")) %>%
  mutate(fr_lay = paste("L",layer_from,sep = "")) %>%
  mutate(fr_grp = paste(fr_node, fr_lay, sep = "_")) %>%
  mutate(to_node = paste("N",node_to,sep = "")) %>%
  mutate(to_lay = paste("L",layer_to,sep = "")) %>%
  mutate(to_grp = paste(to_node, to_lay,  sep = "_"))

pre.mat2 <- pre.mat1 %>%
  ungroup() %>%
  select(fr_grp, to_grp, weight)

pre.st.fr <- pre.mat1 %>%
  ungroup() %>%
  select(layer_from, node_from, fr_grp) %>%
  distinct() %>%
  rename(layer_id=layer_from,
         node_id=node_from,
         st.node.id=fr_grp)

graph1 <- graph_from_data_frame(pre.mat2, directed=F)

mat <- get.adjacency(graph1, sparse = FALSE, attr='weight', type="both")

# Save the matrix:
save(mat, file="mat.Rda")

# Bind ordered node list to state node characteristics
st.node.list <- pre.mat1 %>%
  ungroup() %>%
  select(layer_to, node_to, to_grp) %>%
  distinct() %>%
  rename(layer_id=layer_to,
         node_id=node_to,
         st.node.id=to_grp) %>%
  bind_rows(., pre.st.fr) %>%
  left_join(., deg.str.2k.all.phys, by="node_id") %>%
  select(-strength, -edge_type) %>%
  left_join(., plas_mods.df) %>%
  left_join(., plas_mods.stats)

st.node.order <- as.data.frame(mat) %>%
  rownames_to_column(., var="st.node.id") %>%
  select(st.node.id) %>%
  # Give a rank 
  rownames_to_column(., var="mat.order") %>%
  mutate(mat.order=as.numeric(mat.order))

st.node.list.ord <- st.node.list %>%
  left_join(., st.node.order) %>%
  mutate(layer_id=as.factor(layer_id)) %>%
  arrange(mat.order) %>%
  distinct()

# Save the data
save(st.node.list.ord, file="st.node.list.ord.Rda")

mat.ids <- st.node.list.ord %>%
  select(mat.order, layer_id)

# Save the data
save(mat.ids, file="mat.ids.Rda")


