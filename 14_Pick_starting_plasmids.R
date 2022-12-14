####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 14: Pick starting plasmids for dynamical Gillespie model
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

# Starting files:
# Cow per plasmid created in script 04_Basic_network_statistics.R
load("plas.per.cow.Rda")

# Infomap data created in script 08_Infomap_analysis_full_network.R
load("plas_mods.df.Rda")

# Statistics on infomap modules created in script 08_Infomap_analysis_full_network.R
load("plas_mods.stats.Rda")

# Degree of physical nodes, created in script 04_Basic_network_statistics.R
load("deg.str.2k.all.phys.Rda")

# Ordered list of state nodes, created in script 13_Gillespie_function.R
load("st.node.list.ord.Rda")
####################################################################################################################


####################################################################################################################
# Section 1: RANDOMLY PICK CENTRAL PLASMIDS (BELONGING TO MODULE 1)
####################################################################################################################

# Identify state nodes in module 1 (the biggest module)
mod1 <- plas_mods.df %>%
  distinct() %>%
  filter(module==1) %>%
  left_join(., deg.str.2k.all.phys, by="node_id") %>%
  select(-strength, -edge_type) %>%
  mutate(layer_id=as.factor(layer_id)) %>%
  left_join(., plas_mods.stats) %>%
  left_join(., st.node.list.ord)

# Randomly sample state nodes in module 1
mod1.sampled <- mod1 %>%
  select(st.node.id, node_id) %>%
  #distinct(node_id) %>%
  group_by(node_id) %>%
  sample_n(., 1) %>%
  ungroup() %>%
  sample_n(., 10, replace=FALSE)  %>%
  left_join(., st.node.list.ord) %>%
  mutate(post = 1514 - mat.order) %>%
  arrange(mat.order) %>%
  mutate(pre = mat.order - 1,
         true = 1, 
         post = 1514 - mat.order)

# Save the data
save(mod1.sampled, file="mod1.sampled.Rda")
####################################################################################################################


####################################################################################################################
# Section 2: CREATE STARTING BOOLEAN LIST FOR CENTRAL PLASMIDS
####################################################################################################################
# Create an empty list
bool.centr.obs <- list()

# Create a boolean list for each starting plasmid, with the position of the starting plasmid TRUE and all other
# positions FALSE
for(i in 1:nrow(mod1.sampled)) {
  bt <- c(rep(FALSE, each=mod1.sampled$pre[i]), rep(TRUE,1), rep(FALSE,each=mod1.sampled$post[i]))
  
  bool.centr.obs[[length(bool.centr.obs)+1]]<-bt
  
}

# Assign name of each list as the starting plasmid
names(bool.centr.obs) <- c(paste("centr",mod1.sampled$st.node.id, sep = "_"))

# Save the data
save(bool.centr.obs, file="bool.centr.obs.Rda")
####################################################################################################################


####################################################################################################################
# Section 3: RANDOMLY PICK PERIPHERAL PLASMIDS (BELONGING TO SMALLEST MODULES)
####################################################################################################################

# Identify state nodes in the smallest modules (those with 2 state nodes)
not.mod1 <- plas_mods.df %>%
  distinct() %>%
  left_join(., deg.str.2k.all.phys, by="node_id") %>%
  select(-strength, -edge_type) %>%
  mutate(layer_id=as.factor(layer_id)) %>%
  left_join(., plas_mods.stats) %>%
  left_join(., st.node.list.ord) %>%
  filter(n.state.nodes==2)

# Join module and plasmid metadata to state nodes that are not in the biggest module (module 1)
not.mod1.cows <- plas_mods.df %>%
  distinct() %>%
  filter(module!=1) %>%
  left_join(., deg.str.2k.all.phys, by="node_id") %>%
  select(-strength, -edge_type) %>%
  mutate(layer_id=as.factor(layer_id)) %>%
  left_join(., plas_mods.stats) %>%
  left_join(., st.node.list.ord) %>%
  group_by(layer_id) %>%
  mutate(mod.size=max(n.state.nodes)) %>%
  select(layer_id, mod.size) %>%
  distinct() %>%
  filter(., !(layer_id %in% mod1$layer_id))
  #filter(layer_id !(%in% mod1$layer_id))

# Randomly sample the state nodes in the smallest modules
not.mod1.sampled <- not.mod1 %>%
  select(st.node.id, node_id) %>%
  #distinct(node_id) %>%
  group_by(node_id) %>%
  sample_n(., 1) %>%
  ungroup() %>%
  sample_n(., 10, replace=FALSE)  %>%
  left_join(., st.node.list.ord) %>%
  mutate(pre = mat.order - 1,
         true = 1, 
         post = 1514 - mat.order) %>%
  arrange(mat.order)

# Save the data
save(not.mod1.sampled, file="not.mod1.sampled.Rda")
####################################################################################################################


####################################################################################################################
# Section 4: CREATE STARTING BOOLEAN LIST FOR PERIPHERAL PLASMIDS
####################################################################################################################
# Repeat the process from Section 2 on the peripheral plasmids
# Create an empty list
bool.periph.obs <- list()

# Create a boolean list for each starting plasmid, with the position of the starting plasmid TRUE and all other
# positions FALSE
for(i in 1:nrow(not.mod1.sampled)) {
  bt <- c(rep(FALSE, each=not.mod1.sampled$pre[i]), rep(TRUE,1), rep(FALSE,each=not.mod1.sampled$post[i]))
  
  bool.periph.obs[[length(bool.periph.obs)+1]]<-bt
  
}

# Assign name of each list as the starting plasmid
names(bool.periph.obs) <- c(paste("periph",not.mod1.sampled$st.node.id, sep = "_"))

# Save the data 
save(bool.periph.obs, file="bool.periph.obs.Rda")
