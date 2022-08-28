####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 8: Multilayer network analysis on full network using Infomap
#           
# 
# The following outputs are used for downstream analysis: 
# net.2k.edgelist.Rda, plas.2k.list.Rda, layer.id.table.2k.Rda, plas_mods.df.Rda, plas_mods.stats.Rda
# 
# The following figures are created:
# Figure 1B, simplified network visualization
# 
# Script tested for R version 4.1.1
####################################################################################################################


####################################################################################################################
# SCRIPT SET-UP
####################################################################################################################
# Set working directory to wherever your files are located

# Load necessary packages:
library(tidyverse)
library(infomapecology)
library(igraph)

# Starting files:
# Edgelist with edgeweights created in script 03_Network_setup.R
load("net.dat2k.ew.Rda")

# Metadata file created in script 03_Network_setup.R
load("plasmid.metadat.updated.Rda")

# Plasmids per cow file created in in script 03_Network_setup.R
load("plas.per.cow.Rda")
####################################################################################################################


####################################################################################################################
# Section 1: SET-UP DATA FOR INFOMAP
####################################################################################################################
# Create the necessary dataframes for infomap
# First, a simplified extended edgelist in the right order:
net.2k.edgelist <- net.dat2k.ew %>%
  ungroup() %>%
  select(layer_from, node_from, layer_to, node_to, weight)

# Save the data
save(net.2k.edgelist, file="net.2k.edgelist.Rda")

# Make a metadata file, node id table and layer table for the included plasmids
# Start by making a full list of the included plasmids
plas.2k.list <- net.2k.edgelist %>%
  ungroup() %>%
  select(node_from, node_to) %>%
  gather(.) %>%
  select(-key) %>%
  mutate(node_id = as.integer(value)) %>%
  select(-value) %>%
  left_join(plasmid.metadat.updated, by = "node_id") %>%
  select(node_id, plasmid_length) %>%
  distinct()

# Create layer id table
layer.id.table.2k <- net.2k.edgelist %>%
  ungroup() %>%
  select(layer_from, layer_to) %>%
  gather(.) %>%
  select(value) %>%
  rename(layer_id= value) %>%
  distinct() %>%
  droplevels()

# Save for downstream analysis: 
save(plas.2k.list, file="plas.2k.list.Rda")
save(layer.id.table.2k, file="layer.id.table.2k.Rda")
####################################################################################################################


####################################################################################################################
# Section 2: RUN INFOMAP ON FULL NETWORK
####################################################################################################################
# Making multilayer network with infomapecology package
# Default:
plas.multilay.2k <- create_multilayer_object(extended = net.2k.edgelist, 
                                    nodes = plas.2k.list,  layers = layer.id.table.2k)


#Run infomap
plas_modules <- run_infomap_multilayer(M=plas.multilay.2k,infomap_executable = "Infomap",
                                       flow_model = 'undirected', silent = T, relax = F,
                                       trials = 100, seed = 123, temporal_network = F, run_standalone = T)
####################################################################################################################


####################################################################################################################
# Section 3: STATISTICAL ANALYSIS ON INFOMAP MODULES
####################################################################################################################
# Turn infomap output into a dataframe and reorganize slightly:
plas_mods.df <- as.data.frame(plas_modules$modules) %>%
  select(module, everything())

# Calculate the number of nodes and layers per module:
plas_mods.stats <- plas_mods.df %>%
  group_by(module) %>%
  summarise(n.state.nodes=n(),
            n.phys.nodes = n_distinct(node_id),
            n.layers = n_distinct(layer_id),
            mod.flow=sum(flow)) 

# Save data for downstream analysis:
save(plas_mods.df, file="plas_mods.df.Rda")
save(plas_mods.stats, file="plas_mods.stats.Rda")

# Visualize module characteristics with histograms:
# Distribution of state nodes per module
st.nodes.mod <- ggplot(plas_mods.stats, aes(x=n.state.nodes))+
  geom_histogram(color="black", fill="white") +
  labs(x = "Number of State Nodes", y="Count(Modules)") +
  theme(text = element_text(size=14)) +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))

# Distribution of physical nodes per module
phys.nodes.mod <- ggplot(plas_mods.stats, aes(x=n.phys.nodes))+
  geom_histogram(color="black", fill="white") +
  scale_x_continuous(breaks = seq(2, 13, 3)) +
  labs(x = "Number of Physical Nodes", y="Count (Modules)" ) +
  theme(text = element_text(size=16)) +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))

# Distribution of layers per module
layer.mod <- ggplot(plas_mods.stats, aes(x=n.layers))+
  geom_histogram(color="black", fill="white") +
  scale_x_continuous(breaks = seq(1, 16, 3)) +
  labs(x = "Number of Layers", y="Count (Modules)" ) +
  theme(text = element_text(size=14)) +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
####################################################################################################################


####################################################################################################################
# Section 4: DATA FORMATTING FOR NETWORK VISUALIZATIONS
####################################################################################################################
# Format data for visualizations:
plas_mods.stats.lab <- plas_mods.stats %>%
  mutate(mod.type=ifelse(module == 1, "big", "not big")) 
  
plas_mods2 <- plas_mods.df %>%
  left_join(., plas_mods.stats.lab, by="module") %>%
  select(layer_id,node_id, mod.type) %>%
  filter(mod.type=="big")

plas_mods3 <- net.dat2k.ew %>%
  left_join(., plas_mods2, by=c("node_from"="node_id","layer_from"="layer_id")) %>%
  rename(mod.from = mod.type) %>%
  left_join(., plas_mods2, by=c("node_to"="node_id","layer_to"="layer_id")) %>%
  rename(mod.to = mod.type) %>%
  mutate(inter.big=ifelse(mod.from=="big" & mod.to=="big", "big.to.big", "not big"))

# Biggest module, intra-module edges
big.to.big <- plas_mods3 %>%
  filter(inter.big == "big.to.big")

# Biggest module, inter-module edges
big.to.not.big <- plas_mods3 %>%
  filter(mod.from=="big" & is.na(mod.to) | is.na(mod.from) & mod.to=="big")

# Total edges
tot.big.edges <- big.to.big %>%
  bind_rows(., big.to.not.big)

# To get a list of all plasmids involved in the biggest module:
big.plas_mods.from <- plas_mods3 %>%
  ungroup() %>%
  select(layer_from, node_from, mod.from) %>%
  filter(mod.from=="big") %>%
  dplyr::rename(mod.lab = mod.from,
         layer_id = layer_from,
         node_id = node_from)

big.plas_mods.to <- plas_mods3 %>%
  ungroup() %>%
  select(layer_to, node_to,mod.to) %>%
  filter(mod.to=="big") %>%
  dplyr::rename(mod.lab = mod.to,
         layer_id= layer_to,
         node_id = node_to)

big.mod.plas.list <- big.plas_mods.from %>%
  rbind(., big.plas_mods.to) %>%
  select(-layer_id) %>%
  distinct()

big.mod.plas.list.layers <- big.plas_mods.from %>%
  rbind(., big.plas_mods.to) %>%
  distinct()

not.big.plas_mods.from <- plas_mods3 %>%
  ungroup() %>%
  select(layer_from, mod.from) %>%
  filter(is.na(mod.from)) %>%
  rename(mod.lab = mod.from,
         layer_id = layer_from) %>%
  mutate(mod.lab = replace_na(mod.lab, "not big"))


not.big.plas_mods.to <- plas_mods3 %>%
  ungroup() %>%
  select(layer_to, mod.to) %>%
  filter(is.na(mod.to)) %>%
  rename(mod.lab = mod.to,
         layer_id = layer_to) %>%
  mutate(mod.lab = replace_na(mod.lab, "not big")) %>%
  rbind(., not.big.plas_mods.from) %>%
  distinct()

only.not.big <- not.big.plas_mods.to %>%
  filter(!(layer_id %in% big.mod.plas.list.layers$layer_id))

mod.layer.list <- big.mod.plas.list.layers %>%
  select(-node_id) %>%
  rbind(only.not.big) %>%
  mutate(layer_id = as.factor(layer_id))


# Set up data for visualizations:
# Labelling edges as part of the largest module or not:
inter.for.vis2 <- plas_mods3 %>%
  filter(edge_type=="Inter") %>%
  ungroup() %>%
  mutate(l_fr = paste("L",layer_from,sep = "")) %>%
  mutate(l_to = paste("L",layer_to, sep = "")) %>%
  rowwise() %>%
  mutate(uniq_laypair = paste(sort(c(l_fr, l_to)), collapse = "_")) %>%
  select(-l_fr, -l_to) %>%
  select(uniq_laypair, everything()) %>%
  group_by(uniq_laypair, inter.big) %>%
  mutate(weight = n()) %>%
  mutate(layer_from = as.numeric(layer_from),
         layer_to = as.numeric(layer_to)) %>%
  mutate_at(vars(inter.big), ~replace_na(., "not.big")) %>%
  group_by(uniq_laypair, inter.big) %>%
  slice(1) %>%
  distinct() %>%
  ungroup %>%
  select(layer_from, layer_to, weight, inter.big)

# Node list
inter.list2 <- inter.for.vis2 %>%
  ungroup() %>%
  select(layer_from, layer_to) %>%
  gather() %>%
  select(value) %>%
  distinct() %>%
  rename(layer=value)

# Intra-layer values for visualization:
intra.for.vis2 <- plas_mods3 %>%
  filter(edge_type=="Intra") %>%
  ungroup() %>%
  group_by(layer_from, layer_to) %>%
  mutate(intra.edg = n()) %>%
  rename(layer=layer_from) %>%
  ungroup() %>%
  select(layer, intra.edg) %>%
  bind_rows(.,inter.list2) %>%
  group_by(layer) %>%
  arrange(desc(intra.edg)) %>%
  slice(1) %>%
  distinct() %>%
  mutate_at(vars(intra.edg), ~replace_na(., 0)) %>%
  mutate(layer = as.factor(layer)) %>%
  left_join(., mod.layer.list, by=c("layer"="layer_id")) %>%
  left_join(., plas.per.cow, by=c("layer"="layer_id")) %>%
  distinct()
####################################################################################################################


####################################################################################################################
# Section 5: NETWORK VISUALIZATIONS
####################################################################################################################
# Visualization of the network with interlayer edges and the largest module highlighted
graph2 <- graph_from_data_frame(inter.for.vis2, directed=F, vertices = intra.for.vis2)

set_edge_attr(graph2, "weight", value= inter.for.vis2$weight)
is.weighted(graph2)

V(graph2)$size <- 12+intra.for.vis2$intra.edg
V(graph2)$color <- intra.for.vis2$mod.lab

E(graph2)$color <- as.factor(inter.for.vis2$inter.big)
V(graph2)$color <- as.factor(intra.for.vis2$mod.lab)

E(graph2)$color <- ifelse(inter.for.vis2$inter.big == "big.to.big", "red", "gray")
V(graph2)$color <- ifelse(intra.for.vis2$mod.lab == "big", "red", "gray")

# Test plotting
plot(graph2, edge.width=(E(graph2)$weight*0.1),vertex.label=intra.for.vis2$num.plasmids, 
     vertex.label.font=2, vertex.label.cex=0.9,
     vertex.label.color = "black", layout=layout_with_lgl)


# Create Figure 2 Panel B: Full network plot with largest module highlighted:
png(filename = "big.mod.highlight.plot.labs.png",res = 900, width = 10, height = 10, units = "in") 
{
  plot(graph2, edge.width=(E(graph2)$weight*0.1),vertex.label=intra.for.vis2$num.plasmids, 
       vertex.label.font=2, vertex.label.cex=0.9 , 
       vertex.label.color = "black",layout=layout_with_lgl)
}
dev.off()

# Save as pdf
pdf(file = "big.mod.highlight.plot.pdf", width = 6, height = 6) 
{
  plot(graph2, edge.width=(E(graph2)$weight*0.1),vertex.label=intra.for.vis2$num.plasmids, 
       vertex.label.font=2, vertex.label.cex=0.9 , 
       vertex.label.color = "black",layout=layout_with_lgl)
}
dev.off()
####################################################################################################################


####################################################################################################################
# Section 6: LARGEST MODULE CHARACTERISTICS
####################################################################################################################
# Compare flow of the largest module vs rest of modules in the observed network with a Z-score
z.sc.flow.big <- (max(plas_mods.stats$mod.flow) - mean(plas_mods.stats$mod.flow)) / sd(plas_mods.stats$mod.flow)

# P value
2*pnorm(q=z.sc.flow.big, lower.tail=FALSE)

# Inter- and intra-module edges
big.mod2 <- plas_mods3 %>%
  filter(inter.big == "big.to.big") %>%
  ungroup()

big.mod.intra <- big.mod2 %>%
  count(edge_type=="Intra")

big.mod.inter <- big.mod2 %>%
  count(edge_type=="Inter")
####################################################################################################################