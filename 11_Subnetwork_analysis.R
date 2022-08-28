####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 11: Analysis of subnetworks
# 
# 
# The following outputs are used for downstream analysis: 
# hgt.netdat.et.Rda, distdisp.netdat.et.Rda, recdisp.netdat.et.Rda,
# hgt.node.id.Rda, disp.node.id.Rda, recdisp.node.id.Rda, distdisp.node.id.Rda
# 
# The following figures are created:
# Figure 3A, network visualization of subnetworks
# Figure 3B, histogram of edge-weights
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

# Edgelist prior to edge-weight creation created in script 02_Threshold_sensitivity.R
load("thresh2000.a.Rda")

# Node id table created in script 08_Infomap_analysis_full_network.R
load("plas.2k.list.Rda")

# Infomap file created in script 08_Infomap_analysis_full_network
load("plas_mods.df.Rda")
####################################################################################################################


####################################################################################################################
# Section 1: IDENTIFY BREAK-POINTS IN EDGE-WEIGHTS FOR SUBNETWORKS
####################################################################################################################
# Label edges by subnetwork:
net.dat2k.ew.col <- net.dat2k.ew %>%
  mutate(Subnetwork = case_when(
    weight <= 0.5 ~ "HGT",
    weight > 0.5 & weight < 0.95 ~ "Distant dispersal",
    weight >= 0.95 ~ "Recent dispersal"
  )) %>%
  ungroup() %>%
  mutate(Subnetwork = factor(Subnetwork, levels=c("HGT", "Distant dispersal","Recent dispersal")))

# Visualize the distribution of edge-weights, identify breaks
# Create Figure 3B: Histogram of edge-weights
ew.hist.col <- ggplot(net.dat2k.ew.col, aes(x=weight, fill=Subnetwork))+
  geom_histogram(color="black", breaks=seq(0.2,1,0.05)) +
  scale_fill_manual(values=c("indianred1", "thistle2","slategray2")) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.y=element_text(size=24),
        axis.title.x=element_text(size=24)) +
  labs(x = "Edge Weight", y="Count (Plasmids)" ) +
  geom_vline(xintercept = 0.5, color = "gray59", size=2) +
  geom_vline(xintercept = 0.95, color = "gray59", size=2)

# Save in desired format
ggsave(ew.hist.col, filename="ew.distr.colorbar.png", dpi = 900, width = 8, height = 6, units = "in")
ggsave(ew.hist.col, filename="ew.distr.colorbar.pdf", width = 8, height = 6, units = "in")
####################################################################################################################


####################################################################################################################
# Section 2: RANDOMIZE EDGE-WEIGHT DISTRIBUTION BY SHUFFLING ALIGNMENT LENGTH BETWEEN PLASMIDS
####################################################################################################################
set.seed(123)

# Function to shuffle length of alignment between plasmids
# NOTE: Uses thresh2000.a dataframe, NOT the network edgelist, in order to recalculate the edge-weight 
# with the shuffled alignment length.
shuff.len.fun <- function(x) {
  transform(thresh2000.a, length = sample(length, replace = F) )
}

# Choose how many times to shuffle the data using the function created above.
# Creates a list with each shuffling
shuff.len.list <- (lapply(seq(1000),shuff.len.fun))

# Give each shuffled list a name corresponding to the number of its shuffling
names(shuff.len.list) <- seq(1000)


# Create a new column that ID's the shuffling repetition number of each iteration
# Bind each list together into a single dataframe. 
shuff.len.df <- Map(cbind, shuff.len.list, shuff.rep = names(shuff.len.list)) %>%
  bind_rows()

# Recalculate edge-weights on this reshuffled dataframe
ew.len.shuff <- shuff.len.df %>%
  group_by(uniq_pair, shuff.rep) %>% 
  mutate(align.proportion.query = length/plasmid.length.query,
         align.proportion.sub = length/plasmid.length.subject ) %>%
  rowwise() %>%
  # Calculate the maximum proportion of plasmid length covered by the alignment
  # E.g. the alignment length / length of smaller plasmid
  mutate(min.align.prop = min(align.proportion.query, align.proportion.sub)) %>%
  rowwise() %>%
  mutate(preweight = min.align.prop*(pident/100)) %>%
  group_by(uniq_pair, shuff.rep) %>%   
  mutate(weight = sum(preweight)) %>%
  select(-min.align.prop, -max.align.prop, -qstart, -qend, -sstart, -send, -evalue,-qs.qe.diff,
         -align.proportion.query,-align.proportion.sub) %>%
  select(-preweight) %>%
  mutate(dat.type = "Shuffled") %>%
  distinct()

# Label and bind the observed data and shuffled data
plot.len.dat <- net.dat2k.ew %>%
  mutate(dat.type = "Observed",
         shuff.rep = "0") %>%
  bind_rows(., ew.len.shuff)

# Create Figure S3, comparing the density of edge-weights for the observed data vs the shuffled data:
dens.len.plot <- ggplot(plot.len.dat, aes(x=weight, color=dat.type, fill=dat.type)) + # if you want to scale y, add argument: "y=..scaled.." to aes
  geom_density(alpha=0.1, size=1)+
  labs(x="Edge-weight", y = "Density", color=NULL,fill=NULL) +
  theme_classic() +
  scale_fill_viridis(end =0.5, discrete = T) +
  scale_color_viridis(end =0.5, discrete=T) +
  xlim(0,3) +
  theme(legend.text=element_text(size=18),
        axis.title=element_text(size=20),
        axis.text = element_text(size=16))

# Save in desired format
ggsave(filename="obs.shuff.ew.png", dpi = 1800, width = 8, height = 6, units = "in")
ggsave(filename="obs.shuff.ew.pdf", dpi = 900, width = 8, height = 6, units = "in")

# Compare the distribution of edge-weights between observed and shuffled (alignment length) data
# with a statistical test (Kolmogorov-Smirnov)
shuff.len.ks <- ks.test(net.dat2k.ew$weight, ew.len.shuff$weight)
####################################################################################################################


####################################################################################################################
# Section 3: DIVIDE THE NETWORK INTO SUBNETWORKS BASED ON EDGE-WEIGHTS
####################################################################################################################
# Break up data into 3 subnetworks: HGT, recent dispersal, distant dispersal
# based on the edge-weights

# HGT network (EW < 0.5)
hgt.netdat.et <- net.dat2k.ew %>%
  filter(weight <= 0.5) %>%
  ungroup() %>%
  select(uniq_pair,layer_from, node_from, layer_to, node_to, edge_type, weight)

# Distant dispersal (0.5 < EW < 0.95)
distdisp.netdat.et <- net.dat2k.ew %>%
  filter(weight > 0.5 & weight < 0.95) %>%
  ungroup() %>%
  select(uniq_pair,layer_from, node_from, layer_to, node_to, edge_type, weight)

# Recent dispersal (EW > 0.95)
recdisp.netdat.et <- net.dat2k.ew %>%
  filter(weight >= 0.95) %>%
  ungroup() %>%
  select(uniq_pair,layer_from, node_from, layer_to, node_to, edge_type, weight)

# Save each subnetwork for downstream analysis
save(hgt.netdat.et, file="hgt.netdat.et.Rda")
save(distdisp.netdat.et, file="distdisp.netdat.et.Rda")
save(recdisp.netdat.et, file="recdisp.netdat.et.Rda")
####################################################################################################################


####################################################################################################################
# Section 4: NODE ID TABLE FOR EACH SUBNETWORK
####################################################################################################################
# Node ID table for each subnetwork:
hgt.node.id <- hgt.netdat.et %>%
  ungroup() %>%
  select(node_from, node_to) %>%
  gather(.) %>%
  select(-key) %>%
  mutate(node_id = as.integer(value)) %>%
  select(-value) %>%
  distinct() %>%
  left_join(., plas.2k.list)

recdisp.node.id <- recdisp.netdat.et %>%
  ungroup() %>%
  select(node_from, node_to) %>%
  gather(.) %>%
  select(-key) %>%
  mutate(node_id = as.integer(value)) %>%
  select(-value) %>%
  distinct() %>%
  left_join(., plas.2k.list)

distdisp.node.id <- distdisp.netdat.et %>%
  ungroup() %>%
  select(node_from, node_to) %>%
  gather(.) %>%
  select(-key) %>%
  mutate(node_id = as.integer(value)) %>%
  select(-value) %>%
  distinct() %>%
  left_join(., plas.2k.list)

# Save this data:
save(hgt.node.id, file = "hgt.node.id.Rda")
save(recdisp.node.id, file = "recdisp.node.id.Rda")
save(distdisp.node.id, file = "distdisp.node.id.Rda")
####################################################################################################################


####################################################################################################################
# Section 5: FLOW PER SUBNETWORK
####################################################################################################################
# Join the node id's to the infomap data to get flow per plasmid:
plas.flow <- plas_mods.df %>%
  select(node_id, flow) %>%
  distinct()

# Calculate flow of plasmids in HGT subnetwork
hgt.flow <- hgt.node.id %>%
  left_join(., plas.flow)

sum(hgt.flow$flow)

# Calculate flow of plasmids in recent dispersal subnetwork
rec.flow <- recdisp.node.id %>%
  left_join(., plas.flow)

sum(rec.flow$flow)

# Calculate flow of plasmids in distant dispesal subnetwork
dist.flow <- distdisp.node.id %>%
  left_join(., plas.flow)

sum(dist.flow$flow)
####################################################################################################################


####################################################################################################################
# Section 6: COUNTS OF INTER- TO INTRA-LAYER EDGES IN EACH SUBNETWORK
####################################################################################################################
# HGT
hgt.edg.df <- hgt.netdat.et %>%
  count(uniq_pair, edge_type) %>%
  spread(., edge_type, n, fill = 0, 
         drop = FALSE, sep = NULL)

hgt.intra.rat <- hgt.edg.df %>%
  ungroup() %>%
  count(Intra > 0)

hgt.intra <- 11


hgt.inter.rat <- hgt.edg.df %>%
  ungroup() %>%
  count(Inter > 0)

hgt.inter <- 223

hgt.rat <- hgt.inter/hgt.intra


# Recent dispersal
recdisp.edg.df <- recdisp.netdat.et %>%
  count(uniq_pair, edge_type) %>%
  spread(., edge_type, n, fill = 0, 
         drop = FALSE, sep = NULL)

recdisp.intra.rat <- recdisp.edg.df %>%
  ungroup() %>%
  count(Intra > 0)

recdisp.intra <- 24

recdisp.inter.rat <- recdisp.edg.df %>%
  ungroup() %>%
  count(Inter > 0)

recdisp.inter <- 496

recdisp.rat <- recdisp.inter/recdisp.intra


# Distant dispersal
distdisp.edg.df <- distdisp.netdat.et %>%
  count(uniq_pair, edge_type) %>%
  spread(., edge_type, n, fill = 0, 
         drop = FALSE, sep = NULL)

distdisp.intra.rat <- distdisp.edg.df %>%
  ungroup() %>%
  count(Intra > 0)

distdisp.intra <- 67 

distdisp.inter.rat <- distdisp.edg.df %>%
  ungroup() %>%
  count(Inter > 0)

distdisp.inter <- 1917
####################################################################################################################


####################################################################################################################
# Section 7: IDENTIFY PLASMIDS IN MULTIPLE SUBNETWORKS
####################################################################################################################

# HGT + Recent Dispersal
hgt.recdisp <- hgt.node.id %>%
  inner_join(., recdisp.node.id)

nrow(hgt.recdisp)
# 35 shared plasmids

# HGT + Distant Dispersal
hgt.distdisp <- hgt.node.id %>%
  inner_join(., distdisp.node.id)

nrow(hgt.distdisp)
# 133 shared plasmids

# Recent Dispersal + Distant dispersal
rec.distdisp <- recdisp.node.id %>%
  inner_join(., distdisp.node.id)

nrow(rec.distdisp)
# 293 shared plasmids

# In all three
all.three.subs <- hgt.node.id %>%
  inner_join(., rec.distdisp)

nrow(all.three.subs)
# 24 shared plasmids
####################################################################################################################


####################################################################################################################
# Section 8: VISUALIZE EACH SUBNETWORK WITH IGRAPH
####################################################################################################################
# HGT network:
# Set up interlayer edges for network edge visualization
hgt.inter.vis <- hgt.netdat.et %>%
  filter(edge_type=="Inter") %>%
  ungroup() %>%
  group_by(layer_from, layer_to) %>%
  mutate(weight = n()) %>%
  mutate(layer_from = as.numeric(layer_from),
         layer_to = as.numeric(layer_to)) %>%
  select(layer_from, layer_to, weight) %>%
  distinct() %>%
  ungroup

# Node list
hgt.inter.list <- hgt.inter.vis %>%
  select(layer_from, layer_to) %>%
  gather() %>%
  select(value) %>%
  distinct() %>%
  dplyr::rename(layer=value)

# Intralayer edges for node size
hgt.intra.vis <- hgt.netdat.et %>%
  filter(edge_type=="Intra") %>%
  ungroup() %>%
  group_by(layer_from, layer_to) %>%
  mutate(intra.edg = n()) %>%
  dplyr::rename(layer=layer_from) %>%
  ungroup() %>%
  select(layer, intra.edg) %>%
  bind_rows(., hgt.inter.list) %>%
  distinct() %>%
  group_by(layer) %>%
  arrange(desc(intra.edg)) %>%
  slice(1) %>%
  distinct() %>%
  replace(is.na(.), 0)

# Create a graph
graph.hgt <- graph_from_data_frame(hgt.inter.vis, directed=F, vertices = hgt.intra.vis)

# Set edge weight attribute for visualization
set_edge_attr(graph.hgt, "weight", value= hgt.inter.vis$weight)
is.weighted(graph.hgt)

# Weight node size by number of intralayer edges per node
V(graph.hgt)$size <- 5+hgt.intra.vis$intra.edg

# Plot network
plot(graph.hgt, edge.width=(E(graph.hgt)$weight*0.1),vertex.label=NA, layout=layout_with_lgl, vertex.color="red",edge.color="black")

# Save coordinates of first (HGT) network and apply to all visualizations.
# In this way, each cow will be in the same place for each subnetwork visualization.
Coords <- layout_with_fr(graph.hgt)


# Recent dispersal
# Set up interlayer edges for network edge visualization
rec.inter.vis <- recdisp.netdat.et %>%
  filter(edge_type=="Inter") %>%
  ungroup() %>%
  group_by(layer_from, layer_to) %>%
  mutate(weight = n()) %>%
  mutate(layer_from = as.numeric(layer_from),
         layer_to = as.numeric(layer_to)) %>%
  select(layer_from, layer_to, weight) %>%
  distinct() %>%
  ungroup

# Node list
rec.inter.list <- rec.inter.vis %>%
  select(layer_from, layer_to) %>%
  gather() %>%
  select(value) %>%
  distinct() %>%
  dplyr::rename(layer=value)

# Intralayer edges for node size
rec.intra.vis <- recdisp.netdat.et %>%
  filter(edge_type=="Intra") %>%
  ungroup() %>%
  group_by(layer_from, layer_to) %>%
  mutate(intra.edg = n()) %>%
  dplyr::rename(layer=layer_from) %>%
  ungroup() %>%
  select(layer, intra.edg) %>%
  bind_rows(., rec.inter.list) %>%
  distinct() %>%
  group_by(layer) %>%
  arrange(desc(intra.edg)) %>%
  slice(1) %>%
  distinct() %>%
  replace(is.na(.), 0)

# Create a graph from dataframe
graph.rec <- graph_from_data_frame(rec.inter.vis, directed=F, vertices = rec.intra.vis)

# Set edge attributes by weight
set_edge_attr(graph.rec, "weight", value= rec.inter.vis$weight)
is.weighted(graph.rec)

# Set node attribute
V(graph.rec)$size <- 5+rec.intra.vis$intra.edg

# Plot network
plot(graph.rec, edge.width=(E(graph.rec)$weight*0.4),vertex.label=NA, layout=Coords, vertex.color="slategray2", edge.color="black")


# Distant dispersal
# Set up interlayer edges for network edge visualization
dist.inter.vis <- distdisp.netdat.et %>%
  filter(edge_type=="Inter") %>%
  ungroup() %>%
  group_by(layer_from, layer_to) %>%
  mutate(weight = n()) %>%
  mutate(layer_from = as.numeric(layer_from),
         layer_to = as.numeric(layer_to)) %>%
  select(layer_from, layer_to, weight) %>%
  distinct() %>%
  ungroup

# Node list
dist.inter.list <- dist.inter.vis %>%
  select(layer_from, layer_to) %>%
  gather() %>%
  select(value) %>%
  distinct() %>%
  dplyr::rename(layer=value)

# Intralayer edges for node size
dist.intra.vis <- distdisp.netdat.et %>%
  filter(edge_type=="Intra") %>%
  ungroup() %>%
  group_by(layer_from, layer_to) %>%
  mutate(intra.edg = n()) %>%
  dplyr::rename(layer=layer_from) %>%
  ungroup() %>%
  select(layer, intra.edg) %>%
  bind_rows(., dist.inter.list) %>%
  distinct() %>%
  group_by(layer) %>%
  arrange(desc(intra.edg)) %>%
  slice(1) %>%
  distinct() %>%
  replace(is.na(.), 0)

# Create graph from dataframe
graph.dist <- graph_from_data_frame(dist.inter.vis, directed=F, vertices = dist.intra.vis)

# Set edge weight attribute
set_edge_attr(graph.dist, "weight", value= dist.inter.vis$weight)
is.weighted(graph.dist)

# Intralayer edges for node size
V(graph.dist)$size <- 5+dist.intra.vis$intra.edg

# Plot network
plot(graph.dist, edge.width=(E(graph.dist)$weight*0.4),vertex.label=NA, layout=Coords, vertex.color="thistle2", edge.color="black")

# Create Figure 3A: Subnetwork plots
# (Each subnetwork plotted separately as a pdf)
# Save pdfs:
pdf("hgt.network.plot.pdf", width = 4, height = 4)
{
  plot(graph.hgt, edge.width=(E(graph.hgt)$weight*0.2),vertex.label=NA, layout=Coords, vertex.color="indianred1",edge.color="black")
  
}
dev.off()

pdf("hgt.network.plot.pdf", width = 4, height = 4)
{
  plot(graph.hgt, edge.width=(E(graph.hgt)$weight*0.2),vertex.label=NA, layout=Coords, vertex.color="indianred1",edge.color="black",
       vertex.frame.width=0.1)
  
}
dev.off()

# Recent
pdf("rec.network.plot.pdf", width = 4, height = 4)
{
  plot(graph.rec, edge.width=(E(graph.rec)$weight*0.2),vertex.label=NA, layout=Coords, vertex.color="slategray2", edge.color="black",
       vertex.frame.width=0.1)
}
dev.off()

# Distant
pdf("dist.network.plot.pdf", width = 5, height = 5) 
{
  plot(graph.dist, edge.width=(E(graph.dist)$weight*0.2),vertex.label=NA, layout=Coords, vertex.color="thistle2", edge.color="black",
       vertex.frame.width=0.1)
  
}
dev.off()
####################################################################################################################