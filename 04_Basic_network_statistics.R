####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 4: Basic network statistics (degree centrality, strength, etc) and visualizations
# 
# 
# The following outputs are used for downstream analysis: 
# plas.per.cow.Rda, cows.per.plasmid.Rda, deg.str.2k.all.phys.Rda, real.poten.edges.obs.Rda, layer_links.Rda
# 
# The following figures are created:
# Figure 1A pie chart
# 
# Script tested for R version 4.1.1
####################################################################################################################


####################################################################################################################
# SCRIPT SET-UP
####################################################################################################################
# Set working directory to wherever your files are located

# Load necessary packages:
library(tidyverse)
library(viridis)
library(ggrepel)

# Starting files:
# Edgelist with edgeweights created in script 03_Network_setup.R
load("net.dat2k.ew.Rda")

# List of state nodes (layers - nodes) created in script 03_Network_setup.R
load("plas.2k.Rda")
####################################################################################################################


####################################################################################################################
# Section 1: NETWORK DENSITY
#################################################################################################################### 
# M = # edges
# N = state nodes
M <- nrow(net.dat2k.ew)
N <- nrow(plas.2k)

density <- 2*M/(N*(N-1))
density*100
# 0.239
####################################################################################################################


####################################################################################################################
# Section 2: CALCULATE BASIC STATISTICS ON PLASMID/COW DISTRIBUTION
####################################################################################################################
# Calculate the number of cows per plasmid and the number of plasmids per cow

# Number of plasmids per cow
plas.per.cow <- plas.2k %>%
  ungroup() %>%
  group_by(layer_id) %>%
  summarise(num.plasmids = n())

# Save the data:
save(plas.per.cow, file="plas.per.cow.Rda")

# Number of cows each plasmid found in 
cows.per.plasmid <- plas.2k %>%
  group_by(node_id) %>%
  mutate(num.cows = n_distinct(layer_id)) %>%
  select(node_id, num.cows) %>%
  distinct()

# Save the data:
save(cows.per.plasmid, file="cows.per.plasmid.Rda")
####################################################################################################################


####################################################################################################################
# Section 3: INTER VS. INTRA EDGES - REALIZED VS POTENTIAL
####################################################################################################################
# Count realized inter- and intra-layer edges:
inter.intra.edg.count <- net.dat2k.ew %>%
  group_by(edge_type) %>%
  summarise(n=n()) %>%
  arrange(desc(edge_type)) %>%
  mutate(per.edg = n/2738*100) 

# Create pie chart in Figure 1A:
# Visualize breakdown in edge type with a pie chart with labels for number of edges
edge.pie <- ggplot(inter.intra.edg.count, aes(x="", y=n, fill=factor(edge_type))) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  labs(fill='Edge Type') +
  theme_void() +
  scale_fill_viridis(begin = 0.5, end=1, discrete = T) + 
  geom_label_repel(
    aes(y = n, label = paste0(n, " edges\n", round(per.edg,digits=1),"%")),
    size = 7, nudge_x = 2, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Edge Type")) +
  theme(legend.title=element_text(size=20), legend.text=element_text(size=18))

# Save pie chart from Figure 1A in desired format:
ggsave(edge.pie,filename="inter.intra.pie.png", dpi = 900, width = 8, height = 6, units = "in")
ggsave(edge.pie, filename="edge.pie.pdf", dpi = 900, width = 8, height = 6, units = "in")


# Determine what portion of *potential* intra- and inter-layer edges are realized
# Begin with intra-layer edges by calculating the potential number of intra-layer edges per cow:
intra.inter.possible <- plas.per.cow %>%
  group_by(layer_id) %>%
  mutate(poss.intra.cow = (num.plasmids*(num.plasmids -1)) / 2) %>%
  ungroup() %>%
  mutate(poss.intra = sum(poss.intra.cow)) %>%
  mutate(shuff.rep = 0) %>%
  select(shuff.rep, poss.intra) %>%
  distinct() %>%
  mutate(poss.inter = (1514*(1514-1))/2 - poss.intra)


# Bind the realized edges to the potential and calculate the percent realized and the ratio of potential:realized
real.poten.edges.obs <- intra.inter.possible %>%
  mutate(real.intra = 102,
         real.inter = 2636) %>%
  mutate(per.intra.real = real.intra/poss.intra*100,
         per.inter.real = real.inter/poss.inter*100,
         ratio.real = per.inter.real/per.intra.real)

# Save the data for downstream analysis
save(real.poten.edges.obs, file="real.poten.edges.obs.Rda")
####################################################################################################################


#########################################################################################################################
# Section 4: CALCULATE BASIC NETWORK METRICS: DEGREE
#################################################################################################################### 
# Prepare data
plas.2k.fr.wei <- net.dat2k.ew %>%
  ungroup() %>%
  select(layer_from, node_from, weight, edge_type) %>%
  mutate(layer_id = layer_from,
         node_id = node_from) %>%
  select(-layer_from, -node_from)

plas.2k.to.wei <- net.dat2k.ew %>%
  ungroup() %>%
  select(layer_to, node_to, weight, edge_type) %>%
  mutate(layer_id = layer_to,
         node_id = node_to) %>%
  select(-layer_to, -node_to) 

plas.2k.wei <- bind_rows(plas.2k.fr.wei, plas.2k.to.wei) %>%
  select(layer_id, node_id, weight, edge_type)

# Calculate (total) degree for each physical node:
deg.str.2k.all.phys <- plas.2k.wei %>%
  group_by(node_id) %>% 
  select(-layer_id) %>%
  summarise(degree = n(),
            strength = sum(weight)) %>%
  mutate(edge_type = "All")

# Save this data:
save(deg.str.2k.all.phys, file="deg.str.2k.all.phys.Rda")

# Degree for physical nodes, intra-layer edges:
pre.deg.str.intra.phys <- plas.2k.wei %>%
  filter(edge_type == "Intra") %>%
  group_by(node_id) %>% 
  summarise(degree = n(),
            strength = sum(weight)) %>%
  mutate(edge_type ="Intra")

# Identify nodes with no intra-layer edges
no.intra <- plas.2k %>%
  select(node_id) %>%
  distinct() %>%
  filter(!(node_id %in% pre.deg.str.intra.phys$node_id)) %>%
  mutate(degree = 0,
         strength = 0,
         edge_type = "Intra")

# Combine pre.deg.str.intra.phys and no.intra into a single dataframe
# that includes all plasmids:
deg.str.intra.phys <- pre.deg.str.intra.phys %>%
  bind_rows(., no.intra) 

# Degree for physical nodes, inter-layer edges (physical):
deg.str.inter.phys <- plas.2k.wei %>%
  filter(edge_type == "Inter") %>%
  group_by(node_id) %>% 
  summarise(degree = n(),
            strength = sum(weight)) %>%
  mutate(edge_type="Inter")

# Identify plasmids with no inter layer edges
no.inter <- plas.2k %>%
  select(node_id) %>%
  distinct() %>%
  filter(!(node_id %in% deg.str.inter.phys$node_id)) %>%
  mutate(degree = 0,
         strength = 0,
         edge_type = "Inter")

deg.str.inter.phys2 <- deg.str.inter.phys %>%
  bind_rows(., no.inter)

# Reformat (easier for certain visualizations, etc) by gathering degree and strength and labeling as inter/intra
deg.intra.inter.comp <- deg.str.inter.phys %>%
  bind_rows(., deg.str.intra.phys)

# Second format with separate columns for inter and intra degree / strength
# Edit column names and fill in 0's for nodes without any intra- or inter-layer edges   
deg.intra.inter.reformat <- deg.str.inter.phys %>%
  dplyr::rename(degree.inter = degree,
         strength.inter = strength) %>%
  full_join(., deg.str.intra.phys, by = "node_id") %>%
  dplyr::rename(degree.intra = degree,
         strength.intra = strength) %>%
  select(-edge_type.x, -edge_type.y) %>%
  replace(is.na(.), 0)
####################################################################################################################


####################################################################################################################
# Section 5: CORRELATIONS BETWEEN INTER AND INTRA DEGREE
####################################################################################################################
# Scatterplot comparing intra/inter degree.
# Note outlier
inter.intra.deg.scatter <- ggplot(deg.intra.inter.reformat, aes(x=degree.intra, y=degree.inter)) +
  geom_point() +
  xlab("Intra degree") + 
  ylab("Inter degree") +
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=14)) 

# Correlations between inter and intra degree
# Due to outlier, avoid Pearson's
cor(deg.intra.inter.reformat$degree.inter, deg.intra.inter.reformat$degree.intra, method = "kendall")
cor.test(deg.intra.inter.reformat$degree.inter, deg.intra.inter.reformat$degree.intra, method = "kendall")
####################################################################################################################


####################################################################################################################
# Section 6: CALCULATE LAYER LINKS, VISUALIZE DISTRIBUTION
####################################################################################################################
# Layer links (number of edges to plasmids in other cows): 
layer.link.setup <- net.dat2k.ew %>%
  ungroup() %>%
  select(-uniq_pair) %>%
  filter(node_from != node_to) %>%
  filter(edge_type=="Inter") 

layer.link.fr <- layer.link.setup %>%
  select(layer_from, node_from, weight, edge_type) %>%
  mutate(layer_id = layer_from,
         node_id = node_from) %>%
  select(-layer_from, -node_from)

layer.link.to <- layer.link.setup %>%
  select(layer_to, node_to, weight, edge_type) %>%
  mutate(layer_id = layer_to,
         node_id = node_to) %>%
  select(-layer_to, -node_to)

pre.layer_links <- layer.link.fr %>%
  bind_rows(., layer.link.to) %>%
  select(node_id, layer_id) %>%
  distinct() %>%
  group_by(node_id) %>% 
  summarise(links.to.cows = n())

no.layer.links <- plas.2k %>%
  select(node_id) %>%
  distinct() %>%
  filter(!(node_id %in% pre.layer_links$node_id)) %>%
  mutate(links.to.cows = 0)

layer_links <- pre.layer_links %>%
  rbind(no.layer.links)

# Save the data:
save(layer_links, file = "layer_links.Rda")
####################################################################################################################