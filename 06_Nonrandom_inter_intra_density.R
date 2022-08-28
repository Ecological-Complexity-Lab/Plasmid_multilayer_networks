######################################################################################################################################
# Plasmid rumen network analysis
# 
# Script 6: Comparison of intra- and inter-layer edge density in full observed vs shuffled networks 
#           Determines if density and ratio of realized inter:intra edge density are non-random
# 
# The following figures are created:
# Supplementary Figure S1
#
# Script tested for R version 4.1.1
######################################################################################################################################


######################################################################################################################################
# SCRIPT SET-UP
######################################################################################################################################
# Set working directory to wherever your files are located

# Load necessary packages:
library(tidyverse)
library(ggpubr)
library(viridis)

# Load starting files:
# List of shuffled edgelists created in script 05_Network_shuffling.R
load("shuff.lf.list.Rda")

# Load density of observed data created in script 04_Basic_network_statistics
load("real.poten.edges.obs.Rda") 

######################################################################################################################################
# Section 1: REALIZED POTENTIAL EDGES IN OBSERVED VS SHUFFLED FULL NETWORK
######################################################################################################################################
# Compare the realized potential edges in the observed vs. shuffled FULL network

# Bind the lists of each shuffled dataframe into a single dataframe:
shuff.df <- Map(cbind, shuff.lf.list, shuff.rep = names(shuff.lf.list)) %>%
  bind_rows() %>%
  mutate(shuff.rep = as.integer(shuff.rep)) %>%
  mutate(edge_type=ifelse(layer_from==layer_to,'Intra','Inter')) 

# Count the realized number of intra and inter layer edges per shuffled network
shuff.edg.ct <- shuff.df %>%
  mutate(shuff.rep = as.integer(shuff.rep)) %>%
  group_by(shuff.rep) %>%
  mutate(real.inter = sum(edge_type=="Inter"),
         real.intra = sum(edge_type=="Intra")) %>%
  select(shuff.rep, real.intra, real.inter) %>%
  distinct() 

# To determine the number of plasmids per cow, first combine the "from","to", combine, then remove duplicates:
# State nodes "from"
plas.2k.fr.sh <- shuff.df %>%
  group_by(shuff.rep) %>%
  select(shuff.rep,layer_from, node_from) %>%
  mutate(layer_id = layer_from,
         node_id = node_from) %>%
  select(-layer_from, -node_from)


# State nodes "to"
plas.2k.to.sh <- shuff.df %>%
  group_by(shuff.rep) %>%
  select(shuff.rep,layer_to, node_to) %>%
  mutate(layer_id = layer_to,
         node_id = node_to) %>%
  select(-layer_to, -node_to)


# Combine - All state nodes
plas.2k.sh <- plas.2k.fr.sh %>%
  bind_rows(., plas.2k.to.sh) %>%
  distinct()

#Sanity check - make sure there are 1344 plasmids per shuff.rep
plas.per.rep <- plas.2k.sh %>%
  select(shuff.rep, node_id) %>%
  group_by(shuff.rep) %>%
  summarise(num.plasmids = n_distinct(node_id), .groups = "keep")
# Good, sanity preserved, all have 1344


# Plasmids per cow per shuffling rep:
plas.per.cow.sh <- plas.2k.sh %>%
  group_by(shuff.rep, layer_id) %>%
  summarise(num.plasmids = n(), .groups = "keep")

# Number of cows per plasmid
cows.per.plas.sh <- plas.2k.sh %>%
  group_by(shuff.rep, node_id) %>%
  summarise(num.cows = n(), .groups = "keep")

# State nodes per shuffled network
st.nodes.sh <- plas.2k.sh %>%
  group_by(shuff.rep) %>%
  distinct() %>%
  summarise(st.nodes = n(), .groups = "keep")  

# Join with the plasmids per cow per shuffled network
nodes.cows.sh <- plas.per.cow.sh %>%
  full_join(., st.nodes.sh, by = "shuff.rep")

# Calculate the number of potential intra and inter layer edges per shuffled network
intra.inter.possible.sh <- nodes.cows.sh %>%
  group_by(shuff.rep, layer_id) %>%
  mutate(poss.intra.cow = (num.plasmids*(num.plasmids -1)) / 2) %>%
  ungroup() %>%
  group_by(shuff.rep) %>%
  mutate(poss.intra = sum(poss.intra.cow)) %>%
  group_by(shuff.rep) %>%
  mutate(poss.inter = (st.nodes*(st.nodes-1))/2 - poss.intra) %>%
  distinct()

# Bind the realized and potential inter and intra layer edges for the shuffled networks
real.poten.edges.sh <- shuff.edg.ct %>%
  full_join(., intra.inter.possible.sh, by="shuff.rep") %>%
  mutate(per.intra.real = real.intra/poss.intra*100,
         per.inter.real = real.inter/poss.inter*100,
         ratio.real = per.inter.real/per.intra.real)

# Compare percent of potential edges realized, intra vs. intra:
# First, combine inter/intra realized in one column
# Shuffled networks:
real.poten.edges.sh2 <- real.poten.edges.sh %>%
  rename("Intra-Layer" = per.intra.real, # Note: Previously " Intra-/Inter-Layer Edges (%)"
         "Inter-Layer" = per.inter.real) %>%
  gather(starts_with("Int"), key="edge_type", value="percent.realized")

# Observed network:
real.poten.edges.obs2 <- real.poten.edges.obs %>%
  rename("Intra-Layer" = per.intra.real,
         "Inter-Layer" = per.inter.real) %>%
  gather(starts_with("Int"), key="edge_type", value="percent.realized")

# Create Supplementary Figure S1 Panel A: Comparing density of inter- and intra-layer edges in the observed 
# and shuffled networks:
per.intra.inter.real <- ggplot(real.poten.edges.sh2, aes(x=percent.realized, color=edge_type, fill=edge_type))+
  geom_histogram(alpha = 0.4, position="identity") +
  labs(x = "Density (%)", y="Edges (Count)") + 
  geom_vline(data=real.poten.edges.obs2, aes(xintercept=percent.realized, color=edge_type),
             linetype="dashed", size = 1.5) +
  scale_fill_viridis(end =0.5, discrete = T) +
  scale_color_viridis(end =0.5, discrete=T) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text=element_text(size=11),
        axis.title=element_text(size=15),
        axis.text = element_text(size=14)) +
  theme(plot.margin = margin(1.5, 0.75, 0.5, 0.5, "cm"))


# Create Supplementary Figure S1 Panel B: Comparing the ratio of inter- to intra-layer density in the observed 
# and shuffled networks:
ratio.real.pot.inter.intra <- ggplot(real.poten.edges.sh, aes(x=ratio.real))+
  geom_histogram(color="black", fill="white") +
  labs(x = "Ratio Inter:Intra-Layer Density", y="Edge Ratio (Count)") + 
  geom_vline(data=real.poten.edges.obs, aes(xintercept=ratio.real, color="Observed"),
             linetype="dashed", size = 1.5) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text=element_text(size=14),
        axis.title=element_text(size=15),
        axis.text = element_text(size=14)) +
  theme(plot.margin = margin(1.5, 0.75, 0.5, 0.5, "cm"))


# Combine two panels into Supplementary Figure S1 with ggarrange:
ggarrange(per.intra.inter.real,ratio.real.pot.inter.intra,  
          labels = c("A.", "B."),
          ncol = 2, nrow = 1,font.label = list(size = 19))

# Save in desired format
ggsave(filename = "Supp.full.net.edg.pdf",  width = 12, height = 5, units = "in")

ggsave(filename="Supp.full.net.edg.png",dpi = 900, width = 12, height = 5, units = "in")


# Calculate p-values by counting the number of times the values from the shuffled networks
# are greater or less than the observed value:
p.shuff.real.poten <- real.poten.edges.sh %>%
  mutate(ratio.great = ifelse(ratio.real > real.poten.edges.obs$ratio.real, 1, 0), # Ratio
         intra.real.great = ifelse(per.intra.real > real.poten.edges.obs$per.intra.real, 1, 0), # Intra-layer density
         inter.real.great = ifelse(per.inter.real > real.poten.edges.obs$per.inter.real, 1, 0)) # Inter-layer density

length(which(p.shuff.real.poten$ratio.great == 1)) /1000 # p < 0.001
length(which(p.shuff.real.poten$intra.real.great == 1))/1000 # p < 0.001
length(which(p.shuff.real.poten$inter.real.great == 1))/1000 # p < 0.001
######################################################################################################################################