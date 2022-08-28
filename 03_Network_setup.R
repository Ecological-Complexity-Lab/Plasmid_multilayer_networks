####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 3: Network set-up
# 
# 
# The following outputs are used for downstream analysis: 
# net.dat2k.ew.Rda, plas.2k.Rda, plasmid.2k.metadat.Rda, plas.2k.name.node.id.Rda 
# 
# Script tested for R version 4.1.1
####################################################################################################################


####################################################################################################################
# SCRIPT SET-UP
####################################################################################################################
# Set working directory to wherever your files are located

# Load necessary packages:
library(tidyverse)
library(seqinr)

# Starting files:
# Plasmid metadata file created in script 01_Initial_data_processing.R
load("plasmid.metadat.updated.Rda")

# Plasmid similarity dataframe filtered by plasmid and alignment length created in script 02_Threshold_sensitivity.R
load("thresh2000.a.Rda")
####################################################################################################################


####################################################################################################################
# Section 1: CALCULATE EDGE-WEIGHTS
####################################################################################################################
#Before creating the edgeweights, which requires reducing all pairs to one alignment, calculate the number of 
# alignments between each pair
num.aligns.per.pair <- thresh2000.a %>%
  group_by(layer_from, node_from, layer_to, node_to) %>%
  summarise(num.als = n())

# Calculate the edge-weights
# Keep only relevant columns 
net.dat2k.ew <- thresh2000.a  %>%
  rowwise() %>%
  mutate(preweight = min.align.prop*(pident/100)) %>%
  group_by(uniq_pair) %>%   
  mutate(weight = sum(preweight)) %>%
  select(-min.align.prop, -max.align.prop, -qstart, -qend, -sstart, -send, -evalue,-qs.qe.diff,
         -align.proportion.query,-align.proportion.sub) %>%
  select(-preweight) %>%
  distinct(., uniq_pair, weight, .keep_all=T)

# Save this data:
save(net.dat2k.ew, file="net.dat2k.ew.Rda")

# Count edge types:
layer.type.count <- net.dat2k.ew %>%
  group_by(edge_type) %>%
  count()
####################################################################################################################


####################################################################################################################
# Section 2: CREATE METADATA AND LIST FILES FOR NETWORK PLASMIDS
####################################################################################################################
# Node-from plasmids
plas.2k.fr <- net.dat2k.ew %>%
  ungroup() %>%
  select(layer_from, node_from, weight, edge_type) %>%
  mutate(layer_id = layer_from,
         node_id = node_from) %>%
  select(-layer_from, -node_from)

# Node-to plasmids
plas.2k.to <- net.dat2k.ew %>%
  ungroup() %>%
  select(layer_to, node_to, weight, edge_type) %>%
  mutate(layer_id = layer_to,
         node_id = node_to) %>%
  select(-layer_to, -node_to) 

# All plasmids with their layers
plas.2k <- bind_rows(plas.2k.fr, plas.2k.to) %>%
  select(layer_id, node_id) %>%
  distinct()

# Save this data
save(plas.2k, file="plas.2k.Rda")

# Create metadata for plasmids in network (>2k)
plasmid.2k.metadat <- plasmid.metadat.updated %>%
 filter(node_id %in% plas.2k$node_id) %>%
 mutate(node_id = as.integer(node_id)) %>%
 droplevels()

# Save this data:
save(plasmid.2k.metadat, file = "plasmid.2k.metadat.Rda")

# Link only plasmid node id's to plasmid names:
plas.2k.name.node.id <- plasmid.2k.metadat %>%
  select(plasmid_name, node_id, plasmid_length) %>%
  filter(!is.na(.)) %>%
  distinct()

# Save this data:
save(plas.2k.name.node.id, file="plas.2k.name.node.id.Rda")

# List of node ID's (only) in final data set:
node.ids <- plas.2k %>%
  select(node_id) %>%
  distinct() %>%
  droplevels()
####################################################################################################################


####################################################################################################################
# Section 3: CREATE NEW FASTA FILE FOR 2K PLASMID DATA SET
####################################################################################################################

# Outwrite new fasta file including only the plasmids included in the network:
# Filter full fasta file down to sequences that have passed the filters

library(seqinr)
# Read the fasta file containing the plasmid sequences (created in 01_Initial_data_processing.R)
fast1 <- read.fasta("plasmid.idnums.fasta")

# List of filtered plasmids (created in section above, but load if necessary)
node.ids

# Filter the full fasta (with ID numbers) to 
fasta.filtered <- fast1[c(which(names(fast1) %in% node.ids$node_id))]

View(fasta.filtered)

# Outwrite new fasta
write.fasta(fasta.filtered,names = names(fasta.filtered), 
            file.out = "plasmids.2k.fasta")
####################################################################################################################