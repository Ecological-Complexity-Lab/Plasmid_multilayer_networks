####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 2: Sensitivity analysis of plasmid and alignment length thresholds
# 
# 
# The following outputs are used for downstream analysis: 
# thresh2000.a.Rda
# 
# The following figures are created:
# Supplementary Figure S5A S5B, threshold effects
# 
# Script tested for R version 4.1.1
####################################################################################################################


####################################################################################################################
# SCRIPT SET-UP
####################################################################################################################
# Set working directory to wherever your files are located

# Load necessary packages
library(tidyverse)
library(Hmisc)
library(ggpubr)

# Starting files:
# Network data created in script 01_Initial_Data_Processing.R
load("full.net.dat.Rda")
####################################################################################################################


####################################################################################################################
# SECTION 1: FILTER FULL NETWORK BY PLASMID LENGTH AND ALIGNMENT LENGTH 
####################################################################################################################
# Filter the full.net.dat dataframe by plasmid length in intervals of 500 bp 
# and alignment lengths (20% of plasmid alignment or 20% of the threshold)

# 500 bp plasmid length and alignments of 20% of the shortest plasmid in a pair
thresh500.a <- full.net.dat %>%
    # Filter by plasmid length
    filter(plasmid.length.query >= 500 & plasmid.length.subject >= 500) %>%
    # Calculate the proportion of the plasmid length covered by the alignment for both 
    # plasmid subject and plasmid query:
    mutate(align.proportion.query = length/plasmid.length.query,
           align.proportion.sub = length/plasmid.length.subject ) %>%
    rowwise() %>%
    # Calculate the maximum proportion of plasmid length covered by the alignment
    # E.g. the alignment length / length of smaller plasmid
    mutate(max.align.prop = max(align.proportion.query, align.proportion.sub),
           min.align.prop = min(align.proportion.query, align.proportion.sub)) %>%
    # Filter the alignments by only including those above the defined thresholds. Basically filters
    # out whatever alignments we decide are "too short"
    filter(min.align.prop >= 0.2) %>%
    # Calculate the minimum proportion that an alignment covers
    # E.g. the alignment length / length of the longer plasmid
    # Will be used to calculate edge-weights.
    select(uniq_pair,layer_from, node_from, layer_to,node_to, max.align.prop, min.align.prop, everything()) %>%
    ungroup()
  # 28,606 alignments

# 500 bp filter but alignments must cover 20% of 500 (100 bp)
thresh500.b <- full.net.dat %>%
  # Filter by plasmid length
  filter(plasmid.length.query >= 500 & plasmid.length.subject >= 500) %>%
  # Calculate the proportion of the plasmid length covered by the alignment
  # for both plasmid subject and plasmid query:
  mutate(align.proportion.query = length/plasmid.length.query,
         align.proportion.sub = length/plasmid.length.subject ) %>%
  rowwise() %>%
  # Calculate the maximum proportion of plasmid length covered by the alignment
  # E.g. the alignment length / length of smaller plasmid
  mutate(max.align.prop = max(align.proportion.query, align.proportion.sub),
         min.align.prop = min(align.proportion.query, align.proportion.sub)) %>%
  # Filter the alignments by only including those above the defined thresholds. Basically filters
  # out whatever alignments we decide are "too short"
  filter(length >= 0.2*500) %>%
  # Calculate the minimum proportion that an alignment covers
  # E.g. the alignment length / length of the longer plasmid
  # Will be used to calculate edge-weights.
  select(uniq_pair,layer_from, node_from, layer_to,node_to, max.align.prop, min.align.prop, everything()) %>%
  ungroup()
# 40,231 alignments

# 1000 bp plasmid length and alignments of 20% of the shortest plasmid in a pair
thresh1000.a <- full.net.dat %>%
  # Filter by plasmid length
  filter(plasmid.length.query >= 1000 & plasmid.length.subject >= 1000) %>%
  # Calculate the proportion of the plasmid length covered by the alignment
  # for both plasmid subject and plasmid query:
  mutate(align.proportion.query = length/plasmid.length.query,
         align.proportion.sub = length/plasmid.length.subject ) %>%
  rowwise() %>%
  # Calculate the maximum proportion of plasmid length covered by the alignment
  # E.g. the alignment length / length of smaller plasmid
  mutate(max.align.prop = max(align.proportion.query, align.proportion.sub),
         min.align.prop = min(align.proportion.query, align.proportion.sub)) %>%
  # Filter the alignments by only including those above the defined thresholds. Basically filters
  # out whatever alignments we decide are "too short"
  filter(min.align.prop >= 0.2) %>%
  # Calculate the minimum proportion that an alignment covers
  # E.g. the alignment length / length of the longer plasmid
  # Will be used to calculate edge-weights.
  select(uniq_pair,layer_from, node_from, layer_to,node_to, max.align.prop, min.align.prop, everything()) %>%
  ungroup()
# 12,367 alignments

# 1000 bp filter but alignments must cover 20% of 1000 (200 bp)
thresh1000.b <- full.net.dat %>%
  # Filter by plasmid length
  filter(plasmid.length.query >= 1000 & plasmid.length.subject >= 1000) %>%
  # Calculate the proportion of the plasmid length covered by the alignment
  # for both plasmid subject and plasmid query:
  mutate(align.proportion.query = length/plasmid.length.query,
         align.proportion.sub = length/plasmid.length.subject ) %>%
  rowwise() %>%
  # Calculate the maximum proportion of plasmid length covered by the alignment
  # E.g. the alignment length / length of smaller plasmid
  mutate(max.align.prop = max(align.proportion.query, align.proportion.sub),
         min.align.prop = min(align.proportion.query, align.proportion.sub)) %>%
  # Filter the alignments by only including those above the defined thresholds. Basically filters
  # out whatever alignments we decide are "too short"
  filter(length >= 0.2*1000) %>%
  # Calculate the minimum proportion that an alignment covers
  # E.g. the alignment length / length of the longer plasmid
  # Will be used to calculate edge-weights.
  select(uniq_pair,layer_from, node_from, layer_to,node_to, max.align.prop, min.align.prop, everything()) %>%
  ungroup()
# 14,856 alignments

# 1500 bp plasmid length and alignments of 20% of the shortest plasmid in a pair
thresh1500.a <- full.net.dat %>%
  # Filter by plasmid length
  filter(plasmid.length.query >= 1500 & plasmid.length.subject >= 1500) %>%
  # Calculate the proportion of the plasmid length covered by the alignment
  # for both plasmid subject and plasmid query:
  mutate(align.proportion.query = length/plasmid.length.query,
         align.proportion.sub = length/plasmid.length.subject ) %>%
  rowwise() %>%
  # Calculate the maximum proportion of plasmid length covered by the alignment
  # E.g. the alignment length / length of smaller plasmid
  mutate(max.align.prop = max(align.proportion.query, align.proportion.sub),
         min.align.prop = min(align.proportion.query, align.proportion.sub)) %>%
  # Filter the alignments by only including those above the defined thresholds. Basically filters
  # out whatever alignments we decide are "too short"
  filter(min.align.prop >= 0.2) %>%
  # Calculate the minimum proportion that an alignment covers
  # E.g. the alignment length / length of the longer plasmid
  # Will be used to calculate edge-weights.
  select(uniq_pair,layer_from, node_from, layer_to,node_to, max.align.prop, min.align.prop, everything()) %>%
  ungroup()
#  7,153 alignments

# 1500 bp filter but alignments must cover 20% of 1500 (300 bp)
thresh1500.b <- full.net.dat %>%
  # Filter by plasmid length
  filter(plasmid.length.query >= 1500 & plasmid.length.subject >= 1500) %>%
  # Calculate the proportion of the plasmid length covered by the alignment
  # for both plasmid subject and plasmid query:
  mutate(align.proportion.query = length/plasmid.length.query,
         align.proportion.sub = length/plasmid.length.subject ) %>%
  rowwise() %>%
  # Calculate the maximum proportion of plasmid length covered by the alignment
  # E.g. the alignment length / length of smaller plasmid
  mutate(max.align.prop = max(align.proportion.query, align.proportion.sub),
         min.align.prop = min(align.proportion.query, align.proportion.sub)) %>%
  # Filter the alignments by only including those above the defined thresholds. Basically filters
  # out whatever alignments we decide are "too short"
  filter(length >= 0.2*1500) %>%
  # Calculate the minimum proportion that an alignment covers
  # E.g. the alignment length / length of the longer plasmid
  # Will be used to calculate edge-weights.
  select(uniq_pair,layer_from, node_from, layer_to,node_to, max.align.prop, min.align.prop, everything()) %>%
  ungroup()
# 8,015 alignments

# 2000 bp plasmid length and alignments of 20% of the shortest plasmid in a pair
thresh2000.a <- full.net.dat %>%
  # Filter by plasmid length
  filter(plasmid.length.query >= 2000 & plasmid.length.subject >= 2000) %>%
  # Calculate the proportion of the plasmid length covered by the alignment
  # for both plasmid subject and plasmid query:
  mutate(align.proportion.query = length/plasmid.length.query,
         align.proportion.sub = length/plasmid.length.subject ) %>%
  rowwise() %>%
  # Calculate the maximum proportion of plasmid length covered by the alignment
  # E.g. the alignment length / length of smaller plasmid
  mutate(max.align.prop = max(align.proportion.query, align.proportion.sub),
         min.align.prop = min(align.proportion.query, align.proportion.sub)) %>%
  # Filter the alignments by only including those above the defined thresholds. Basically filters
  # out whatever alignments we decide are "too short"
  filter(min.align.prop >= 0.2) %>%
  # Calculate the minimum proportion that an alignment covers
  # E.g. the alignment length / length of the longer plasmid
  # Will be used to calculate edge-weights.
  select(uniq_pair,layer_from, node_from, layer_to,node_to, max.align.prop, min.align.prop, everything()) %>%
  ungroup()
# 2,844 alignments


# 2000 bp filter but alignments must cover 20% of 2000 (400 bp)
thresh2000.b <- full.net.dat %>%
  # Filter by plasmid length
  filter(plasmid.length.query >= 2000 & plasmid.length.subject >= 2000) %>%
  # Calculate the proportion of the plasmid length covered by the alignment
  # for both plasmid subject and plasmid query:
  mutate(align.proportion.query = length/plasmid.length.query,
         align.proportion.sub = length/plasmid.length.subject ) %>%
  rowwise() %>%
  # Calculate the maximum proportion of plasmid length covered by the alignment
  # E.g. the alignment length / length of smaller plasmid
  mutate(max.align.prop = max(align.proportion.query, align.proportion.sub),
         min.align.prop = min(align.proportion.query, align.proportion.sub)) %>%
  # Filter the alignments by only including those above the defined thresholds. Basically filters
  # out whatever alignments we decide are "too short"
  filter(length >= 0.2*2000) %>%
  # Calculate the minimum proportion that an alignment covers
  # E.g. the alignment length / length of the longer plasmid
  # Will be used to calculate edge-weights.
  select(uniq_pair,layer_from, node_from, layer_to,node_to, max.align.prop, min.align.prop, everything()) %>%
  ungroup()
# 3,255 alignments

# 2500 bp plasmid length and alignments of 20% of the shortest plasmid in a pair
thresh2500.a <- full.net.dat %>%
  # Filter by plasmid length
  filter(plasmid.length.query >= 2500 & plasmid.length.subject >= 2500) %>%
  # Calculate the proportion of the plasmid length covered by the alignment
  # for both plasmid subject and plasmid query:
  mutate(align.proportion.query = length/plasmid.length.query,
         align.proportion.sub = length/plasmid.length.subject ) %>%
  rowwise() %>%
  # Calculate the maximum proportion of plasmid length covered by the alignment
  # E.g. the alignment length / length of smaller plasmid
  mutate(max.align.prop = max(align.proportion.query, align.proportion.sub),
         min.align.prop = min(align.proportion.query, align.proportion.sub)) %>%
  # Filter the alignments by only including those above the defined thresholds. Basically filters
  # out whatever alignments we decide are "too short"
  filter(min.align.prop >= 0.2) %>%
  # Calculate the minimum proportion that an alignment covers
  # E.g. the alignment length / length of the longer plasmid
  # Will be used to calculate edge-weights.
  select(uniq_pair,layer_from, node_from, layer_to,node_to, max.align.prop, min.align.prop, everything()) %>%
  ungroup()
# 1,446 alignments

# 2500 bp filter but alignments must cover 20% of 2500 (500 bp)
thresh2500.b <- full.net.dat %>%
  # Filter by plasmid length
  filter(plasmid.length.query >= 2500 & plasmid.length.subject >= 2500) %>%
  # Calculate the proportion of the plasmid length covered by the alignment
  # for both plasmid subject and plasmid query:
  mutate(align.proportion.query = length/plasmid.length.query,
         align.proportion.sub = length/plasmid.length.subject ) %>%
  rowwise() %>%
  # Calculate the maximum proportion of plasmid length covered by the alignment
  # E.g. the alignment length / length of smaller plasmid
  mutate(max.align.prop = max(align.proportion.query, align.proportion.sub),
         min.align.prop = min(align.proportion.query, align.proportion.sub)) %>%
  # Filter the alignments by only including those above the defined thresholds. Basically filters
  # out whatever alignments we decide are "too short"
  filter(length >= 0.2*2500) %>%
  # Calculate the minimum proportion that an alignment covers
  # E.g. the alignment length / length of the longer plasmid
  # Will be used to calculate edge-weights.
  select(uniq_pair,layer_from, node_from, layer_to,node_to, max.align.prop, min.align.prop, everything()) %>%
  ungroup()
# 1,729 alignments
  
# 3000 bp plasmid length and alignments of 20% of the shortest plasmid in a pair
thresh3000.a <- full.net.dat %>%
  # Filter by plasmid length
  filter(plasmid.length.query >= 3000 & plasmid.length.subject >= 3000) %>%
  # Calculate the proportion of the plasmid length covered by the alignment
  # for both plasmid subject and plasmid query:
  mutate(align.proportion.query = length/plasmid.length.query,
         align.proportion.sub = length/plasmid.length.subject ) %>%
  rowwise() %>%
  # Calculate the maximum proportion of plasmid length covered by the alignment
  # E.g. the alignment length / length of smaller plasmid
  mutate(max.align.prop = max(align.proportion.query, align.proportion.sub),
         min.align.prop = min(align.proportion.query, align.proportion.sub)) %>%
  # Filter the alignments by only including those above the defined thresholds. Basically filters
  # out whatever alignments we decide are "too short"
  filter(min.align.prop >= 0.2) %>%
  # Calculate the minimum proportion that an alignment covers
  # E.g. the alignment length / length of the longer plasmid
  # Will be used to calculate edge-weights.
  select(uniq_pair,layer_from, node_from, layer_to,node_to, max.align.prop, min.align.prop, everything()) %>%
  ungroup()
# 710 alignments

# 3000 bp filter but alignments must cover 20% of 3000 (600 bp)
thresh3000.b <- full.net.dat %>%
  # Filter by plasmid length
  filter(plasmid.length.query >= 3000 & plasmid.length.subject >= 3000) %>%
  # Calculate the proportion of the plasmid length covered by the alignment
  # for both plasmid subject and plasmid query:
  mutate(align.proportion.query = length/plasmid.length.query,
         align.proportion.sub = length/plasmid.length.subject ) %>%
  rowwise() %>%
  # Calculate the maximum proportion of plasmid length covered by the alignment
  # E.g. the alignment length / length of smaller plasmid
  mutate(max.align.prop = max(align.proportion.query, align.proportion.sub),
         min.align.prop = min(align.proportion.query, align.proportion.sub)) %>%
  # Filter the alignments by only including those above the defined thresholds. Basically filters
  # out whatever alignments we decide are "too short"
  filter(length >= 0.2*3000) %>%
  # Calculate the minimum proportion that an alignment covers
  # E.g. the alignment length / length of the longer plasmid
  # Will be used to calculate edge-weights.
  select(uniq_pair,layer_from, node_from, layer_to,node_to, max.align.prop, min.align.prop, everything()) %>%
  ungroup()
# 894 alignments
####################################################################################################################

####################################################################################################################
# SECTION 2: CALCULATE THE EFFECTS OF EACH THRESHOLD ON NUMBER PLASMIDS AND ALIGNMENTS RETAINED IN DATA-SET
####################################################################################################################
# First make a list of the dataframes, using llist from Hmisc package preserving the names of each dataframe:
thresh.list <- llist(thresh500.a, thresh500.b,thresh1000.a, thresh1000.b,
                     thresh1500.a, thresh1500.b, thresh2000.a, thresh2000.b,
                     thresh2500.a, thresh2500.b, thresh3000.a, thresh3000.b)

# Specify names in list
thresh.df.names <- as.data.frame(names(thresh.list))

# Calculate the total number of alignments (nrows) per threshold
thresh.tot.align <- as.data.frame(sapply(thresh.list, nrow))

# Function to calculate the number of plasmids retained per threshold:
num.plas.fun <- function(x) {
  x %>%
    select(starts_with("node")) %>%
    gather(.,key=node_type, value=node_id) %>%
    select(-node_type) %>%
    distinct() %>%
    nrow()
}

# Apply the function:
thresh.num.plas <- as.data.frame(sapply(thresh.list, num.plas.fun))


# Functions to count the number of intra vs. inter-layer edges at each threshold level:
# Intra first:
intra.count <- function(x) { 
  x %>%
    select(layer_from, node_from, layer_to, node_to, edge_type) %>%
    distinct() %>%
    filter(edge_type == "Intra") %>%
    nrow()
}

# Apply the function:
thresh.intra <- as.data.frame(sapply(thresh.list, intra.count))

# Repeat for inter:
inter.count <- function(x) { 
  x %>%
    select(layer_from, node_from, layer_to, node_to, edge_type) %>%
    distinct() %>%
    filter(edge_type == "Inter") %>%
    nrow()
}

# Apply the function:
thresh.inter <- as.data.frame(sapply(thresh.list, inter.count))

# Put the information for each threshold level into a single table
thresh.df <- bind_cols(thresh.num.plas,thresh.tot.align, thresh.intra, thresh.inter)

# Set column names:
thresh.col.names <- as.vector(c("num_plasmids","total_alignments",
                      "intra","inter"))

# Apply column names to the dataframe:
thresh.df <- thresh.df %>% 
  rename_with(~thresh.col.names)

# Minor edits/format to dataframe:
thresh.df1 <- thresh.df %>%
  rownames_to_column(., var = "thresh") %>%
  mutate(min.leng = as.numeric(str_extract(thresh, "\\d+"))) %>%
  mutate(align.type = str_sub(thresh,-1),
         align.type = ifelse(align.type == "a", "20% shorter plasmid", "20% threshold"))

# Gather for plotting:
thresh.df.intra.inter <- thresh.df1 %>%
  dplyr::rename(Inter = inter, Intra = intra, Total = total_alignments) %>%
  pivot_longer(cols=c("Total", "Inter", "Intra"), names_to = "edge.type", values_to = "alignment.count")
####################################################################################################################


####################################################################################################################
# SECTION 3: VISUALIZATION OF THRESHOLD EFFECTS 
####################################################################################################################
# Create Supplementary Figure S5:
# Panel A:
plot.thresh <- ggplot(data = thresh.df1) +
  geom_line(aes(x=min.leng, y=num_plasmids, color=align.type), size=1) +
  labs(x = "Minimum Plasmid Length", y="Plasmids (Count)") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text=element_text(size=11),
        axis.title=element_text(size=16),
        axis.text = element_text(size=13)) +
  theme(plot.margin = margin(1.5, 0.75, 0.5, 0.5, "cm"))

# Can save intermediate plot if desired
#ggsave(plot.thresh ,filename="thresh.plas.count.png", dpi = 900, width = 10, height = 4, units = "in")

# Panel B:
plot.thresh.intra.inter <- ggplot(data=thresh.df.intra.inter) +
  geom_line(aes(x=min.leng, y=alignment.count, color=align.type, linetype=edge.type),size=1.25) +
  scale_linetype_manual(values=c("solid","twodash", "dotted")) +
  labs(x = "Minimum Plasmid Length", y="Alignments (Count)") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text=element_text(size=11),
        axis.title=element_text(size=16),
        axis.text = element_text(size=13)) +
  theme(plot.margin = margin(1.5, 0.75, 0.5, 0.5, "cm"))

# Can save intermediate plot if desired
#ggsave(plot.thresh.intra.inter,filename="length.thresh.intra.inter.png", dpi = 900, width = 10, height = 4, units = "in")
  
# Combine plots with ggarrange:
thresh.panel.plot <- ggarrange(plot.thresh, plot.thresh.intra.inter,
                                        labels = c("A.", "B."),
                                        ncol = 1, nrow = 2,
                                        common.legend = FALSE)

# Save Figure S5 in desired format: 
ggsave(thresh.panel.plot,filename="thresh.panel.png", dpi = 900, width = 10, height = 8, units = "in")

ggsave(thresh.panel.plot,filename="thresh.panel.pdf", dpi = 900, width = 10, height = 8, units = "in")
####################################################################################################################


####################################################################################################################
# SECTION 4: SAVE THE NETWORK DATA  
####################################################################################################################
# Save the data at the chosen threshold for further research:
# 2000 bp length, alignment coverage of 20% of the shortest plasmid in the pair:
save(thresh2000.a, file="thresh2000.a.Rda")
####################################################################################################################