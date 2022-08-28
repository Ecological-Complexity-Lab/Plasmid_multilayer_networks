####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 10: Analyze the number of shared modules between cows, including comparison to shuffled networks.
#
# The following figures are created:
# Figure 2A, module barplot
# Figure 2B, shared modules heatmap
# Figure S2, shared plasmids heatmap
# 
# Script tested for R version 4.1.1
####################################################################################################################
# Set working directory to wherever your files are located

# Load necessary packages:
library(tidyverse)
library(vegan)
library(ggplot2)
library(reshape)
library(viridis)
library(corrplot)
library(ggpubr)

# Starting files:
# Infomap file created in script 08_Infomap_analysis_full_network.R
load("plas_mods.df.Rda")

# Plasmid-layer (state node) dataframe created in script 03_Network_setup.R
load("plas.2k.Rda")

# Infomap files for shuffled networks created in script 09_Network_structure_randomizations
load("plas_mods.lf.df.Rda")

# Summary statistics for full network created in script 09_Network_structure_randomizations
load("plas_mods.stats.overall.Rda")
####################################################################################################################


####################################################################################################################
# Section 1: BASIC STATISTICS ON MODULE DISTRIBUTION
####################################################################################################################
# Count number of unique modules per cow (layer)
mods.per.cow <- plas_mods.df %>%
  group_by(layer_id) %>%
  summarise(tot.mods = n_distinct(module)) %>%
  mutate(layer_id = as.factor(layer_id))

# Create Figure 2, Panel A: Bar plot of reordered from lowest number of modules to highest
mods.per.cow.bar.reord <- ggplot(data=mods.per.cow, aes(x=reorder(layer_id, tot.mods), y=tot.mods)) +
  geom_bar(stat="identity", color="black", fill="white") + 
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 13),
        axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20)) +
  theme_classic() +
  xlab("Cow") + 
  ylab("Modules")

# Save in desired format:
ggsave(mods.per.cow.bar.reord,filename = "mods.per.cow.bar.reord.png",dpi = 900, height = 8, width = 8, unit="in")
ggsave(mods.per.cow.bar.reord,filename = "mods.per.cow.bar.reord.pdf",dpi = 900, height = 8, width = 8, unit="in")

# Cows per module
cows.per.mod <- plas_mods.df %>%
  group_by(module) %>%
  summarise(tot.cows = n_distinct(layer_id))
####################################################################################################################


####################################################################################################################
# Section 2: SET UP AND VISUALIZATION OF MATRIX FOR ANALYZING SHARED MODULES BETWEEN COWS
####################################################################################################################
# Make a module layer-by-cow matrix (all modules)
cow.comp <- plas_mods.df %>%
  select(., module, layer_id) %>%
  group_by(module, layer_id) %>%
  distinct() %>%
  mutate(occurence = 1) %>%
  spread(., module, occurence, fill = 0, 
         drop = FALSE) %>%
  column_to_rownames(., var="layer_id")

# Matrix of number of shared modules between cows
cow.comp.mat <- t(as.matrix(cow.comp))

# Cross product to get the shared modules on a cow-by-cow matrix
cross.obs <- crossprod(cow.comp.mat)

# Reorder the matrix by descending number of shared modules
mat.order <- mods.per.cow %>%
  arrange(tot.mods)

mat.ord.vec <- as.vector(mat.order$layer_id)

mat.order2 <- mods.per.cow %>%
  arrange(desc(tot.mods))

mat.ord.vec2 <- as.vector(mat.order2$layer_id)

# Reorder
cross.obs[,mat.ord.vec]

# Melt the matrix and calculate the percentage of shared modules
shared.mods.perc <- cross.obs %>%
  melt.matrix() %>%
  mutate(X1=as.factor(X1)) %>%
  mutate(X2=as.factor(X2)) %>%
  left_join(., mods.per.cow, by=c("X1"="layer_id")) %>%
  mutate(perc.mod.shared = value/tot.mods*100)
  

# Create Figure 2, Panel B: Heat map of shared modules between cows  
shared.mods.perc.plot <- ggplot(data = shared.mods.perc, aes(factor(X1, level=mat.ord.vec), factor(X2, level=mat.ord.vec), fill = perc.mod.shared))+
  geom_tile() +
  scale_fill_viridis(option="inferno", direction=-1) +
  theme_minimal() + 
  xlab("Cow") + ylab("Cow") + labs(fill = "Percent Modules Shared") +
  coord_equal() +
  theme(plot.margin = margin(1, 0.75, 0.5, 0.5, "cm"))

# Save in desired format:
ggsave(shared.mods.perc.plot,filename="shared.mods.perc.plot.png", dpi = 900, width = 8, height = 6, units = "in")
ggsave(filename="shared.mods.perc.plot.pdf", width = 8, height = 6)
####################################################################################################################


####################################################################################################################
# Section 3: ANALYSIS OF SHARED MODULES BETWEEN COWS
####################################################################################################################
# Remove the diagonal (sharing within the same cow, is obviously 100%):
shared.mods.perc.no.diag <- shared.mods.perc %>%
  filter(X1 != X2)

# Examine the hub of sharing:
shared.mods.big.hub <- shared.mods.perc.no.diag %>%
  filter(X1 %in% c("8", "9", "10","11"))

# Identify cows that share 100% of their modules with another cow:
shared.mod.100 <- shared.mods.perc.no.diag %>%
  filter(perc.mod.shared == 100)

# Check for correlation mean percent of modules shared per cow and the total number of modules per cow
shared.mods2 <- shared.mods.perc.no.diag %>%
  group_by(X1) %>%
  mutate(mean.perc.shared = mean(perc.mod.shared)) %>%
  select(X1, tot.mods, mean.perc.shared) %>%
  distinct

cor.test(shared.mods2$tot.mods, shared.mods2$mean.perc.shared, method='kendall')
####################################################################################################################


####################################################################################################################
# Section 4: SET UP AND VISUALIZATION OF MATRIX FOR ANALYZING SHARED PLASMIDS BETWEEN COWS
####################################################################################################################
# Set up matrix:
cow.by.plas.mat <- plas.2k %>%
  distinct() %>%
  mutate(occurence = 1) %>%
  spread(., node_id, occurence, fill = 0, 
         drop = FALSE) %>%
  column_to_rownames(., var="layer_id")

# Create matrix of number of shared modules between cows
cow.cow.plas.mat <- t(as.matrix(cow.by.plas.mat))

# Cross product to get the shared modules on a cow-by-cow matrix
cow.cow.plas.cross <- crossprod(cow.cow.plas.mat)

# Compare percentage of shared plasmids between cows
# Plasmids per cow:
plas.per.cow <- plas.2k %>%
  group_by(layer_id) %>%
  summarise(tot.plas = n_distinct(node_id)) %>%
  mutate(layer_id = as.factor(layer_id))

# Join cross.prod to total number of plasmids per cow
shared.plas.perc <- cow.cow.plas.cross %>%
  melt.matrix() %>%
  mutate(X1=as.factor(X1)) %>%
  mutate(X2=as.factor(X2)) %>%
  left_join(., plas.per.cow, by=c("X1"="layer_id")) %>%
  mutate(perc.plas.shared = value/tot.plas*100)


# Create Supplementary Figure S2: Heatmap visualization of number of plasmids shared between cows
# Use the same order as for the shared modules
shared.plas.perc.plot <- ggplot(data = shared.plas.perc, aes(factor(X1, level=mat.ord.vec), factor(X2, level=mat.ord.vec), fill = perc.plas.shared))+
  geom_tile() +
  scale_fill_viridis(option="inferno", direction=-1) +
  theme_minimal() + 
  xlab("Cow") + ylab("Cow") + labs(fill = "Percent Plasmids Shared") +
  coord_equal()
####################################################################################################################


####################################################################################################################
# Section 5: ANALYSIS OF SHARED PLASMIDS BETWEEN COWS, COMPARISON TO MODULES
####################################################################################################################
# Again, remove the diagonal (sharing within the same cow, is obviously 100%):
shared.plas.mods <- shared.plas.perc %>%
  right_join(., shared.mods.perc, by=c("X1","X2")) %>%
  # Eliminate self-comparison, e.g. cow 1 - cow 1 (diagonal)
  mutate(X1 = as.numeric(X1),
         X2 = as.numeric(X2)) %>%
  filter(!(X1 == X2))

# Correlation test between shared modules and shared plasmids:
cor.test(shared.plas.mods$perc.mod.shared, shared.plas.mods$perc.plas.shared, method='kendall')
####################################################################################################################


####################################################################################################################
# Section 6: CALCULATE NUMBER OF SHARED MODULES IN SHUFFLED NETWORKS
####################################################################################################################
# Create a cow by module matrix for each shuffled network
shuff.cow.comp <- plas_mods.lf.df %>%
  map(~select(., module,layer_id)) %>%
  map(~group_by(., module, layer_id)) %>%
  map(~distinct(.)) %>%
  map(~mutate(., occurence = 1)) %>%
  map(~spread(., module, occurence, fill = 0, 
              drop = FALSE))  %>%
  map(~drop_na(.)) %>%
  map(~remove_rownames(.)) %>%
  map(~column_to_rownames(., var="layer_id"))

# Create a matrix of number of shared modules between cows for each shuffled network
cow.comp.mat.shuff <- lapply(shuff.cow.comp, function(x) t(as.matrix(x)))

# Cross product to get the shared modules on a cow-by-cow matrix for each shuffled network
cross.obs.shuff <- lapply(cow.comp.mat.shuff, function(x) crossprod(x)) 

# Melt it to get a dataframe of shared modules for each pair of cows in each shuffled network
pre.shared.mods.perc.shuff <- lapply(cross.obs.shuff, function(x) x %>%
                                       melt.matrix() %>%
                                       mutate(X1=as.factor(X1)) %>%
                                       mutate(X2=as.factor(X2)) ) 

# Turn the lists of shared modules for each cow pair for each shuffled network into a single data frame:
pre.shared.mods.perc.shuff2 <- Map(cbind, pre.shared.mods.perc.shuff, shuff.rep = names(pre.shared.mods.perc.shuff)) %>%
  bind_rows()

# Number of modules per layer per shuffled network:
mods.per.cow.shuff <- lapply(plas_mods.lf.df, function(x) x %>%
                               group_by(layer_id) %>%
                               summarise(tot.mods = n_distinct(module)) %>%
                               mutate(layer_id = as.factor(layer_id)))

# Bind together the shared modules per cow pair per shuffled network and the total number of modules per cow:
mods.per.cow.shuff2 <- Map(cbind, mods.per.cow.shuff, shuff.rep = names(mods.per.cow.shuff)) %>%
  bind_rows()


# Calculate the shared percent of modules for each cow pair in each shuffled network 
shared.mods.perc.shuff <- pre.shared.mods.perc.shuff2 %>%
  left_join(., mods.per.cow.shuff2, by=c("shuff.rep"= "shuff.rep", "X1"="layer_id")) %>%
  mutate(perc.mod.shared = value/tot.mods*100)
# COOL, I DID IT :D
####################################################################################################################


####################################################################################################################
# Section 7: COMPARE NUMBER OF SHARED MODULES BETWEEN COWS IN OBSERVED VS SHUFFLED NETWORKS
####################################################################################################################
# Calculate p-values for each pair using Z-scores
z.sc.perc.shared <- shared.mods.perc.shuff %>%
  filter(X1 != X2) %>%
  group_by(X1, X2) %>%
  mutate(mean.mods.shared = mean(perc.mod.shared),
         sd.mods.shared = sd(perc.mod.shared),
         X1 = as.numeric(X1),
         X2 = as.numeric(X2)) %>%
  select(X1, X2, mean.mods.shared,sd.mods.shared) %>%
  distinct() %>%
  left_join(., shared.plas.mods, by=c("X1","X2")) %>%
  rowwise() %>%
  mutate(z.sc = (perc.mod.shared - mean.mods.shared)/sd.mods.shared,
         abs.z = abs(z.sc),
         p.val = abs(2*pnorm(q=abs.z, lower.tail=F)))

# Count how many pairs of cows share significantly more modules than random:
sign.greater <- z.sc.perc.shared %>%
  filter(z.sc >= 1.96)

# Count how many pairs of cows share significantly less modules than random:
sign.less <- z.sc.perc.shared %>%
  filter(z.sc <= -1.96)

# Count how many pairs of cows share the "same" (non-significant) number of modules as random
not.sig <- z.sc.perc.shared %>%
  filter(z.sc > -1.96 & z.sc < 1.96)
####################################################################################################################