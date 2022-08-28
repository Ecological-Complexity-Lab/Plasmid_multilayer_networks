####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 9: Randomizaton of layer identity (state nodes) to determine if the multilayer
#            network structure and characteristics are non-random
# 
# The following outputs are used for downstream analysis: 
# plas_mods.lf.df.Rda, plas_mods.lf.stats.Rda, shuff.lf.df.Rda, 
# shuff.lf.df2.Rda, plas_mods.stats.overall.Rda
#
# The following figures are created:
# Figure 1C, comparison of modules between observed and shuffled networks
# Figure 1D, comparison of flow in largest modules between observed and shuffled networks
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
library(ggpubr)
library(viridis)

# Starting files:
# Shuffled network created in script 06_Network_shuffling
load("shuff.lf.list.Rda")

# Node id table created in script 08_Infomap_analysis_full_network.R
load("plas.2k.list.Rda")

# Infomap files created in script 08_Infomap_analysis_full_network.R
load("plas_mods.df.Rda")
load("plas_mods.stats.Rda")
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


####################################################################################################################
# Section 2: RUN INFOMAP ON SHUFFLED NETWORKS
####################################################################################################################
# Making multilayer network on each of the shuffled networks in the list with infomapecology package
plas.multilay.lf.list <- lapply(shuff.lf.list,function(x) create_multilayer_object(extended = x, 
                                                                                   nodes = plas.2k.list,  layers = layer.id.table.2k))

#Run infomap on each of the shuffled networks
plas_modules.lf <- lapply(plas.multilay.lf.list, function(x) run_infomap_multilayer(M=x, relax = F, 
                                                                                          flow_model = 'undirected', silent = T, 
                                                                                          trials = 100, seed = 123, temporal_network = F))

# Save the data:
save(plas_modules.lf, file="plas_modules.lf.Rda")
####################################################################################################################


####################################################################################################################
# Section 3: STATISTICAL ANALYSIS ON INFOMAP MODULES
####################################################################################################################
# Turn the infomap outputs into a dataframe and reorganize slightly:
plas_mods.lf.df <- lapply(plas_modules.lf, function (x) as.data.frame(x$modules) %>%
  select(module, everything()))

# Save the data
save(plas_mods.lf.df, file="plas_mods.lf.df.Rda")

# Calculate the number of nodes and layers and the flow per module:
plas_mods.lf.stats <- lapply(plas_mods.lf.df, function(x) x %>%
  group_by(module) %>%
  summarise(n.state.nodes=n(),
            n.phys.nodes = n_distinct(node_id),
            n.layers = n_distinct(layer_id),
            mod.flow=sum(flow)) )

# Save the data
save(plas_mods.lf.stats, file="plas_mods.lf.stats.Rda")

# Create a new column that ID's the shuffling repetition number of each iteration
# Bind each list together into a single dataframe. 
shuff.lf.df <- Map(cbind, plas_mods.lf.stats, shuff.rep = names(shuff.lf.list)) %>%
  bind_rows()

# Save the data
save(shuff.lf.df, file = "shuff.lf.df.Rda")

# Summarize statistics per shuffled network
shuff.lf.df2 <- shuff.lf.df %>%
  group_by(shuff.rep) %>%
  summarise(n.mods = n_distinct(module),
         n.phys.nodes = mean(n.phys.nodes),
         n.state.nodes = mean(n.state.nodes), 
         n.layers = mean(n.layers),
         n.mod.flow = mean(mod.flow)) 

# Save the data
save(shuff.lf.df2, file="shuff.lf.df2.Rda")
####################################################################################################################


####################################################################################################################
# Section 4: COMPARISONS OF MEAN CHARACTERISTICS PER MODULE IN OBSERVED VS SHUFFLED NETWORKS
####################################################################################################################
# Plot number of modules, average number of state nodes, physical nodes, and layers per module 
# against the observed values

# Calculate measures for the entire observed network for comparison to shuffled networks
plas_mods.stats.overall <- plas_mods.stats %>%
  summarise(n.mods = n_distinct(module),
            n.phys.nodes = mean(n.phys.nodes),
            n.state.nodes = mean(n.state.nodes), 
            n.layers = mean(n.layers),
            n.flow = mean(mod.flow))

# Save the data
save(plas_mods.stats.overall, file= "plas_mods.stats.overall.Rda")

# Plot Figure 1, Panel C: Number of modules in observed vs shuffled networks:
num.mods <- ggplot(shuff.lf.df2, aes(x=n.mods))+
  geom_histogram(color="black", fill="white") +
  labs(x = "Number of Modules", y="Shuffled Networks") + 
  geom_vline(data=plas_mods.stats.overall, aes(xintercept=n.mods, color="Observed"),
         linetype="dashed", size = 1.5) + 
  labs(color="") +
  theme(legend.position = "none") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20)) 

# Save in desired format:
ggsave(filename="num.mods.png", dpi = 1800, width = 8, height = 6, units = "in")
ggsave(num.mods, filename="num.mods.pdf", dpi = 900, width = 8, height = 6, units = "in")


# Compare the structure of the observed network to the shuffled networks.
# Calculate a p-value by counting the number of times the random result is different than the observed
# number of modules and mean number of state nodes, physical nodes, and layers per module:
p.shuff <- shuff.lf.df2 %>%
  mutate(mods.great = ifelse(n.mods > plas_mods.stats.overall$n.mods, 1, 0),
         phys.great = ifelse(n.phys.nodes > plas_mods.stats.overall$n.phys.nodes, 1, 0),
         st.great = ifelse(n.state.nodes > plas_mods.stats.overall$n.state.nodes, 1, 0),
         lay.great = ifelse(n.layers > plas_mods.stats.overall$n.layers, 1, 0),
         flow.great = ifelse(n.layers > plas_mods.stats.overall$n.layers, 1, 0))

length(which(p.shuff$mods.great == 1))
length(which(p.shuff$st.great == 1))
length(which(p.shuff$phys.great == 1))
length(which(p.shuff$lay.great == 1))
length(which(p.shuff$flow.great == 1))
####################################################################################################################


####################################################################################################################
# Section 5: COMPARISON OF THE SINGLE LARGEST MODULE IN OBSERVED VS SHUFFLED
####################################################################################################################
# Compare the number of physical nodes and flow of the biggest module in the observed network to physical nodes  
# and flow in the biggest module in each of the shuffled networks: 
biggest.1mod.comparison.dat <- shuff.lf.df %>%
  group_by(shuff.rep) %>%
  slice_max(., n.state.nodes, n=1, with_ties = FALSE)

biggest.1mod.obs <- plas_mods.stats %>%
  slice_max(., n.state.nodes, n=1, with_ties = FALSE)

# Compare the size and flow of the largest module in the observed network to the shuffled networks.
# Calculate a p value by counting the number of times the random result is different than the observed
# Modules, state nodes, physical nodes, and layers:
p.shuff.biggest.1mod <- biggest.1mod.comparison.dat %>%
  mutate(phys.great = ifelse(n.phys.nodes >= biggest.1mod.obs$n.phys.nodes, 1, 0),
         st.great = ifelse(n.state.nodes > biggest.1mod.obs$n.state.nodes, 1, 0),
         lay.great = ifelse(n.layers > biggest.1mod.obs$n.layers, 1, 0),
         mod.fl.great = ifelse(mod.flow > biggest.1mod.obs$mod.flow, 1, 0))

length(which(p.shuff.biggest.1mod$st.great == 1))
length(which(p.shuff.biggest.1mod$phys.great == 1)) /1000
length(which(p.shuff.biggest.1mod$lay.great == 1))
length(which(p.shuff.biggest.1mod$mod.fl.great == 1)) / 1000


# Create Figure 1, Panel D: Flow in the largest module of the observed vs. shuffled networks:
biggest.1mod.flow <- ggplot(biggest.1mod.comparison.dat, aes(x=mod.flow))+
  geom_histogram(color="black", fill="white") +
  labs(x = "Module Flow", y="Shuffled Networks") + 
  geom_vline(data=biggest.1mod.obs, aes(xintercept=mod.flow, color="Observed"),
             linetype="dashed", size = 1.5) + 
  labs(color="") +
  theme(legend.position = "none") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20)) 
ggsave(filename = "biggest.1mod.flow.png",dpi = 1800, height = 8, width = 6, unit="in")
####################################################################################################################


####################################################################################################################
# Section 6: LAYER TYPES AND FLOW IN LARGEST MODULE
####################################################################################################################
# Compare intermodule and intramodule flow in the largest module:
flow.observ <- net.2k.edgelist %>%
  left_join(., plas_mods.df, by = c("node_from" = "node_id", "layer_from" = "layer_id")) %>%
  mutate(module.from = module) %>%
  select(-module)

flow.observ2 <- flow.observ %>%
  left_join(., plas_mods.df, by = c("node_to" = "node_id", "layer_to" = "layer_id")) %>%
  mutate(module.to=module) %>%
  select(-module)

flow.observ3 <- flow.observ2 %>%
  mutate(mod.flow.type=ifelse(module.from==module.to,'Intramod','Intermod'))

flow.observ4 <- flow.observ3 %>%
  group_by(module.from, module.to, mod.flow.type) %>%
  mutate(flow.strength = sum(weight)) %>%
  select(module.from, module.to, mod.flow.type, flow.strength)

flow.observ5 <- flow.observ4 %>%
  ungroup() %>%
  distinct() 

flow.observ.intramod <- flow.observ5 %>%
  filter(mod.flow.type == "Intramod") %>%
  mutate(module=module.from) %>%
  select(module,everything(),-module.from, -module.to)

flow.observ.intermod <- flow.observ5 %>%
  filter(mod.flow.type == "Intermod") %>%
  gather(module.from:module.to, key = key, value = module) %>%
  select(module, flow.strength, mod.flow.type) %>%
  group_by(module) %>%
  mutate(flow.strength = sum(flow.strength)) %>%
  distinct() %>%
  ungroup()

flow.observ.final <- flow.observ.intramod %>%
  bind_rows(., flow.observ.intermod)

# Flow for module 1 - 155.8 intramodular flow, 2.1 intermodular flow

# Inter vs. intralayer edges in the biggest module
big.mod <- flow.observ3 %>%
  filter(module.from==1 | module.to ==1) %>%
  distinct() %>%
  mutate(edge_type=ifelse(layer_from==layer_to,'Intra','Inter'))

bigmod.edgetype.count <- big.mod %>%
  group_by(edge_type) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n))   

bigmod.mod.flow.count <- big.mod %>%
  group_by(mod.flow.type) %>%
  summarise(n = n()) %>%
  mutate(Freq = n/sum(n))
####################################################################################################################