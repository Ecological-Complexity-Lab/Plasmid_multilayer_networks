####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 12: Effect of plasmid traits (mobility and antimicrobial resistance) on network structure
# 
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
# File linking plasmid node ID's to original names, created in script 03_Network_setup.R
load("plas.2k.name.node.id.Rda")

# Degree of physical nodes, created in script 04_Basic_network_statistics.R
load("deg.str.2k.all.phys.Rda")

# Cow per plasmid created in script 04_Basic_network_statistics.R
load("cows.per.plasmid.Rda")

# Infomap data created in script 08_Infomap_analysis_full_network.R
load("plas_mods.df.Rda")

# Subnetwork data created in script 11_Subnetwork_analysis.R
load("hgt.netdat.et.Rda")
load("distdisp.netdat.et.Rda") 
load("recdisp.netdat.et.Rda")
####################################################################################################################


####################################################################################################################
# Section 1: FORMATTING THE STARTING FILES
####################################################################################################################
# Read the file with the annotations:
ann.dat1 <- read.csv("nr_annotations_min1000.csv", header=TRUE,na.strings=c("","NA"))

# Minor formatting of the file:
ann.dat2 <- ann.dat1 %>%
  # Get rid of prefix before node name
  separate(gene, c("prefix", "gene"), "gene_\\d+_\\d+_\\d+_") %>%
  select(-prefix) %>%
  dplyr::rename(plasmid_name = gene)

# Filter to keep only data from the plasmids included in our analysis (1344 plasmids)
ann.dat3 <- plas.2k.name.node.id %>% 
  left_join(., ann.dat2, by="plasmid_name") %>%
  select(node_id, plasmid_length, final) %>%
  distinct() %>%
  mutate(presence = ifelse(final==0, 0, 1))
####################################################################################################################


####################################################################################################################
# Section 2: FORMATTING THE STARTING FILES FOR MOB GENE ANALYSIS
####################################################################################################################
# Identify plasmids with mob-related genes
mob.types1 <- ann.dat3 %>%
  filter(grepl("mob|PRE|recombinase|relaxase", final)) %>%
  pivot_wider(names_from = final, values_from = presence, values_fill=0)

# Plasmids without mob-related genes:
no.mob <- ann.dat3 %>%
  filter(!(node_id %in% mob.types1$node_id)) %>%
  select(node_id, plasmid_length) %>%
  distinct() %>%
  # Create mob variables to correspond with mob data:
  mutate(PRE = 0,
         mob = 0,
         recombinase = 0,
         relaxase =0)

# Label plasmids as having mob genes or not (presence / absence)
# Use both binary and categorical variable (useful for plotting)
mob.dat.an <- mob.types1 %>%
  bind_rows(., no.mob) %>%
  # Make category for mob and conj:
  mutate(mob.pres = ifelse(mob > 0, "Present","Absent"))
####################################################################################################################


####################################################################################################################
# Section 3: JOIN MOB DATA TO RELEVANT DATA SETS
####################################################################################################################
# Join mob data to relevant data sets:
# Degree data: 
mob.deg.df <- deg.str.2k.all.phys %>%
  select(node_id, degree, strength) %>%
  ungroup() %>%
  full_join(., mob.dat.an, by="node_id")

# Plasmid-cow distribution data:
cows.per.plasmid.mob <- cows.per.plasmid %>%
  left_join(., mob.dat.an)

# Infomap data, plasmid-level:
plas_flow.mob <- plas_mods.df %>%
  select(node_id, flow) %>%
  distinct() %>%
  left_join(., mob.dat.an)

# Infomap data, module-level:
plas_mods.mob <- plas_mods.df %>%
  left_join(., mob.dat.an, by="node_id") %>%
  group_by(module) %>%
  summarise(n.state.nodes=n(),
            n.phys.nodes = n_distinct(node_id),
            n.layers = n_distinct(layer_id),
            mean.size=mean(plasmid_length),
            mob.sum = sum(mob),
            mod.flow=sum(flow)) %>%
  mutate(mob.pres.cat = ifelse(mob.sum > 0, "Present", "Absent"),
         mob.pres.num = ifelse(mob.sum > 0, 1,0),
         prop.mob = mob.sum/n.state.nodes)
####################################################################################################################


####################################################################################################################
# Section 4: COMPARE EFFECT OF MOB GENES ON PLASMID-LEVEL MEASURES
####################################################################################################################
# Test for effect of mob genes on degree centrality:
# Wilcoxon on degree:
wilcox.test(mob.deg.df$degree ~ mob.deg.df$mob, paired=F) 

# Compare the median degree for plasmids with and without mob genes
mob.deg.df%>%
  group_by(mob)%>% 
  summarise(Mean = mean(degree),Median=median(degree)) 

# Test for effect of mob genes on the number of cows a plasmid is found in
wilcox.test(cows.per.plasmid.mob$num.cows ~ cows.per.plasmid.mob$mob.pres)
# Not significant

# Compare the median number of plasmids with/without mob genes cows are found in
cows.per.plasmid.mob%>%
  group_by(mob.pres)%>% 
  summarise(Median=median(num.cows))

# Test for effect of mob genes on flow per plasmid:
wilcox.test(plas_flow.mob$flow ~ plas_flow.mob$mob, paired=F) 

# Compare median flow for plasmids with and without mob genes
plas_flow.mob%>%
  group_by(mob)%>% 
  summarise(Median=median(flow))
####################################################################################################################


####################################################################################################################
# Section 5: COMPARE EFFECT OF MOB GENES ON MODULE-LEVEL MEASURES
####################################################################################################################
# Test if mob presence influences module size and flow
#State nodes
wilcox.test(plas_mods.mob$n.state.nodes ~ plas_mods.mob$mob.pres.cat,paired=F)

# Compare median state nodes for modules with and without mob genes
plas_mods.mob%>%
  group_by(mob.pres.cat)%>% 
  summarise(Median=median(n.state.nodes)) 


# Physical nodes 
wilcox.test(plas_mods.mob$n.phys.nodes ~ plas_mods.mob$mob.pres.cat, paired=F)

# Compare median physical nodes for modules with and without mob genes
plas_mods.mob%>%
  group_by(mob.pres.cat)%>% 
  summarise(Median=median(n.phys.nodes)) 

# Layers
wilcox.test(plas_mods.mob$n.layers ~ plas_mods.mob$mob.pres.cat)

# Compare layers for modules with and without mob genes
plas_mods.mob%>%
  group_by(mob.pres.cat)%>% 
  summarise(Median=median(n.layers), Mean=mean(n.layers)) 

# Flow:
wilcox.test(plas_mods.mob$mod.flow ~ plas_mods.mob$mob.pres.cat)

# Compare median flow for modules with and without mob genes
plas_mods.mob%>%
  group_by(mob.pres.cat)%>% 
  summarise(Median=median(mod.flow))

# Median values for each characteristic:
median.mods <- plas_mods.mob%>%
  group_by(mob.pres.cat)%>% 
  summarise(median.st=median(n.state.nodes),
            median.phys=median(n.phys.nodes),
            median.layers=median(n.layers),
    median.flow=median(mod.flow))
####################################################################################################################


####################################################################################################################
# Section 6: DATA FORMATTING FOR AMR GENES
####################################################################################################################
# Identify plasmids with AMR-related genes
amr.types1 <- ann.dat3 %>%
  filter(grepl("penicillin-binding protein|beta-lactamase|tetracycline resistance|tetracycline resistance protein", final)) %>%
  pivot_wider(names_from = final, values_from = presence, values_fill=0) %>%
  mutate(tet.res = ifelse(`tetracycline resistance` > 0 | `tetracycline resistance protein` > 0, 1, 0),
         # Note confirm that "beta-lactamase" are in columns 3 and 4
         beta.lact = ifelse(.[[3]] > 0 | .[[4]] > 0, 1, 0),
         res="Present") %>%
  dplyr::rename(pen.bind = `penicillin-binding protein`) %>%
  select(-starts_with("tetracycline"), -starts_with("beta-"))

# Plasmids without AMR-related genes:
no.amr <- ann.dat3 %>%
  filter(!(node_id %in% amr.types1$node_id)) %>%
  select(node_id, plasmid_length) %>%
  distinct() %>%
  # Create mob variables to correspond with mob data:
  mutate(pen.bind = 0,
         tet.res = 0,
         beta.lact = 0,
         res ="Absent")

# Bind data with and without mob genes:
amr.dat.an <- amr.types1 %>%
  bind_rows(., no.amr)
####################################################################################################################


####################################################################################################################
# Section 7: IDENTIFY COWS AND MODULES CONTAINTING AMR PLASMIDS
####################################################################################################################
# Join AMR data to module data in order to check which cows and modules AMR plasmids are in
plas_mods.amr <- plas_mods.df %>%
  left_join(., amr.dat.an, by="node_id")

# Modules with any resistance:
# Manual check to identify cows and modules (since there are only 22 state nodes with
# any of the three resistance types)
plas_mods.amr2 <- plas_mods.amr %>%
  filter(res=="Present")
####################################################################################################################


####################################################################################################################
# Section 8: AMR PLASMID SUBNETWORKS
####################################################################################################################
# Beta lactam resistant plasmids:
beta.list <- c(108, 289, 3952, 4543, 5983, 6702, 7146, 7573)

# Recent dispersal, plasmids from 
beta.recdisp.fr <- amr.types1 %>%
  filter(beta.lact==1) %>%
  left_join(., recdisp.netdat.et, by = c("node_id"="node_from"), keep=TRUE) 

# Recent dispersal, plasmids to 
beta.recdisp.to <- amr.types1 %>%
  filter(beta.lact==1) %>%
  left_join(., recdisp.netdat.et, by = c("node_id"="node_to"), keep = TRUE)

# Join together
beta.recdisp <- beta.recdisp.fr %>%
  bind_rows(beta.recdisp.to) %>%
  distinct() %>%
  filter(node_from %in% beta.list & node_to %in% beta.list) %>%
  distinct() %>%
  select(node_from, node_to, weight, everything())

# Distant dispersal, plasmids from
beta.distdisp.fr <- amr.types1 %>%
  filter(beta.lact==1) %>%
  left_join(., distdisp.netdat.et, by = c("node_id"="node_from"), keep=TRUE) 

# Distant dispersal, plasmids to
beta.distdisp.to <- amr.types1 %>%
  filter(beta.lact==1) %>%
  left_join(., distdisp.netdat.et, by = c("node_id"="node_to"), keep = TRUE)

# Join together
beta.distdisp <- beta.distdisp.fr %>%
  bind_rows(beta.distdisp.to) %>%
  distinct() %>%
  filter(node_from %in% beta.list & node_to %in% beta.list) %>%
  distinct() %>%
  select(node_from, node_to, weight, everything())

# HGT, plasmids from
beta.hgt.fr <- amr.types1 %>%
  filter(beta.lact==1) %>%
  left_join(., hgt.netdat.et, by = c("node_id"="node_from"), keep=TRUE) 

# HGT, plasmids to
beta.hgt.to <- amr.types1 %>%
  filter(beta.lact==1) %>%
  left_join(., hgt.netdat.et, by = c("node_id"="node_to"), keep = TRUE)

# Join together
beta.hgt <- beta.hgt.fr %>%
  bind_rows(beta.hgt.to) %>%
  distinct() %>%
  filter(node_from %in% beta.list & node_to %in% beta.list) %>%
  distinct() %>%
  select(node_from, node_to, weight, everything())

# Repeat for tetracycline resistance
# Tetracycline resistant plasmids
tet.list <- c(1293, 3066, 3455)

# Recent dispersal, plasmids from
tet.recdisp.fr <- amr.types1 %>%
  filter(tet.res==1) %>%
  left_join(., recdisp.netdat.et, by = c("node_id"="node_from"), keep=TRUE) 

# Recent dispersal, plasmids to
tet.recdisp.to <- amr.types1 %>%
  filter(tet.res==1) %>%
  left_join(., recdisp.netdat.et, by = c("node_id"="node_to"), keep = TRUE)

# Join together
tet.recdisp <- tet.recdisp.fr %>%
  bind_rows(tet.recdisp.to) %>%
  distinct() %>%
  filter(node_from %in% tet.list & node_to %in% tet.list) %>%
  distinct() %>%
  select(node_from, node_to, weight, everything())

# Distant dispersal, plasmids from
tet.distdisp.fr <- amr.types1 %>%
  filter(tet.res==1) %>%
  left_join(., distdisp.netdat.et, by = c("node_id"="node_from"), keep=TRUE) 

# Distant dispersal, plasmids to
tet.distdisp.to <- amr.types1 %>%
  filter(tet.res==1) %>%
  left_join(., distdisp.netdat.et, by = c("node_id"="node_to"), keep = TRUE)

# Join together
tet.distdisp <- tet.distdisp.fr %>%
  bind_rows(tet.distdisp.to) %>%
  distinct() %>%
  filter(node_from %in% tet.list & node_to %in% tet.list) %>%
  distinct() %>%
  select(node_from, node_to, weight, everything())

# HGT, plasmids from
tet.hgt.fr <- amr.types1 %>%
  filter(tet.res==1) %>%
  left_join(., hgt.netdat.et, by = c("node_id"="node_from"), keep=TRUE)

# HGT plasmids to
tet.hgt.to <- amr.types1 %>%
  filter(tet.res==1) %>%
  left_join(., hgt.netdat.et, by = c("node_id"="node_to"), keep = TRUE)

# Join together
tet.hgt <- tet.hgt.fr %>%
  bind_rows(tet.hgt.to) %>%
  distinct() %>%
  filter(node_from %in% tet.list & node_to %in% tet.list) %>%
  distinct() %>%
  select(node_from, node_to, weight, everything())
####################################################################################################################