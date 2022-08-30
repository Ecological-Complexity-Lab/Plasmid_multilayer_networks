####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 16: Run statistics on output of Gillespie dynamical model on central plasmids (using BGU HPC)
# 
# 
# Script tested for R version 4.1.1, run on the BGU HPC
####################################################################################################################

#! /gpfs0/shai/projects/R4/R-4.0.3/bin/Rscript
.libPaths("/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/library")
print(.libPaths())
print(sessionInfo())

# Load the necessary libraries
library(tidyverse)

# Load the data needed to run the script (upload to folder in which you are running this script on the HPC)
# Outputs from script 15_Gillespie_central_hpc.R
load("step3.hi.Rda")
load("pars.reps.df.Rda")


# Join the list of dataframes into one dataframe and then join by the df.id column in the 
# parameters dataframe in order to join the parameters to the results of each simulation
sim.df.hi <- do.call(rbind.data.frame, step3.hi) %>%
  left_join(., pars.reps.df, by="df.id") %>%
  select(-df.id) %>%
  mutate(contact_loss = paste(contact_rate,loss_rate, sep = "_"))

# Calculate the mean number of cows infected per time-step across simulations
sim.df.hi.mean <- sim.df.hi %>%
  group_by(time.step, contact_loss, loss_rate, contact_loss, groups=T) %>%
  summarise(mean.gene=mean(with.gene))

# Save the outputs
save(sim.df.hi, file="sim.df.hi.Rda")

save(sim.df.hi.mean, file="sim.df.hi.mean.Rda")