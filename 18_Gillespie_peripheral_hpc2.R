####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 18: Run statistics on output of Gillespie dynamical model on peripheral plasmids (using BGU HPC)
# 
# 
# Script tested for R version 4.1.1
####################################################################################################################

#! /gpfs0/shai/projects/R4/R-4.0.3/bin/Rscript
.libPaths("/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/library")
print(.libPaths())
print(sessionInfo())

library(tidyverse)

# Outputs from script 17
load("step3.Rda")
load("pars.reps.df.Rda")


# Join the list of dataframes into one dataframe and then join by the df.id column in the 
# parameters dataframe in order to join the parameters to the results of each simulation
sim.df.low <- do.call(rbind.data.frame, step3) %>%
  left_join(., pars.reps.df, by="df.id") %>%
  select(-df.id) %>%
  mutate(contact_loss = paste(contact_rate,loss_rate, sep = "_"))

# Mean cows infected per time-step across simulations
sim.df.low.mean <- sim.df.low %>%
  group_by(time.step, contact_loss, loss_rate, contact_loss, groups=T) %>%
  summarise(mean.gene=mean(with.gene))

# Save outputs
save(sim.df.low, file="sim.df.low.Rda")

save(sim.df.low.mean, file="sim.df.low.mean.Rda")
