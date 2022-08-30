####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 17: Run Gillespie model for the peripheral plasmids, using HPC at BGU
# 
# 
# Script tested for R version 4.1.1,  run on the BGU HPC
####################################################################################################################

#! /gpfs0/shai/projects/R4/R-4.0.3/bin/Rscript
.libPaths("/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/library")
print(.libPaths())
print(sessionInfo())

# Load the necessary libraries
library(tidyverse)

# Load the data needed to run the script (upload to folder in which you are running this script on the HPC)
# Matrix, created in 13_Gillespie_model_setup.R
load("mat.Rda")

# Matrix id's, created in 13_Gillespie_model_setup.R
load("mat.ids.Rda")

# Plasmids to run model on, created in 14_Pick_starting_plasmids.R
load("not.mod1.sampled.Rda")
load("bool.periph.obs.Rda")

# Create the Gillespie function: 
gillespie <- function(boolean, # Your state nodes with or without gene
                      sim.matrix, # Matrix defining the flow
                      pars, # Parameters of the model
                      time_v){ #And a vector of times when things happen
  # Rates
  p <- pars$success_probability # For gene transfer
  cont <- pars$contact_rate # Rate at which plasmids contact (NOTE that I see this as a total rate, not a per capita one)
  loss <- pars$loss_rate # Gene loss per capita rate (NOTE that this changes when we gain or loose genes)
  
  # States: note that the state of the system is encoded in the boolean but it
  # will be handy to track the number of genes in the (meta)population.
  G <- sum(boolean)
  
  n_plasmids <- length(boolean) # Needed for later
  
  # Propensities
  r1 <- function(){cont}  # Two plasmids collide
  r2 <- function(){loss * G} # Losing a gene
  
  # Vector of propensities
  rates <- c(r1(), r2())
  
  # Names of events
  events <- c("f1", "f2")
  # Events: each one actualizes the corresponding states and propensities.
  f1 <- function(){ #Two plasmids collide
    x1 <- sample.int(n_plasmids, 1)
    x2 <- sample.int(n_plasmids, 1)
    
    if(boolean[x1] != boolean[x2]){
      gain <- runif(1) < (p * sim.matrix[x1, x2])
      if (gain) boolean[x1] <<- T; boolean[x2] <<- T; G <<- G + 1; rates[2] <<- r2() 
    }
  }
  f2 <- function(){ #A gene is loss at random
    index <- which(boolean)
    lost <- sample.int(G, 1)
    boolean[index[lost]] <<- F
    G <<- G - 1
    rates[2] <<- r2()
    rm(index); rm(lost)
  }
  
  
  # Initialize variables
  t <- time_v[1]
  out <- matrix(nrow = n_plasmids, ncol = length(time_v))
  out[, 1] <- boolean
  
  # Note: a little fraction of time is out of sync. Correct when everything else is settled.
  # Size of error depends on contact rate (higher contact rate = lower error)
  
  # Start simulations
  for(i in 2:length(time_v)){
    while (t < time_v[i]){
      
      # Sampling
      dt <- stats::rexp(1, sum(rates))
      
      z <- sample.int(length(rates), 1, prob = rates)
      do.call(events[z], args = list())
      
      t <- t + dt
    }
    out[, i] <- boolean #saving the states
  }
  
  out
}


# Create a for-loop to run the above Gillespie function for each starting plasmid and each parameter

# Set up the parameters


# Set up the parameters
# Success rate:
success <- 1

# Contact rates between plasmids:
cr.list <- c(10, 100, 1000)

# Loss rate of gene:
lr.list <- c(0, 0.01, 0.1)

# Number of repetitions
reps.ls <- c(seq(1, 300, 1))

# Create a dataframe of the combination of each parameters:
pars.df <- expand.grid(success, cr.list, lr.list)

# Set column names of the parameters dataframe
colnames(pars.df) <- c("success_probability", "contact_rate", "loss_rate") 

# Combine the parameters together in a single dataframe 
pars.reps.df <- expand.grid(success, cr.list, lr.list,reps.ls, names(bool.periph.obs)) %>%
  rownames_to_column(., var="df.id")

# Set column names of the parameters dataframe
colnames(pars.reps.df) <- c("df.id","success_probability", "contact_rate", "loss_rate", "sim.rep", "plasmid.rep") 


# Actual for-loop for running the model on each starting plasmid with each parameter combination
# Create an empty list where you will store the results:
ResultList.low<-list()

# Loop over each starting plasmid, parameter combination, and iteration of the simulation (300) 
# Results of each iteration of the loop will populate ResultList
for(i in 1:length(bool.periph.obs)){
  for(j in 1:length(reps.ls))
  {
    for(k in 1:nrow(pars.df))
    {
      boolean <- bool.periph.obs[[i]]
      # Run the Gillespie model with given parameters for 1000 time-steps, save result as a dataframe
      temp <- as.data.frame(gillespie(boolean = boolean, 
                                      sim.matrix = mat, 
                                      pars = list(success_probability = pars.df[k,1], contact_rate = pars.df[k,2], loss_rate = pars.df[k,3]), 
                                      time_v = 0:1000)) 
      
      resList<-(temp)
      
      # Append results to the model list
      ResultList.low[[length(ResultList.low)+1]]<-resList
    }
  }
}

# Save the Result list
save(ResultList.low, file= "ResultList.low.Rda")

# Some organizing / housekeeping of the results for each iteration of the model
step1 <- lapply(ResultList.low, function(x) as.data.frame(x) %>%
                  rownames_to_column(., var="mat.order") %>%
                  mutate(mat.order = as.double(mat.order)) %>%
                  left_join(., mat.ids) %>%
                  select(-mat.order) %>%
                  select(layer_id, everything()) %>%
                  group_by(layer_id) %>% 
                  summarise(across(starts_with("V"), sum)) %>%
                  column_to_rownames(., var="layer_id") %>%
                  mutate_at(vars(starts_with("V")), ~as.logical(.)))

# Sum the number of COWS with the gene in each time step
step2 <- lapply(step1, function(x) as.data.frame(colSums(x)) %>%
                  mutate(time.step = 1:length(x)) %>%
                  dplyr::rename(with.gene=`colSums(x)`))

# Name each list
names(step2) <- c(1:length(step2))

# Match each list to the parameters of the list, in pars.reps.df
# Each row with a set of parameters is assigned a df.id
pars.reps.df

# Create a new variable in each dataframe with a df.id, corresponding to the dataframe number in the list
step3 <- map2(step2, names(step2), ~ mutate(.x, df.id = .y))

# Save outputs
save(step2, file = "step2.Rda")
save(step3, file = "step3.Rda")
