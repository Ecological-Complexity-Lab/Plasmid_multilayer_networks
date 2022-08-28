####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 19: Plot results of Gillespie dynamical models
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
library(ggplot2)
library(viridis)

# Load data processed on the server: 
# Central plasmids
load("sim.df.hi.Rda")
load("sim.df.hi.mean.Rda")

# Peripheral plasmids:
load("sim.df.low.Rda")
load("sim.df.low.mean.Rda")

####################################################################################################################
# Section 1: PLOT FOR CENTRAL PLASMIDS
####################################################################################################################

pd <- position_dodge(0.4)
plot1.high <- ggplot(sim.df.hi.mean, aes(x=time.step, y=mean.gene, color = contact_loss)) +
  geom_point(size=1.75, position = pd) +
  geom_line(size=1.25, position = pd) +
  scale_colour_manual(values = c("#f3ee27","#fad824","#feba2c",
                                 "#f2645c","#de4968","#c03a76",
                                 "#a1307e","#6e1e81","#2f1163"),
                      name = "",
                      labels = c("Low contact, no loss", "Low contact, intermediate loss", "Low contact, high loss",
                                 "Intermediate contact, no loss","Intermediate contact, intermediate loss", "Intermediate contact, high loss",
                                 "High contact, no loss","High contact, intermediate loss","High contact, high loss")) + 
  theme_classic() +
  xlab("Time step") +
  ylab("Cows with gene (mean)") +
  theme(legend.position="bottom",
        legend.text=element_text(size=11),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  guides(colour = guide_legend(nrow = 3))

# Export
ggsave(plot1.high, filename="plot1.high.png", dpi = 1800, width = 9, height = 6, units = "in")
ggsave(plot1.high, filename="plot1.high.pdf", width = 9, height = 6, units = "in")
####################################################################################################################


####################################################################################################################
# Section 2: PLOT FOR PERIPHERAL PLASMIDS
####################################################################################################################

pd <- position_dodge(0.4)
plot1.low <- ggplot(sim.df.low.mean, aes(x=time.step, y=mean.gene, color = contact_loss)) +
  geom_point(size=1.75, position = pd) +
  geom_line(size=1.25, position = pd) +
  scale_colour_manual(values = c("#f3ee27","#fad824","#feba2c",
                                 "#f2645c","#de4968","#c03a76",
                                 "#a1307e","#6e1e81","#2f1163"),
                      name = "",
                      labels = c("Low contact, no loss", "Low contact, intermediate loss", "Low contact, high loss",
                                 "Intermediate contact, no loss","Intermediate contact, intermediate loss", "Intermediate contact, high loss",
                                 "High contact, no loss","High contact, intermediate loss","High contact, high loss")) + 
  theme_classic() +
  xlab("Time step") +
  ylab("Cows with gene (mean)") +
  theme(legend.position="bottom",
        legend.text=element_text(size=11),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  guides(colour = guide_legend(nrow = 3))

# Save plots
ggsave(plot1.low, filename="plot1.low.png", dpi = 1800, width = 9, height = 6, units = "in")
ggsave(plot1.low, filename="plot1.low.pdf", width = 9, height = 6, units = "in")
####################################################################################################################


####################################################################################################################
# Section 3: RESULTS OF GILLESPIE DYNAMICAL MODELS
####################################################################################################################

# Time to getting to all 21, central plasmids
last.step.high <- sim.df.hi %>%
  filter(with.gene == 21) %>%
  ungroup() %>%
  group_by(sim.rep, plasmid.rep,contact_loss, .groups=T) %>%
  slice_min(time.step)
#slice_head()
#summarise(min.time.step=min(time.step))  

last.step.high.test <- sim.df.hi %>%
  filter(with.gene == 21 & contact_loss=="100_0") %>%
  ungroup() %>%
  group_by(sim.rep, plasmid.rep, .groups=T) %>%
  slice_min(time.step)

last.step.high2 <- last.step.high %>%
  group_by(contact_loss) %>%
  mutate(mean.time.step = mean(time.step),
         per.reps.21 = n()/3000*100) %>%
  ungroup() %>%
  select(contact_loss, mean.time.step, per.reps.21) %>%
  distinct()

last.step.high <- sim.df.hi %>%
  filter(with.gene == 21) %>%
  ungroup() %>%
  group_by(sim.rep, plasmid.rep,contact_loss, .groups=T) %>%
  slice_min(time.step)
#slice_head()
#summarise(min.time.step=min(time.step))  


# Repeat for peripheral plasmids
last.step.low <- sim.df.low %>%
  filter(with.gene == 21) %>%
  ungroup() %>%
  group_by(sim.rep, plasmid.rep,contact_loss, .groups=T) %>%
  slice_min(time.step)

last.step.low.not21 <- sim.df.low %>%
  filter(with.gene != 21) %>%
  ungroup() %>%
  group_by(sim.rep, plasmid.rep,contact_loss, .groups=T) %>%
  slice_max(time.step)

last.step.low2 <- last.step.low %>%
  group_by(contact_loss) %>%
  mutate(mean.time.step = mean(time.step),
         per.reps.21 = n()/3000*100) %>%
  ungroup() %>%
  select(contact_loss, mean.time.step, per.reps.21) %>%
  distinct()

# Mean for last step
last.step.low3 <- sim.df.low %>%
  filter(time.step == 300) %>%
  group_by(contact_loss) %>%
  mutate(mean.inf.300 = mean(with.gene))%>%
  select(contact_loss, mean.inf.300) %>%
  distinct()

last.step.low.under21 <- sim.df.low %>%
  filter(with.gene != 21) %>%
  ungroup() %>%
  group_by(sim.rep, plasmid.rep, .groups=T) %>%
  slice_max(time.step)

last.step.infected <- sim.df.low %>%
  ungroup() %>%
  group_by(sim.rep, plasmid.rep, .group=T) %>%
  slice_max(time.step)
