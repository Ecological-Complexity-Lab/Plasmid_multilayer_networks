####################################################################################################################
# Plasmid rumen network analysis
# 
# Script 7: Analyze the skew and potential power law distribution of degree
# 
#
# The following figures are created:
# Supplementary Figure S6
# 
# Script tested for R version 4.1.1
####################################################################################################################


####################################################################################################################
# SCRIPT SET-UP
####################################################################################################################
# Set working directory to wherever your files are located
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(moments)
library(igraph)
library(poweRlaw)


# Starting files:
# Degree of physical nodes, created in script 04_Basic_network_statistics.R
load("deg.str.2k.all.phys.Rda")

# Plasmid-to-cow links, created in script 04_Basic_network_statistics.R
load("layer_links.Rda")
####################################################################################################################


####################################################################################################################
# Section 1: SKEWNESS AND POWER LAW FOR DEGREE & LAYER LINKS
####################################################################################################################
# Calculate skewness:
skewness(deg.str.2k.all.phys$degree)
skewness(layer_links$links.to.cows)

# Using fit_power_law in igraph
# For degree (physical nodes):
pow.deg <- fit_power_law(deg.str.2k.all.phys$degree, implementation = "plfit")

# alpha 2.747697, xmin = 8, logLik=-394.4605, KS = 0.0549, KS.p = 0.8281631
# Note: Only 130 plasmids with degree >= 8

# Using poweRlaw package, compare the fit of different distributions:
# Power law
x.deg <- deg.str.2k.all.phys$degree
m1.deg = displ$new(x.deg) 
m1.deg$setXmin(estimate_xmin(m1.deg))

# Compare Power law to discrete exponential
m2.deg = disexp$new(x.deg) 
m2.deg$setXmin(m1.deg$getXmin()) 
est2.deg = estimate_pars(m2.deg) 
m2.deg$setPars(est2.deg$pars)

comp1.deg = compare_distributions(m1.deg, m2.deg)


# Compare Power law to discrete log-normal
m3.deg = dislnorm$new(x.deg) 
m3.deg$setXmin(m1.deg$getXmin()) 
est3.deg = estimate_pars(m3.deg) 
m3.deg$setPars(est3.deg$pars)

comp2.deg = compare_distributions(m1.deg, m3.deg)


# Compare Power law to discrete Poisson:
m4.deg = dispois$new(x.deg) 
m4.deg$setXmin(m1.deg$getXmin()) 
est4.deg = estimate_pars(m4.deg) 
m4.deg$setPars(est4.deg$pars)

comp3.deg = compare_distributions(m1.deg, m4.deg)
####################################################################################################################


####################################################################################################################
# Section 2: DEGREE HISTOGRAM AND LOG-LOG PLOT
####################################################################################################################
# Create Supplementary Figure S6: Degree distribution and log-log degree-frequency plot

# Count frequency of each value of degree:
deg.count <- plyr::count(deg.str.2k.all.phys$degree) %>%
  mutate(p.k = freq/1344,
         log.pk = log(p.k))

# Supplementary Figure S6, Panel A: Histogram of degree:
deg.plot <- ggplot(deg.str.2k.all.phys, aes(x=degree))+
  geom_histogram(color="black", fill="white") +
  theme_classic() +
  labs(x = "Degree (Physical nodes)", y="Plasmids (Count)") +
  theme(axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.y=element_text(size=22),
        axis.title.x=element_text(size=22)) +
  theme(plot.margin = margin(1.5, 0.75, 0.5, 0.5, "cm")) 


# Supplementary Figure S6, Panel B:Log-log plot of degree vs. probability of each degree:
log.deg.plot <- ggplot(deg.count, aes(x=log(x), y=log(p.k))) +
  geom_point(size = 4) +
  #geom_abline(intercept = -2, slope = -1) +
  xlab("k") +
  ylab("P(k)") +
  theme_classic() + 
  theme(axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.y=element_text(size=22),
        axis.title.x=element_text(size=22)) +
  theme(plot.margin = margin(1.5, 0.75, 0.5, 0.5, "cm")) +
  scale_y_continuous(breaks = c(-9.21, -6.9, -4.6,-2.3,0), 
                   labels = c("0.0001", "0.001", "0.01", "0.1","1"))

# Plot Panels A and B together with ggarrange:
ggarrange(deg.plot, log.deg.plot, 
                   labels = c("A.", "B."),
                  ncol = 1, nrow = 2,font.label = list(size = 26))

# Save plots in desired format:
ggsave(filename = "deg.hist.log.log.png",  width = 8, height = 10, units = "in")

ggsave(filename = "deg.hist.log.log.pdf",  width = 8, height = 10, units = "in")
####################################################################################################################