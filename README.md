# General information about code and data

The repository contains data and scripts to carry out all analysis in Shapiro et al., "Multilayer Networks of Plasmid Genetic Similarity Reveal Potential Pathways of Gene Transmission" (in revision, currently in BioRxiv [https://www.biorxiv.org/content/10.1101/2022.09.08.507140v1](https://www.biorxiv.org/content/10.1101/2022.09.08.507140v1)

All scripts tested in R version 4.1.1

This repository contains all files needed to run analysis on plasmid rumen networks. To begin, download and unzip the `StartingFiles.zip` folder and run analyses sequentially starting with script 01_Initial_data_processing.

Scripts can be started independently using the Rdata saved in the `Rdata_Intermediate_files.zip` folder. Rdata files are saved to a folder corresponding to the name of the script in which they were created, e.g. all Rdata outputs created in 03_Network_setup.R that are used downstream can be found in the folder Rdata_Intermediate_files/Outputs_Script03_Network_setup.zip.

Each script is sub-divided to sections. All the scripts are fully commented.

All scripts can be run on a local computer, besides the transmission model (scripts 15-18), which requires use of HPC and was run on the Ben-Gurion University server.

# Highlighted files

For convenience, we include a few result files that give direct information without needing to go through the whole analysis.

* `'plas.2k.name.node.id.csv`: A table with plasmid nam,es and ids
* `plasmid.2k.metadat.csv`: A table describing which plasmid occurs in which cow.
* `net.dat2k.ew.csv`: A table of edge lists (the network). In an extended edge list format
* `plas_mods.df.csv`: A table with the module assignment of each plasmid.
* `ann.metadat.v1.csv`: A table with module assignments and ORF annotations.
* `KOs_in_plasmids.txt`: A table with KEGG orthologies.

