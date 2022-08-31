# Shapiro_etal_Plasmid_Networks
Data and scripts to carry out all analysis in Shapiro et al., "Multilayer Networks of Plasmid Genetic Similarity Reveal Potential Pathways of Gene Transmission"

All scripts tested for R version 4.1.1

This repository contains all files needed to run analysis on plasmid rumen networks. To begin, download and unzip the `StartingFiles.zip` folder and run analyses sequentially starting with script 01_Initial_data_processing.

Scripts can be started independently using the Rdata saved in the `Rdata_Intermediate_files.zip` folder. Rdata files are saved to a folder corresponding to the name of the script in which they were created, e.g. all Rdata outputs created in 03_Network_setup.R that are used downstream can be found in the folder Rdata_Intermediate_files/Outputs_Script03_Network_setup.zip.

Each script is sub-divided to sections. All the scripts are fully commented.

All scripts can be run on a local computer, besides the transmission model (scripts 15-18), which requires use of HPC and was run on the Ben-Gurion University server.
