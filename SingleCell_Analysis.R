setwd("/Volumes/Archivio_MacNicola/Progetti UNIVERSITA/Transcriptomics_project_SC")
getwd

# Data Loading & Setup

library(dplyr)

#install.packages("Seurat")
library(Seurat)
library(patchwork)
library(ggplot2)

# Load PanglaoDB data
load("SRA701877_SRS3279692.sparse.RData")

# Mini-Challenge: Clean up gene names
# Remove the ENSEMBL ID (everything after the underscore)

