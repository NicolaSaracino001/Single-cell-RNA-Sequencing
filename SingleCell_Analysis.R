setwd("/Volumes/Archivio_MacNicola/Progetti UNIVERSITA/Transcriptomics_project_SC")
getwd
library(dplyr)
#install.packages("Seurat")
library(Seurat)
library(patchwork)
library(ggplot2)
load("SRA701877_SRS3279692.sparse.RData")