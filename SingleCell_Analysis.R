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
sm@Dimnames[[1]] <- sub("[_].*","", sm@Dimnames[[1]])
head(sm)

# Handle duplicate gene names
# Check if there are duplicates (if it returns a number > 0, there are duplicate)
anyDuplicated(rownames(sm))

# Make gene names unique, as required by the professor 
rownames(sm) <- make.unique(rownames(sm))

# Initialize the Seurat object 
# We name our object "PanIslets_Nicola"
PanIslets_Nicola <- CreateSeuratObject(counts = sm,
                                project = "PanIslets", min.cells = 3,
                                min.features = 200)

# Take a look at the newly created abject to see how many cells and genes we have 
PanIslets_Nicola
head(PanIslets_Nicola)