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

###################### QUALITY CONTROL (QC) ################################
# Calculate the % of mitochondrial contamination, we use the MT prefix
grep("^MT-",rownames(PanIslets_Nicola),value = TRUE)
# The [[ operator can add columns to object metadata. This is a great place to store additional info/data
PanIslets_Nicola[["percent.mt"]] <- PercentageFeatureSet(PanIslets_Nicola, pattern = "^MT-")

# We also calculate ribosomal protein genes as a control 
grep("^RP[LS]",rownames(PanIslets_Nicola),value = TRUE) # ribosomal protein genes
PanIslets_Nicola[["percent.rbp"]] <- PercentageFeatureSet(PanIslets_Nicola, pattern = "^RP[LS]")

# Check the first 5 rows of the metadata to verify the new columns 
head(PanIslets_Nicola@meta.data, 5)

# Visualize QC metrics as violin plots
# First visualization (with dots)
VlnPlot(PanIslets_Nicola, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), ncol = 4)

# Second visualization (cleaner look without dots)
VlnPlot(PanIslets_Nicola, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), 
        ncol = 4, 
        pt.size = 0, 
        cols = c("#78c679")) 

# Visualize feauture-feauture relationships with ggplot2
meta_data <- PanIslets_Nicola@meta.data
#plot1 <- FeatureScatter(PanIslets, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(PanIslets, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Calculate Pearson correlations
cor1 <- cor(meta_data$nCount_RNA, meta_data$percent.mt, method = "pearson")
cor2 <- cor(meta_data$nCount_RNA, meta_data$nFeature_RNA, method = "pearson")
cor3 <- cor(meta_data$nCount_RNA, meta_data$percent.rbp, method = "pearson")

# Format the correlation test
cor_label1 <- paste0("Pearson's r = ", round(cor1, 2))
cor_label2 <- paste0("Pearson's r = ", round(cor2, 2))
cor_label3 <- paste0("Pearson's r = ", round(cor3, 2))

# Plot 1: nCount_RNA vs percent.mt
plot1 <- ggplot(meta_data, aes(x = nCount_RNA, y = percent.mt)) +
  geom_point(color = "#addd8e", alpha = 0.6, size = 1) +
  theme_minimal() +
  labs(title = "nCount_RNA vs percent.mt", subtitle = cor_label1,
       x = "nCount_RNA", y = "percent.mt") +
  theme(
    plot.title = element_text(face = "bold"),           # Bold title
  )

# Plot 2: nCount_RNA vs nFeature_RNA
plot2 <- ggplot(meta_data, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(color = "#78c679", alpha = 0.6, size = 1) +
  theme_minimal() +
  labs(title = "nCount_RNA vs nFeature_RNA", subtitle = cor_label2,
       x = "nCount_RNA", y = "nFeature_RNA") +
  theme(
    plot.title = element_text(face = "bold"),           # Bold title
  )

# Plot 3: nCount_RNA vs percent.rbp
plot3 <-  ggplot(meta_data, aes(x = nCount_RNA, y = percent.rbp)) +
  geom_point(color = "#31a354", alpha = 0.6, size = 1) +
  theme_minimal() +
  labs(title = "nCount_RNA vs percent.rbp", subtitle = cor_label3,
       x = "nCount_RNA", y = "percent.rbp") +
  theme(
    plot.title = element_text(face = "bold"),           # Bold title
  )

# Display the three scatter plots together using patchwork
plot1 + plo2 + plot3

# Filter the dataset based on QC thresholds 
# Store the initial number of cells to calculate how many we discard 
initial_cells <- ncol(PanIslets_Nicola)
cat("Initial number of cells:", initial_cells, "\n")

# Subset the object based on thresholds 
# (Note: keeping cells with features > 100, features < 5600, and MT % < 20)
PanIslets_Nicola <- subset(PanIslets, subset = nFeature_RNA > 100 & nFeature_RNA < 5600 & percent.mt < 20)
PanIslets_Nicola

# Calculate and print the discarded cells for the presentation
final_cells <- ncol(PanIslets_Nicola)
discarded_cells <- initial_cells - final_cells

cat("Final number of cells:",final_cells, "\n")
cat("Number of dicarded cells:", disdarded_cells,"\n")