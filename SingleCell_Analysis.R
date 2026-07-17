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
plot1 + plot2 + plot3

# Filter the dataset based on QC thresholds 
# Store the initial number of cells to calculate how many we discard 
initial_cells <- ncol(PanIslets_Nicola)
cat("Initial number of cells:", initial_cells, "\n")

# Subset the object based on thresholds 
# (Note: keeping cells with features > 100, features < 5600, and MT % < 20)
PanIslets_Nicola <- subset(PanIslets_Nicola, subset = nFeature_RNA > 100 & nFeature_RNA < 5600 & percent.mt < 20)
PanIslets_Nicola

# Calculate and print the discarded cells for the presentation
final_cells <- ncol(PanIslets_Nicola)
discarded_cells <- initial_cells - final_cells

cat("Final number of cells:",final_cells, "\n")
cat("Number of dicarded cells:", discarded_cells,"\n")    

############## NORMALIZATION #####################
# Normalization 
# 10x data are usually transformed into counts per 10,000 reads. 
# The final expression estimate used for downstream analyses is the log of normalized count 
PanIslets_Nicola <- NormalizeData(PanIslets_Nicola, normalization.method = "LogNormalize", scale.factor = 10000)
PanIslets_Nicola@assays #25245 genes and 3060 cells

# Calculate average gene expression across all cells and show the top 50 
gene.expression <- apply(PanIslets_Nicola[["RNA"]]$data,1,mean)
gene.expression <- sort(gene.expression, decreasing = TRUE)

# Top 50 genes that have the highest mean expression across our cells
head(gene.expression, n = 50)

# Plot expression of a couple of highly expressed housekeeping genes
VlnPlot(PanIslets, features = c("RPS2","GAPDH"))

# Cell Cycle Scoring
cc.genes.updated.2019
PanIslets_Nicola <- CellCycleScoring(PanIslets_Nicola, 
                                     s.features = cc.genes.updated.2019$s.genes, 
                                     g2m.features = cc.genes.updated.2019$g2m.genes, 
                                     set.ident = TRUE) 

# Each cell is a point in a n-dimensional space, where n is the number of genes considered 
# The closer two points, the more similar are the transcriptomes of the corresponding cells

# Feature Selection (Highly Variable Genes)
# The 'vst' method estimates the mean-variance relationship and chooses the 2000 most variable genes
PanIslets_Nicola <- FindVariableFeatures(PanIslets_Nicola, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(PanIslets_Nicola), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(PanIslets_Nicola)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
MostVarGenes <- plot1 + plot2
MostVarGenes + ggtitle("Top 10 Variable Features in PanIslets")

# Scaling the Data
# A linear transformation that is a standard pre-processing step prior to dimensional reduction techniques like PCA
all.genes <- rownames(PanIslets_Nicola)
PanIslets_Nicola <- ScaleData(PanIslets_Nicola, features = all.genes)

# 