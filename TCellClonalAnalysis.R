#  Install Core Packages
# install.packages("BiocManager")
# install.packages("Seurat")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("patchwork")
# install.packages("remotes")
# BiocManager::install("scRepertoire")
# BiocManager::install(version = "3.20")
# 
# Install Bioconductor Packages
# BiocManager::install("IRanges")
# 
# # remotes::install_github("carmonalab/scGate", ref="v1.6.2")
# remotes::install_github("carmonalab/STACAS")
# remotes::install_github("carmonalab/ProjecTILs")

# Load Packages to Verify Installation
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(IRanges)
library(remotes)
library(BiocManager)
library(ProjecTILs)
library(scRepertoire)
library(stringr)


#===============================Data Preprocessing==============================

# Define the location of the data
data_dirs <- list(
  "X_112_BL" = "data/X_112_BL/sample_filtered_feature_bc_matrix",
  "X_112_9M" = "data/X_112_9M/sample_filtered_feature_bc_matrix",
  "X_113_BL" = "data/X_113_BL/sample_filtered_feature_bc_matrix",
  "X_113_9M" = "data/X_113_9M/sample_filtered_feature_bc_matrix"
)


# Create a list to store the Seurat objects
seurat_objects <- list()

# Iterate over each directory to read the data
for (dir in names(data_dirs)) {
  data <- Read10X(data.dir = data_dirs[[dir]])
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = data, project = "HannanLab", min.cells = 10, min.features = 200)
  
  # Add metadata
  seurat_obj <- AddMetaData(seurat_obj, c(dir), col.name = 'sample')
  patient_id <- substr(dir, 1, nchar(dir) - 3)
  seurat_obj <- AddMetaData(seurat_obj, c(patient_id), col.name = 'patient')

  if (grepl("BL", dir)) {
    seurat_obj <- AddMetaData(seurat_obj, c('BL'), col.name = 'treatment')
  } else {
    seurat_obj <- AddMetaData(seurat_obj, c('9M'), col.name = 'treatment')
  }
  
  # Add to list of Seurat objects
  seurat_objects[[dir]] <- seurat_obj
}

# Read TCR data for all four datasets
tcr_data <- list()
for (dir in names(data_dirs)) {
  tcr_data[[dir]] <- read.csv(paste0("data/", dir, "/clonotypes.csv"), sep = ",")
  
  # Check whether there is a alpha-beta pair
  alpha <- str_count(tcr_data[[dir]]$cdr3s_aa, pattern = "TRA:")
  beta <- str_count(tcr_data[[dir]]$cdr3s_aa, pattern = "TRB:")
  
  tcr_data[[dir]]$has_alphabeta <- ifelse(alpha > 0 & beta > 0, 1, 0)
  tcr_data[[dir]]$has_singlepair <- ifelse(alpha == 1 & beta == 1, 1, 0)
}


# Merge TCR data with Seurat objects
for (dir in names(seurat_objects)) {
  airr_data <- read.csv(paste0("data/", dir, "/airr_rearrangement.tsv"), sep = "\t")
  airr_data <- merge(airr_data, tcr_data[[dir]], by.x = "clone_id", by.y = 'clonotype_id', all.x = TRUE, all.y = TRUE)
  
  # Remove duplicates based on barcode and keep only the first clonotype
  unique_airr_data <- airr_data[!duplicated(airr_data$cell_id), ]
  
  # Change the rownames to barcode
  row.names(unique_airr_data) <- unique_airr_data$cell_id
  
  # Add to Seurat object
  seurat_objects[[dir]] <- AddMetaData(seurat_objects[[dir]], unique_airr_data)
}


# Make unique IDs for TCRs per patient
for (dir in names(seurat_objects)) {
  seurat_objects[[dir]]$cdr3s_pat <- NA
  inds <- !is.na(seurat_objects[[dir]]$cdr3s_aa)
  seurat_objects[[dir]]$cdr3s_pat[inds] <- paste0(seurat_objects[[dir]]$cdr3s_aa[inds],
                                                  "_", seurat_objects[[dir]]$patient[inds])
}


# Merge all Seurat objects
# Initialize the first Seurat object as the starting point
merged_seurat <- seurat_objects[[1]]

# Loop over the remaining Seurat objects and merge them one by one
for (i in 2:length(seurat_objects)) {
  merged_seurat <- merge(merged_seurat, seurat_objects[[i]],
                         add.cell.ids = c(names(seurat_objects)[1], names(seurat_objects)[i]))
}


# Subset on pre-treatment samples
table(merged_seurat$treatment, merged_seurat$patient)



#================================ProjecTILs analysis============================

## Load reference CD8 T map
ref.file <- "CD8T_human_ref_v1.rds"
if (!file.exists(ref.file)) {
  download.file("https://figshare.com/ndownloader/files/38921366", destfile = ref.file)
}
ref.cd8 <- load.reference.map(ref.file)  # Use variable ref.file directly

# Visualize the reference CD8 T cell subtypes
DimPlot(ref.cd8, 
        cols = ref.cd8@misc$atlas.palette, 
        label = TRUE, 
        raster = TRUE) +
  theme(aspect.ratio = 1) + 
  ggtitle("CD8 T reference") + 
  NoLegend()

# Optional: Save the plot
ggsave("plots/CD8T_reference_umap.pdf", height = 4, width = 5)



# Key marker genes for the reference subtypes
fig.height <- 7
fig.width <- 4
DefaultAssay(ref.cd8) <- "RNA"

genes <- c("SELL", "TCF7", "LEF1", "CCR7", "S1PR1", "LMNA", "IL7R", "GZMK", "FGFBP2",
           "FCGR3A", "XCL1", "XCL2", "CD200", "CRTAM", "GNG4", "TOX", "PDCD1", "HAVCR2",
           "GZMB", "PRF1", "LAG3", "KLRB1", "TRAV1-2")


genes_2 = c("TCF7", "CCR7", "PDCD1", "TNFRSF9")


# Plot violin plots for key marker genes
VlnPlot(ref.cd8, 
        features = genes, 
        stack = TRUE, 
        flip = TRUE, 
        fill.by = "ident", 
        cols = ref.cd8@misc$atlas.palette) + 
  NoLegend()

# Optional: Save the violin plot
ggsave("plots/CD8T_ref_subset_violins.pdf", height = 12, width = 6)

# Switch to the integrated assay
DefaultAssay(ref.cd8) <- "integrated"


# Subset merged Seurat object based on treatment
merged_seurat_BL <- subset(merged_seurat, subset = treatment == "BL")
merged_seurat_9M <- subset(merged_seurat, subset = treatment == "9M")

# Split Seurat objects by sample
merged_seurat_list_BL <- SplitObject(merged_seurat_BL, split.by = "patient")
merged_seurat_list_9M <- SplitObject(merged_seurat_9M, split.by = "patient")

# Run ProjecTILs analysis for the BL timepoint (using 6 cores)
merged_seurat_projected_BL <- Run.ProjecTILs(merged_seurat_list_BL, ref = ref.cd8, ncores = 6)

# Run ProjecTILs analysis for the 9M timepoint (optional, can be added if needed)
merged_seurat_projected_9M <- Run.ProjecTILs(merged_seurat_list_9M, ref = ref.cd8, ncores = 6)

# Marker genes for expression profile verification
fig.height <- 8
fig.width <- 14

genes4radar <- c("CD4", "CD8A", "TCF7", "CCR7", "IL7R", "LMNA", "GZMA", "GZMK", "FGFBP2",
                 "XCL1", "CRTAM", "TOX", "PDCD1", "HAVCR2", "PRF1", "GZMB", "TRAV1-2", "KLRB1")

samples <- c("X_112_BL", "X_112_9M", "X_113_BL", "X_113_9M")
for (query_sample in samples) {
  # Extract the timepoint (e.g., "BL" or "9M") from the sample name
  time <- gsub(".*_(BL|9M)$", "\\1", query_sample)  
  seurat_obj <- get(paste0("merged_seurat_projected_", time))
  patient_id <- substr(query_sample, 1, nchar(dir) - 3)
  # Access the query data (seurat_obj$query)
  query_data <- seurat_obj[[patient_id]]
  
  plot_title <- paste("Radar Plot for", query_sample)
  plot_filename <- paste("plots/radars_projected_", query_sample, ".pdf", sep = "")
  
  plot.states.radar(ref.cd8, query = query_data, genes4radar = genes4radar) 
  
  # Save radar plot dynamically (adjust file name for each timepoint)
  ggsave(plot_filename, height = 11, width = 16)
}



##What is the subtype composition in different patients? Visualize projection and composition for individual patients:

# query_sample <- "X_112_BL"
# a <- plot.projection(ref.cd8, query = merged_seurat_projected_BL[[query_sample]],
#                      pointsize = 0.1, linesize = 0.3, raster=T) + NoLegend() +
#   ggtitle(query_sample)
# b <- plot.statepred.composition(ref.cd8, query = merged_seurat_projected_BL[[query_sample]]) 
# a | b
# plot_title <- paste("plots/single_patient_composition Plot for", query_sample, ".pdf")
# ggsave(plot_title, height=4, width=8)

samples <- c("X_112_BL", "X_112_9M", "X_113_BL", "X_113_9M")
for (query_sample in samples) {
  # Extract the timepoint (e.g., "BL" or "9M") from the sample name
  time <- gsub(".*_(BL|9M)$", "\\1", query_sample)
  
  seurat_obj <- get(paste0("merged_seurat_projected_", time))
  # Access the query data (seurat_obj$query)
  print(seurat_obj)
  patient_id <- substr(query_sample, 1, nchar(dir) - 3)
  print(patient_id)
  query_data <- seurat_obj[[patient_id]]
  
  a <- plot.projection(ref.cd8, query = query_data,
                       pointsize = 0.1, linesize = 0.3, raster=T) + NoLegend() +
    ggtitle(query_sample)
  
  b <- plot.statepred.composition(ref.cd8, query = query_data) 
  a | b
  plot_title <- paste("plots/single_patient_composition Plot for", query_sample, ".pdf")
  ggsave(plot_title, height=4, width=8)
}


###Remove samples with very few cells. A minimum number of cells is required for robust analyses.
# Define figure size
fig.height <- 7
fig.width <- 14

# Get the number of cells (columns) for each sample in BL and 9M datasets
sizes_BL <- as.vector(lapply(merged_seurat_projected_BL, ncol))
sizes_9M <- as.vector(lapply(merged_seurat_projected_9M, ncol))

# Keep only samples with more than 100 cells for both BL and 9M datasets
keep_BL <- names(sizes_BL)[sizes_BL > 100]
keep_9M <- names(sizes_9M)[sizes_9M > 100]

# Filter the Seurat objects
merged_seurat_projected_BL <- merged_seurat_projected_BL[keep_BL]
merged_seurat_projected_9M <- merged_seurat_projected_9M[keep_9M]

# Combine both datasets into one list
# merged_seurat_projected <- c(merged_seurat_projected_BL, merged_seurat_projected_9M)
merged_seurat_projected = merged_seurat_projected_BL
treatment = "_BL"


# Create the plots for each remaining sample in both BL and 9M datasets
plots <- lapply(names(merged_seurat_projected), function(sample_name) {
  plot.statepred.composition(ref.cd8, query = merged_seurat_projected[[sample_name]],
                             metric = "Percent") + ggtitle(paste0(sample_name, treatment))+ theme(legend.position = "none")
})
plots[[1]] <- plots[[1]] + theme(legend.position = "right")  # Add legend to the first plot
# Name the plots using the sample names
names(plots) <- names(merged_seurat_projected)

# Sort by the fraction of Tex cells for both BL and 9M datasets
tex.fraction <- lapply(merged_seurat_projected, function(sample) {
  # Calculate the fraction of Tex cells in each sample
  cluster_table <- table(sample$functional.cluster)
  tex_fraction <- unname(cluster_table["CD8.TEX"] / sum(cluster_table))
  tex_fraction
})

# Convert tex.fraction to a named vector for easy sorting
tex.fraction <- unlist(tex.fraction)

# Sort the plots by tex.fraction (descending)
plots <- plots[names(sort(tex.fraction, decreasing = TRUE))]

# Combine and display the plots in a single column
combined_plot <- wrap_plots(plots, ncol = length(names(merged_seurat_projected)))

# Display the combined plot
combined_plot

# Optionally, save the plot as a PDF
ggsave("plots/all_patient_composition_BL.pdf", combined_plot, height = fig.height, width = fig.width)


## Merge data into a single object to enable analyses across patients
# Define figure size
fig.width <- 12
fig.height <- 8

# Merge Seurat objects into a single object
merged_hannan_projected <- Reduce(merge.Seurat.embeddings, merged_seurat_projected)

# Set the identity class to "functional.cluster"
Idents(merged_hannan_projected) <- "functional.cluster"

# Set the default assay to "RNA"
DefaultAssay(merged_hannan_projected) <- "RNA"


# Clonotype analysis

## First, calculate frequencies over the total of cells with a TCR chain in the same sample:
# Define figure size
fig.width <- 5
fig.height <- 3

# Step 1: Calculate frequencies of clonotypes across all samples
freqs <- lapply(merged_seurat_projected, function(x) {
  # Calculate the frequency of each clonotype, normalized by total number of TCR cells
  table(x$cdr3s_pat) / sum(!is.na(x$cdr3s_pat))
})

# Step 2: Combine the frequency tables into a single vector
freqs <- Reduce(c, freqs)

# Step 3: Sort the frequencies in decreasing order
sorted <- sort(freqs, decreasing = TRUE)

# Step 4: Get the top 6 largest clones
largest_clones <- head(sorted, 6)

# Optional: View the result
print(largest_clones)


### Plot individual clones in the reference space:
fig.width=10
fig.height=6
plots <- list()
for (i in 1:length(largest_clones)) {
  
  ctype <- names(largest_clones)[i]
  patient <- gsub("^[^_]+_(.*)$", "\\1", ctype)
  
  
  
  cells_in_clone <- which(merged_hannan_projected[["cdr3s_pat"]]==ctype)
  plots[[i]] <- plot.projection(ref.cd8, merged_hannan_projected[,cells_in_clone],
                                pointsize = 1, linesize = 0.3, raster=T) + 
    ggtitle(sprintf("%s - Freq %.1f %%", patient, 100*largest_clones[i])) + NoLegend()
  print(patient)
}
wrap_plots(plots, ncol = 3)
ggsave("plots/largest_clones_umap_BL.pdf", height=7, width=12)

## Clonal expansion
##Summarize clonal sizes:

cutoffs <- c(0, 1, 5, 20, 50)
names(cutoffs) <- c("Single (=1)","Small (<=5)","Medium (<=20)","Large (<=50)","Hyper (>50)")
call <- "cdr3s_pat"
tab <- table(merged_hannan_projected[[call]])
merged_hannan_projected$cloneSize <- NA
merged_hannan_projected$cloneType <- NA

for (r in 1:nrow(merged_hannan_projected[[]])) {
  e <- merged_hannan_projected@meta.data[r,call]
  if (!is.na(e) & e %in% names(tab)) {
    merged_hannan_projected@meta.data[r,"cloneSize"] <- tab[e]
  }
}

for (i in seq_along(cutoffs)) {
  merged_hannan_projected$cloneType[merged_hannan_projected$cloneSize > cutoffs[i]] <- names(cutoffs[i])
}
merged_hannan_projected$cloneType <- factor(merged_hannan_projected$cloneType, levels= rev(names(cutoffs)))

merged_hannan_projected$Frequency <- merged_hannan_projected$cloneSize

table(merged_hannan_projected$cloneType, merged_hannan_projected$cloneSize)


clonalDiversity(merged_hannan_projected, cloneCall = "cdr3s_pat") +
  scale_color_manual(values = ref.cd8@misc$atlas.palette)

clonalProportion(merged_hannan_projected, cloneCall = "cdr3s_pat")

clonalHomeostasis(merged_hannan_projected, cloneCall = "cdr3s_pat")

clonalCompare(merged_hannan_projected, 
                  top.clones = 5,
                  cloneCall="cdr3s_pat", 
                  graph = "alluvial")

clonalDiversity(merged_hannan_projected, cloneCall = "cdr3s_pat") +
  scale_color_manual(values = ref.cd8@misc$atlas.palette)


clonalOverlay(merged_hannan_projected, 
              reduction = "umap", 
              cut.category = "cloneSize",
              cutpoint = 50, 
              bins = 10, 
              facet.by = "sample") + 
  guides(color = "none")


#ggraph needs to be loaded due to issues with ggplot
library(ggraph)

#No Identity filter
clonalNetwork(merged_hannan_projected, 
              reduction = "umap", 
              group.by = "functional.cluster",
              filter.clones = NULL,
              filter.identity = NULL,
              cloneCall = "cdr3s_pat")



shared.clones <- clonalNetwork(merged_hannan_projected, 
                               reduction = "umap", 
                               group.by = "functional.cluster",
                               cloneCall = "cdr3s_aa", 
                               exportClones = TRUE)
head(shared.clones)

merged_hannan_projected <- highlightClones(merged_hannan_projected, 
                                 cloneCall= "cdr3s_aa", 
                                 sequence = c("TRB:CASSSTGYTNEKLFF;TRA:CGTDPTATDKLIF", 
                                              "TRB:CSAPAAGGLETQYF;TRA:CAVPSGGSYIPTF;TRA:CVVEGSSNDYKLSF"))

Seurat::DimPlot(merged_hannan_projected, group.by = "highlight") + 
  ggplot2::theme(plot.title = element_blank())


###A useful visualization is the occupied single-cell repertoire. 
##It measures the number of cells, and their clonal expansion type, for each of the cell subtypes:

merged_hannan_projected$functional.cluster <- gsub("\\.", "_", merged_hannan_projected$functional.cluster)
# Ensure that 'functional.cluster' is treated as a factor
# Ensure functional.cluster is a factor
merged_hannan_projected$functional.cluster <- factor(merged_hannan_projected$functional.cluster)

# Check the levels (optional)
levels(merged_hannan_projected$functional.cluster)

# Create the plot with custom colors for each functional cluster
# clonalOccupy(merged_hannan_projected, x.axis = "functional.cluster")

library(circlize)
library(scales)


large <- subset(merged_hannan_projected, subset = clonesize >= 5)

circles <- getCirclize(large,
                       cloneCall = "cdr3s_pat",
                       group.by = "functional.cluster",
                       proportion = TRUE)

#Just assigning the normal colors to each cluster
grid.cols <- hue_pal()(length(unique(merged_hannan_projected$functional.cluster)))
names(grid.cols) <- unique(merged_hannan_projected$functional.cluster)

#Graphing the chord diagram
chordDiagram(circles, self.link = 1, grid.col = grid.cols)


## Clonotype bias

##Find clones with significant enrichment for specific cell subtypes:

# ggsave("plots/clonotype_bias.pdf", height=3, width=5)

large <- subset(merged_hannan_projected, subset = cloneSize >= 5)

biased <- clonalBias(large, cloneCall = "cdr3s_pat", split.by = "patient",
                        group.by = "functional.cluster", min.expand = 10, exportTable = TRUE)

most.biased <- biased[order(biased$Z.score, decreasing = TRUE),]
head(most.biased)


clonalBias(merged_hannan_projected, cloneCall = "cdr3s_pat", split.by = "patient",
              group.by = "functional.cluster", min.expand = 10, n.boots = 20) +
  scale_color_manual(values = ref.cd8@misc$atlas.palette)

ggsave("plots/clonotype_bias_BL.pdf", height=3, width=5)


### See where the most "biased" clonotypes are found on the reference.

plots <- list()
for (i in 1:6) {
  
  ctype <- most.biased[i, "Clone"]
  patient <- most.biased[i, "Sample"]
  size <- most.biased[i, "ncells"]
  zscore <- most.biased[i, "Z.score"]
  cells_in_clone <- which(merged_hannan_projected[["cdr3s_pat"]]==ctype)
  title <- sprintf("Clone %s, size %s, Z=%.1f", i, size, zscore)
  
  plots[[i]] <- plot.projection(ref.cd8, merged_hannan_projected[,cells_in_clone],
                                pointsize = 1, linesize = 0.3, raster=T) + 
    ggtitle(title) + NoLegend()
}

wrap_plots(plots, ncol = 3)

ggsave("plots/biased_clones_umap_BL.pdf", height=7, width=12)


