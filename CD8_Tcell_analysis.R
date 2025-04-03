# Load Packages to Verify Installation

install.packages("BiocManager")
install.packages("remotes")
remotes::install_github("carmonalab/STACAS")
remotes::install_github("carmonalab/ProjecTILs")

BiocManager::install("scRepertoire", force = TRUE)
BiocManager::install("scater", force = TRUE)
BiocManager::install("IRanges")
install.packages("circlize")


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
library(scater)
library(circlize)
library(scales)
library(ggraph)

#===============================Data Preprocessing==============================
datadir="data/"
cached_seurat <- sprintf("%s/HannanLab_data_seurat_3_10.rds", datadir)
if(file.exists(cached_seurat)){
  merged_seurat <- readRDS(cached_seurat)
}else {
  
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
    seurat_obj <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 10)
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
  sample_names <- names(seurat_objects)
  merged_seurat  <- merge(x = seurat_objects[[1]],
                          y = seurat_objects[2:length(seurat_objects)], 
                          add.cell.ids = sample_names)
  #Save object to disk
  saveRDS(merged_seurat, cached_seurat)
}

# Normalize the data using LogNormalize
# merged_seurat <- NormalizeData(merged_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
###Read TCR Informations
S1 = read.csv("data/X_112_BL/filtered_contig_annotations.csv", sep = ",")
S2 = read.csv("data/X_112_9M/filtered_contig_annotations.csv", sep = ",")
S3 = read.csv("data/X_113_BL/filtered_contig_annotations.csv", sep = ",")
S4 = read.csv("data/X_113_9M/filtered_contig_annotations.csv", sep = ",")
contig_list <- list(S1, S2, S3, S4)
head(contig_list[[1]])

### Combining Contigs into Clones
combined_TCR <- combineTCR(contig_list, samples = c("X_112_BL", "X_112_9M", 
                                                    "X_113_BL", "X_113_9M"),
                           removeNA = T,
                           removeMulti = F
                           )
### addVariable
# combined_TCR <- addVariable(combined_TCR, 
#                             variable.name = "time_frame", 
#                             variables = rep(c("BL", "9M"), 2))
# combined_TCR = lapply(combined_TCR, function(x){
#   x$cdr3s_pat = paste(x$CTaa, x$sample, sep = "_"); x})
# 
# head(combined_TCR[[1]])

# Extract cell barcodes from Seurat object
seurat_barcodes <- colnames(merged_seurat)

# Extract cell barcodes from TCR data (assuming `combined_TCR` contains cell barcodes)
# If your TCR data has different naming, replace 'barcode' with the correct column name
tcr_barcodes <- unique(unlist(lapply(combined_TCR, function(x) x$barcode)))
# Subset the Seurat object to retain only cells present in both datasets
merged_seurat_filtered <- merged_seurat[, seurat_barcodes %in% tcr_barcodes]

#Combine TCR data and Seurat Expression data
seurat_sc_combined <- combineExpression(combined_TCR, 
                                        merged_seurat_filtered, 
                                        # cloneCall="gene", 
                                        proportion = FALSE, 
                                        cloneSize=c(Single=1, Small=5, Medium=20,
                                                    Large=100, Hyperexpanded=500))

seurat_sc_combined@meta.data$cdr3s_pat = paste(seurat_sc_combined@meta.data$CTaa,
                                               seurat_sc_combined@meta.data$sample, 
                                               sep = "_")
seurat_sc_combined@meta.data$cdr3s_pat[is.na(seurat_sc_combined@meta.data$CTaa)] <- NA

#remove files after combining
rm(merged_seurat)
rm(merged_seurat_filtered)
rm(S1)
rm(S2)
rm(S3)
rm(S4)
rm(contig_list)

#===========================================ProjectIILS======================================================================
#----------------------------Load reference CD8 T map---------------------------
ref.file <- "CD8T_human_ref_v1.rds"
# ref.file = "CD4T_human_ref_v1.rds"

if (!file.exists(ref.file)) {
  download.file("https://doi.org/10.6084/m9.figshare.23608308", 
                destfile = ref.file)
}
ref.cd8 <- load.reference.map(ref.file)  # Use variable ref.file directly


fig.height = 5
fig.width = 6
# Visualize the reference CD8 T cell subtypes
DimPlot(ref.cd8, 
        cols = ref.cd8@misc$atlas.palette, 
        label = TRUE, 
        raster = TRUE) +
        theme(aspect.ratio = 1) + 
        ggtitle("CD8 T reference") + 
        NoLegend()

fig_name = "plots/CD8T_reference_umap.png"
# Optional: Save the plot
ggsave(fig_name, height = fig.height, width = fig.width, dpi = 300)
#-------------------------------------------------------------------------------
#-----------------------Key marker genes for the reference subtypes-------------
fig.height = 7
fig.width = 6
DefaultAssay(ref.cd8) = "RNA"

genes <- c("SELL", "PTPRC", "TCF7", "LEF1", "CCR7", "TNFRSF9", "S1PR1", "LMNA", "IL7R", "GZMK", "FGFBP2",
           "FCGR3A", "XCL1", "XCL2", "CD200", "CRTAM", "GNG4", "TOX", "PDCD1", "HAVCR2",
           "GZMB", "PRF1", "LAG3", "KLRB1", "TRAV1-2")

# Plot violin plots for key marker genes
VlnPlot(ref.cd8, 
        features = genes, 
        stack = TRUE, 
        flip = TRUE, 
        fill.by = "ident", 
        cols = ref.cd8@misc$atlas.palette) + 
  NoLegend()
fig_name = "plots/CD8T_ref_subset_violins.png"
ggsave(fig_name, height = fig.height, width = fig.width, , dpi = 300)
#-------------------------------------------------------------------------------
#-----------------------Run ProjecTILS and split by sample----------------------

seurat_projected_all = Run.ProjecTILs(seurat_sc_combined, ref = ref.cd8, 
                                      split.by = "sample", ncores = 6)

#-------------------------------------------------------------------------------

#------------------Change the factor order for other plots----------------------
lv_order = c( "Single (0 < X <= 1)", "Small (1 < X <= 5)", "Medium (5 < X <= 20)", 
              "Large (20 < X <= 100)", "Hyperexpanded (100 < X <= 500)")
seurat_projected_all$cloneSize = factor(seurat_projected_all$cloneSize, 
                                        levels = lv_order)
cluster_order = c("CD8.NaiveLike", "CD8.CM", "CD8.EM", "CD8.TEMRA", "CD8.TPEX",
                  "CD8.TEX", "CD8.MAIT")
seurat_projected_all$functional.cluster = factor(seurat_projected_all$functional.cluster, 
                                                 levels = cluster_order)


#-------------------------------------------------------------------------------
#----Compare Marker gene expression profiles of query data with reference map---
fig.height <- 10
fig.width <- 20

# genes4radar <- c("CD4", "CD8A", "PTPRC", "TCF7", "CCR7", "IL7R", "LMNA", "GZMA", "GZMK", "FGFBP2",
#                  "XCL1", "CRTAM", "TOX", "PDCD1", "HAVCR2", "PRF1", "GZMB", "TRAV1-2",
#                  "KLRB1")

genes4radar = c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "TCF7", "CCR7", "PDCD1", 
                "TNFRSF9", "GZMK", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "TOX")

#-------------------------------------------------------------------------------
#------Separate seurat object based on patient and plot sample wise radar ------

patient_id = "X_113"
patient_id = "all"

if(patient_id == "all"){
  query_list = SplitObject(seurat_projected_all, split.by = "sample")
  } else{
  query_patient <- subset(seurat_projected_all, subset = patient == patient_id)
  query_list = SplitObject(query_patient, split.by = "sample")
  }
  
plot.states.radar(ref.cd8, query_list, genes4radar = genes4radar, min.cells = 5)
plot_filename <- paste0("plots/radars_projected_", patient_id, ".png", sep = "")
ggsave(plot_filename, height = fig.height, width = fig.width, dpi = 300)

#-------------------------------------------------------------------------------
#---------------------Cell Subtype Composition of each sample-------------------
fig.height <- 6
fig.width <- 10

query_sample <- "X_113_BL"
query_list = SplitObject(seurat_projected_all, split.by = "sample")

a = plot.projection(ref.cd8, query = query_list[[query_sample]], ref.alpha = 0.5, ref.size=0.6,
                     pointsize = 1, alpha = 1, linesize = 0.2, raster=T) + NoLegend() +
                     ggtitle(query_sample)
b = plot.statepred.composition(ref.cd8, query = query_list[[query_sample]]) 

a | b
fig_name =  paste0("plots/single_patient_composition_", query_sample, ".png")
ggsave(fig_name, height = fig.height, width = fig.width, dpi = 300)

#-------------------------------------------------------------------------------

#---------------Plot cell type composition all samples together ----------------
fig.height <- 7
fig.width <- 14

query_list = SplitObject(seurat_projected_all, split.by = "sample")
calc_metric = "percent" #percent count
#if percentage
if (calc_metric == "percent"){
  y_limit_range <- c(0, 55)
} else {
  #if count
  y_limit_range <- c(0, 2600) #or change accordingly
}
plots = lapply(names(query_list), function(x) {
                    plot.statepred.composition(ref.cd8, query = query_list[[x]],
                                               labels.col="functional.cluster",
                                               metric = calc_metric ) + 
                                               ggtitle(x) + ylim(y_limit_range) +
                                               theme(legend.position = "none") 
                                              })
plots[[1]] = plots[[1]] + theme(legend.position = "left")  # Add legend to the first plot
names(plots) <- names(query_list)
combined_plot <- wrap_plots(plots, ncol = length(names(query_list)), sharey = TRUE)
# Display the combined plot
combined_plot
fig_name =  paste0("plots/all_patient_composition_", calc_metric, ".png")
ggsave(fig_name, combined_plot, height = fig.height, width = fig.width, dpi = 300)

#-------------------------------------------------------------------------------

#------------------Merge lists of object to obtain a single object--------------

# split objects into each sample first
query_list = SplitObject(seurat_projected_all, split.by = "sample")
# merge into a single object 
merged_projected = Reduce(merge.Seurat.embeddings, query_list)

lv_order = c("Hyperexpanded (100 < X <= 500)", "Large (20 < X <= 100)", 
             "Medium (5 < X <= 20)", "Small (1 < X <= 5)", "Single (0 < X <= 1)")
merged_projected$cloneSize = factor(merged_projected$cloneSize, 
                                        levels = lv_order)

merged_projected$functional.cluster = factor(merged_projected$functional.cluster, 
                                                 levels = cluster_order)
# set identification column as functional cluster
Idents(merged_projected) = "functional.cluster"

#-------------------------------------------------------------------------------

#============================================Clonal Analysis=====================================================================

#---------------largest clones per patient--------------------------------------
patient_id = "X_113"
projected_patient = subset(seurat_projected_all, subset = patient == patient_id)
query_list = SplitObject(projected_patient, split.by = "sample")

freqs_patient = lapply(query_list, function(x){
  table(x$cdr3s_pat)/sum(!is.na(x$cdr3s_pat))
})

freqs_patient = Reduce(c, freqs_patient)
sorted_patient = sort(freqs_patient, decreasing = TRUE)
largest_clones = head(sorted_patient, 6)
largest_clones

### Plot n largest clones into the reference map
fig.height = 7
fig.width = 14
plots = list()

for(i in 1:length(largest_clones)){
  ctype = names(largest_clones)[i]
  patient <- substr(ctype, nchar(ctype) - 7, nchar(ctype))
  cells = which(merged_projected[["cdr3s_pat"]] == ctype)
  plots[[i]] = plot.projection(ref.cd8, merged_projected[, cells],
                               pointsize = 1, linesize = 0.3, raster=T) + 
                               ggtitle(sprintf("%s - Freq %.1f %%", 
                               patient, 100*largest_clones[i])) + 
                               theme(legend.position = "none")
                               }
plots[[1]] <- plots[[1]] + theme(legend.position = "left")  # Add legend to the first plot
combined_plot = wrap_plots(plots, ncol = 3)
combined_plot
fig_name = paste0("plots/largest_clones_in_patient_", patient_id, ".png")
ggsave(fig_name, combined_plot, height = fig.height, width = fig.width, dpi = 300)

#-------------------------------------------------------------------------------

#-------------------largest clone per sample -----------------------------------

sample_id = "X_113_9M"
projected_sample = subset(seurat_projected_all, subset = sample == sample_id)
query_list = SplitObject(projected_sample, split.by = "sample")

freqs_sample = lapply(query_list, function(x){
  table(x$cdr3s_pat)/sum(!is.na(x$cdr3s_pat))
})

freqs_sample = Reduce(c, freqs_sample)
sorted_sample = sort(freqs_sample, decreasing = TRUE)
largest_clones = head(sorted_sample, 10)
largest_clones

# Plot n largest clones into the reference map
fig.height = 7
fig.width = 14
plots = list()

for(i in 1:length(largest_clones)){
  ctype = names(largest_clones)[i]
  patient <- substr(ctype, nchar(ctype) - 7, nchar(ctype))
  cells = which(merged_projected[["cdr3s_pat"]] == ctype)
  plots[[i]] = plot.projection(ref.cd8, merged_projected[, cells],
                               pointsize = 1, linesize = 0.3, raster=T) + 
    ggtitle(sprintf("%s - Freq %.1f %%", 
                    patient, 100*largest_clones[i])) + 
    theme(legend.position = "none")
}
plots[[1]] <- plots[[1]] + theme(legend.position = "left")  # Add legend to the first plot
combined_plot = wrap_plots(plots, ncol = 3)
combined_plot
fig_name = paste0("plots/largest_clones_in_sample_", sample_id, ".png")
ggsave(fig_name, combined_plot, height = fig.height, width = fig.width, dpi = 300)

#-------------------------------------------------------------------------------

#-------------------Most expanded shared clones between samples-----------------

sample_id_1 = "X_113_BL"
sample_id_2 = "X_113_9M"

# Subset the two samples
projected_sample_1 = subset(seurat_projected_all, subset = sample == sample_id_1)
projected_sample_2 = subset(seurat_projected_all, subset = sample == sample_id_2)

# Split the samples by "sample"
query_list_1 = SplitObject(projected_sample_1, split.by = "sample")
query_list_2 = SplitObject(projected_sample_2, split.by = "sample")

# Calculate the frequencies of clones for both samples
freqs_sample_1 = lapply(query_list_1, function(x){
  table(x$CTaa)/sum(!is.na(x$CTaa))
})
freqs_sample_2 = lapply(query_list_2, function(x){
  table(x$CTaa)/sum(!is.na(x$CTaa))
})

# Combine frequencies from the list
freqs_sample_1 = Reduce(c, freqs_sample_1)
freqs_sample_2 = Reduce(c, freqs_sample_2)

# Sort the frequencies for both samples in decreasing order
sorted_sample_1 = sort(freqs_sample_1, decreasing = TRUE)
sorted_sample_2 = sort(freqs_sample_2, decreasing = TRUE)

# Extract the top 10 largest clones for both samples
largest_clones_1 = head(sorted_sample_1, length(sorted_sample_1))
largest_clones_2 = head(sorted_sample_2, length(sorted_sample_2))

# Find shared clones between the two samples
shared_clones = intersect(names(largest_clones_1), names(largest_clones_2))

# Display the shared clones
shared_clones

# Create a table with common largest clones and their frequencies in both samples
shared_clones_table = data.frame(
  Clone = shared_clones,
  largest_clones_1 = as.numeric(largest_clones_1[shared_clones]),
  largest_clones_2 = as.numeric(largest_clones_2[shared_clones])
)
columns = c("Clone", sample_id_1, sample_id_2)


# Display the table of common largest clones
print(shared_clones_table)

colnames(shared_clones_table) = columns
print(shared_clones_table)
# Save the shared_clones_table to a CSV file
write.csv(shared_clones_table, file = "Results/shared_clones_table_X_113.csv", row.names = FALSE)


#-------------------------------------------------------------------------------


#-------------------Clonal overlay of all the samples --------------------------
fig.height <- 7
fig.width <- 10
clonalOverlay(merged_projected, 
              reduction = "umap", 
              cutpoint = 1, 
              bins = 20, 
              facet.by = "sample") + guides(color = "none") +
              ggtitle("Clonal Overlay of Expanded Cells All Sample") +
              scale_color_manual(values = ref.cd8@misc$atlas.palette)

fig_name = "plots/clonal_overlay_expanded_cells_all.png"
ggsave(fig_name, height = fig.height, width = fig.width, dpi = 300)

#-------------------------------------------------------------------------------

#-----Clonal Occupy/expansion by subtype All sample plot combined sample--------
fig.height = 7
fig.width = 16


lv_order = c("Hyperexpanded (100 < X <= 500)", "Large (20 < X <= 100)", 
             "Medium (5 < X <= 20)", "Small (1 < X <= 5)","Single (0 < X <= 1)")
seurat_projected_all$cloneSize = factor(seurat_projected_all$cloneSize, 
                                        levels = lv_order)

query_list = SplitObject(seurat_projected_all, split.by = "sample")

prop_metric = "count"
if (prop_metric == "count"){
  y_limit_range <- c(0, 2600)
  prop_flag = FALSE
} else {
  prop_flag = TRUE
  y_limit_range <- c(0, 1) #or change accordingly
  }
plots = lapply(names(query_list), function(x) {
  data_subset <- subset(seurat_projected_all, subset = sample == x)
  
  clonalOccupy(data_subset,
               x.axis = "functional.cluster", proportion=prop_flag,
               palette = "heat") + ggtitle(x) + ylim(y_limit_range) +
               theme(legend.position = "none", 
                     axis.text.x = element_text(angle = 45, hjust = 1)) +
               scale_x_discrete(drop = FALSE)
})
# Add legend to the first plot
plots[[1]] <- plots[[1]] + theme(legend.position = "left")
names(plots) <- names(query_list)
combined_plot <- wrap_plots(plots, ncol = length(names(query_list)), sharey = TRUE, sharex = TRUE)
combined_plot
fig_name = paste0("plots/occupied_repertorire_by_all_sample_functional_cluster_", prop_metric, ".png")
ggsave(fig_name, combined_plot, height = fig.height, width = fig.width, dpi = 300)

#-------------------------------------------------------------------------------

#------------------Clonal Occupy/expansion sample wise-------------------------- 
fig.height = 6
fig.width = 8
plot = clonalOccupy(seurat_projected_all,
             x.axis = "sample", proportion=TRUE,
             palette = "heat") + ggtitle("All sample Colonal Occupy/expansaion proportion") + 
             theme(legend.position = "right", 
             axis.text.x = element_text(angle = 45, hjust = 1)) 
plot
fig_name = paste("plots/occupied_repertorire_by_all_sample_proportion.png")
ggsave(fig_name, plot, height = fig.height, width = fig.width, dpi = 300)

#-------------------------------------------------------------------------------

#----------------Most Proliferative Clones for Each Patient---------------------

#Identify cycling cells
##Change the merged projected object here.
patient_id = "X_112"
merged_projected_patient = subset(seurat_projected_all, subset = patient == patient_id)
merged_projected_patient$is_cycling = ifelse((merged_projected_patient$cycling.score.G1_S > 0.1 |
                                              merged_projected_patient$cycling.score.G2_M > 0.1), 
                                              yes = "Proliferating", no = "Resting")
clonotypes = table(merged_projected_patient$cdr3s_pat)
expanded = names(clonotypes)[clonotypes>=2]
frequency_proliferating = sapply(expanded, function(x){
  sub = subset(merged_projected_patient[[]], subset = cdr3s_pat == x)
  sum(sub$is_cycling == "Proliferating")/nrow(sub) #changing ncol to nrow. does not make sense to divide by number of features
  #print(sum(sub$is_cycling == "Proliferating")/nrow(sub))
})

sorted <- sort(frequency_proliferating, decreasing = T)
top_proliferating = head(sorted, 6)
top_proliferating
#-------------------------------------------------------------------------------

#---------------- See frequency of top 6 clones---------------------------------
fig.height=6
fig.width=10

df <- data.frame(top_proliferating)
df$clone <- rownames(df)
ggplot(df, aes(x=reorder(clone, top_proliferating), y=top_proliferating)) +
      geom_bar(stat="identity", fill="black") +coord_flip() + theme_light() +
      ggtitle(sprintf("Most Proliferative in %s Patient", patient_id))
fig_name = paste0("plots/most_prolif_clones_patient_", patient_id, "_barplot.png")
ggsave(fig_name, height=fig.height, width=fig.width, dpi = 300)

#-------------------------------------------------------------------------------

#----------See Most Proliferative clones on the reference space (UMAP)----------
fig.height=6
fig.width=16

plots <- list()
for (ctype in names(top_proliferating)) {
  sample_id = substr(ctype, nchar(ctype) - 7, nchar(ctype))
  freq <- top_proliferating[ctype]
  title <- sprintf("%s %.0f %% proliferating", sample_id, 100*freq)
  cells_in_clone <- which(merged_projected_patient[["cdr3s_pat"]]==ctype)
  plots[[ctype]] <- plot.projection(ref.cd8, merged_projected_patient[,cells_in_clone],
                                    pointsize = 1, linesize = 0.3, raster=T) + 
    ggtitle(title) + NoLegend()
}
wrap_plots(plots, ncol = 3)
fig_name = paste("plots/most_prolif_clones_umap_", patient_id, "_patient.png")
ggsave(fig_name, height=fig.height, width=fig.width, dpi = 300)

#-------------------------------------------------------------------------------

#------------------Most Proliferative Clones for Each Sample--------------------
#Identify cycling cells
##Change the merged projected object here.

sample_id = "X_113_BL"
merged_projected_sample = subset(seurat_projected_all, subset = sample == sample_id)
merged_projected_sample$is_cycling = ifelse((merged_projected_sample$cycling.score.G1_S > 0.1 |
                                               merged_projected_sample$cycling.score.G2_M > 0.1), 
                                             yes = "Proliferating", no = "Resting")
clonotypes = table(merged_projected_sample$CTaa)
expanded = names(clonotypes)[clonotypes>=2]
frequency_proliferating = sapply(expanded, function(x){
  sub = subset(merged_projected_sample[[]], subset = CTaa == x)
  sum(sub$is_cycling == "Proliferating")/nrow(sub) #changing ncol to nrow. does not make sense to divide by number of features
  #print(sum(sub$is_cycling == "Proliferating")/nrow(sub))
})

sorted <- sort(frequency_proliferating, decreasing = T)
top_proliferating = head(sorted, 6)
top_proliferating

#See frequency of top 6 clones
fig.height=3
fig.width=8

df <- data.frame(top_proliferating)
df$clone <- rownames(df)
ggplot(df, aes(x=reorder(clone, top_proliferating), y=top_proliferating)) +
  geom_bar(stat="identity", fill="black") +coord_flip() + theme_light() +
  ggtitle(sample_id)
fig_name = paste0("plots/most_prolif_clones_sample_", sample_id, "_barplot.png")
ggsave(fig_name, height=fig.height, width=fig.width, dpi = 300)


#----------See Most Proliferative clones on the reference space (UMAP)----------
fig.height=6
fig.width=16

plots <- list()
for (ctype in names(top_proliferating)) {
  sample_id = substr(ctype, nchar(ctype) - 7, nchar(ctype))
  freq <- top_proliferating[ctype]
  title <- sprintf("%s %.0f %% proliferating", sample_id, 100*freq)
  cells_in_clone <- which(merged_projected_sample[["cdr3s_pat"]]==ctype)
  plots[[ctype]] <- plot.projection(ref.cd8, merged_projected_sample[,cells_in_clone],
                                    pointsize = 1, linesize = 0.3, raster=T) + 
    ggtitle(title) + NoLegend()
}
wrap_plots(plots, ncol = 3)
fig_name = paste("plots/most_prolif_clones_umap_", sample_id, "_sample.png")
ggsave(fig_name, height=fig.height, width=fig.width, dpi = 300)

#-------------------------------------------------------------------------------

#============== Clonal Sharing Between T-cell Subtypes==========================

#------------Clonal overlap between cell subtypes Each Patient -----------------
fig.height=4
fig.width=6

patient_id = "X_113"
method = "raw"
query = subset(merged_projected, subset = patient == patient_id)
clonalOverlap(query, 
              cloneCall = "CTaa", method = method) + ggtitle(patient_id) + 
              theme(axis.text.x = element_text(angle = 45, hjust = 1))
fig_name = paste0("plots/", method, "_clonal_overlap_heatmap_", patient_id, ".png")
ggsave(fig_name, height=fig.height, width=fig.width, dpi = 300)

#-------------------------------------------------------------------------------

#----------------------------Clonal overlap Each Sample ------------------------
fig.height=4
fig.width=6

sample_id = "X_112_BL"
method = "raw"
query = subset(merged_projected, subset = sample == sample_id)
clonalOverlap(query, cloneCall = "CTaa", method = method) + ggtitle(sample_id) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
fig_name = paste0("plots/", method, "_clonal_overlap_heatmap_sample_", sample_id, ".png")
ggsave(fig_name, height=fig.height, width=fig.width, dpi = 300)

#-------------------------------------------------------------------------------

#---------------Subset connectivity using circos plots Each Patient-------------
patient_id = "X_112"
query = subset(merged_projected, subset = patient == patient_id)
circles <- getCirclize(query, 
                       cloneCall = "CTaa",
                       group.by = "functional.cluster",
                       proportion = FALSE
)

#Just assigning the normal colors to each cluster
grid.cols <- hue_pal()(length(unique(merged_projected$functional.cluster)))
names(grid.cols) <- unique(merged_projected$functional.cluster)

#Graphing the chord diagram
fig_name = paste0("plots/circos_subtype_", patient_id, ".png")
pdf(fig_name, width = 5, height = 5)
chordDiagram(circles, self.link = 1, grid.col = grid.cols)
dev.off()

### Get only largest clones
large <- subset(query, subset = clonalFrequency >= 5)
circles <- getCirclize(large,
                       cloneCall = "cdr3s_pat",
                       group.by = "functional.cluster",
                       proportion = FALSE
)

#Just assigning the normal colors to each cluster
grid.cols <- hue_pal()(length(unique(merged_projected$functional.cluster)))
names(grid.cols) <- unique(merged_projected$functional.cluster)

fig_name = paste0("plots/circos_subtype_largest_", patient_id, ".png")
pdf(fig_name, width = 5, height = 5)
circlize::chordDiagram(circles,
                       self.link = 1, 
                       grid.col = grid.cols)
dev.off()

#===============================================================================
#-----------------------------Quantifying Clonal Bias Each Patient--------------
##clonalBias
patient_id = "X_113"

fig.height=6
fig.width=10

query = subset(merged_projected, subset = patient == patient_id)
query$cloneSize = factor(query$cloneSize, levels = lv_order)


clonalBias(query, 
           cloneCall = "CTaa", 
           split.by = "patient", 
           group.by = "functional.cluster",
           n.boots = 20, 
           min.expand = 10) +
           ggtitle(paste0("Patient: ", patient_id))


fig_name = paste0("plots/clonotype_bias_patient_", patient_id, ".png")
ggsave(fig_name, height=fig.height, width=fig.width, dpi = 300)

### Getting most biased clones out and plot into the reference map

large = subset(query, subset = clonalFrequency >= 5)

# large$cloneSize = large$clonalFrequency
biased <- clonalBias(large, cloneCall = "CTaa", split.by = "patient",
                     group.by = "functional.cluster", min.expand = 10, exportTable = TRUE)
most.biased <- biased[order(biased$Z.score, decreasing = TRUE),]
head(most.biased)

##### See where the most "biased" clonotypes are found on the reference.

fig.height=8
fig.width=16

plots <- list()
for (i in 1:6) {
  ctype <- most.biased[i, "Clone"]
  patient <- most.biased[i, "Sample"]
  print(ctype)
  size <- most.biased[i, "ncells"]
  zscore <- most.biased[i, "Z.score"]
  cells_in_clone <- which(merged_projected[["CTaa"]]==ctype)
  title <- sprintf("Clone %s, size %s, Z=%.1f", i, size, zscore)
  plots[[i]] <- plot.projection(ref.cd8, merged_projected[,cells_in_clone],
                                pointsize = 1, linesize = 0.3, raster=T) + 
    ggtitle(title) + theme(legend.position = "none")
}
plots[[1]] <- plots[[1]] + theme(legend.position = "left")  # Add legend to the first plot
wrap_plots(plots, ncol = 3)
fig_name = paste0("plots/biased_clones_umap_", patient_id, ".png")
ggsave(fig_name, height=fig.height, width=fig.width, dpi = 300)



#-------------------------------------------------------------------------------
#-----------------------------Quantifying Clonal Bias Each Sample---------------
##clonalBias
sample_id = "X_113_BL"

fig.height=6
fig.width=10

query = subset(merged_projected, subset = sample == sample_id)
min_clonal_freq = 10

# query$cloneSize = factor(query$cloneSize, levels = lv_order)
# clonalBias(query, 
#            cloneCall = "CTaa", 
#            split.by = "sample", 
#            group.by = "functional.cluster",
#            n.boots = 20, 
#            min.expand = min_clonal_freq) + 
#            scale_fill_manual(values = ref.cd8@misc$atlas.palette) +
#           ggtitle(paste0("Sample: ", sample_id)) + xlab("Clonal Size (N cells)") + 
#   ylab("Clonal Bias") + labs(size = "Clone Size/Frequency", fill = "Top State") 
#   

biased <- clonalBias(query, cloneCall = "CTaa",
                     split.by = "sample",
                     group.by = "functional.cluster",
                     min.expand = min_clonal_freq,
                     exportTable = TRUE)

# Convert the cloneSize to a factor with the desired labels
biased$dotSize <- factor(biased$cloneSize,
                         levels = rev(lv_order))

# # Now, plot using the dotSize as the categorical variable for size
ggplot(biased, aes(x=ncells, y=bias)) +
  geom_point(aes(fill=factor(Top_state), size = cloneSize),
             shape = 21, stroke = 0.25) +
  scale_size_manual(values = c(3, 5, 7, 9, 11)) +  # Adjust size values manually for each category
  scale_fill_manual(values = ref.cd8@misc$atlas.palette) +
  theme_classic() + ggtitle(paste0("Sample: ", sample_id)) +
  labs(size = "Clone Size/Frequency", fill = "Top State") +
  xlab("Clonal Size (N cells)") +
  ylab("Clonal Bias")

#CAVNDIRTYKYIF_CASRHLGVGPYNEQFF 331 this is hyperexpanded except in 1 of the occurance

fig_name = paste0("plots/clonotype_bias_sample_", sample_id, "_v2.png")
ggsave(fig_name, height=fig.height, width=fig.width, dpi = 300)

### Getting most biased clones out and plot into the reference map

large = subset(query, subset = clonalFrequency >= 5)

biased <- clonalBias(large, cloneCall = "CTaa", 
                     split.by = "sample",
                     group.by = "functional.cluster",
                     min.expand = min_clonal_freq, 
                     exportTable = TRUE)
most.biased <- biased[order(biased$Z.score, decreasing = TRUE),]
head(most.biased)

##### See where the most "biased" clonotypes are found on the reference.

fig.height=8
fig.width=16

plots <- list()
for (i in 1:6) {
  ctype <- most.biased[i, "Clone"]
  patient <- most.biased[i, "Sample"]
  print(ctype)
  size <- most.biased[i, "ncells"]
  zscore <- most.biased[i, "Z.score"]
  cells_in_clone <- which(merged_projected[["CTaa"]]==ctype)
  title <- sprintf("Clone %s, size %s, Z=%.1f", i, size, zscore)
  plots[[i]] <- plot.projection(ref.cd8, merged_projected[,cells_in_clone],
                                pointsize = 1, linesize = 0.3, raster=T) + 
    ggtitle(title) + theme(legend.position = "none")
}
plots[[1]] <- plots[[1]] + theme(legend.position = "left")  # Add legend to the first plot
wrap_plots(plots, ncol = 3)
fig_name = paste0("plots/biased_clones_umap_sample_", sample_id, ".png")
ggsave(fig_name, height=fig.height, width=fig.width, dpi = 300)
#===============================================================================

#clonal Network and Shared Clones


# patient_id = "X_113"
sample_id = "X_112_9M"
#Levels: CD8.CM CD8.EM CD8.MAIT CD8.NaiveLike CD8.TEMRA CD8.TEX

fig.height=6
fig.width=10
clonecall_var = "CTaa"
filter_identity = NULL
# patient_id = 'all'
# query = subset(merged_projected, subset = patient == patient_id)
query = subset(merged_projected, subset = sample == sample_id)


clonalNetwork(query, 
              reduction = "umap", 
              group.by = "functional.cluster",
              filter.clones = NULL,
              filter.identity = filter_identity,
              cloneCall = "CTaa")+
              scale_color_manual(values = ref.cd8@misc$atlas.palette)



# fig_name = paste0("plots/clonal_network_", patient_id, ".pdf")

fig_name = paste0("plots/clonal_network_", sample_id, ".png")

ggsave(fig_name, height=fig.height, width=fig.width, dpi = 300)

## Shared Clones


fig.height=6
fig.width=16
clonecall_var = "CTaa"
sample_id = "X_112_BL"

patient_id = "X_113"
#Levels: CD8.CM CD8.EM CD8.MAIT CD8.NaiveLike CD8.TEMRA CD8.TEX

# query = subset(merged_projected, subset = sample == sample_id)


query = subset(merged_projected, subset = patient == patient_id)

#query = merged_projected

shared_clones = clonalNetwork(query, reduction = "umap", 
                                  group.by = "functional.cluster",
                                  cloneCall = clonecall_var, exportClones = TRUE)

first_6_clones <- c(shared_clones$clone[1:10])
other_clones <- c(shared_clones$clone[11:length(shared_clones$clone)])
# Print the result
print(first_6_clones)
###### Highlight Clones
query <- highlightClones(query, 
                         cloneCall= clonecall_var, 
                         sequence = first_6_clones)

Seurat::DimPlot(query, group.by = "highlight", split.by = 'sample',
                sizes.highlight = 1,
                pt.size = .5,
                alpha = 1,
                stroke.size = 1) + 
  theme(plot.title = element_blank()) +  scale_color_viridis(discrete = TRUE)

fig_name = paste0("plots/shared_clone_highlight_", patient_id, "_v2.png")
# fig_name = paste0("plots/shared_clone_highlight_", sample_id, "_v2.png")
# fig_name = paste0("plots/shared_clone_highlight_all_sample.png")
ggsave(fig_name, height=fig.height, width=fig.width, dpi = 300)


#============================================End Clonal Analysis=============================================
