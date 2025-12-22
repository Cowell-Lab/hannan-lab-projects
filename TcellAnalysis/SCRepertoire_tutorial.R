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
library(scater)
library(circlize)
library(scales)
#BiocManager::install("scRepertoire", force = TRUE)
# BiocManager::install("scater", force = TRUE)

#===============================Data Preprocessing==============================
datadir="data/"
cached_seurat <- sprintf("%s/HannanLab_alldata_seurat.rds", datadir)
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
    seurat_obj <- CreateSeuratObject(counts = data, project = "HannanLab", 
                                     min.cells = 3, min.features = 100)
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
merged_seurat <- NormalizeData(merged_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

###Read TCR Informations
S1 = read.csv("data/X_112_BL/filtered_contig_annotations.csv", sep = ",")
S2 = read.csv("data/X_112_9M/filtered_contig_annotations.csv", sep = ",")
S3 = read.csv("data/X_113_BL/filtered_contig_annotations.csv", sep = ",")
S4 = read.csv("data/X_113_9M/filtered_contig_annotations.csv", sep = ",")
contig_list <- list(S1, S2, S3, S4)
head(contig_list[[1]])

### Combining Contigs into Clones
combined_TCR <- combineTCR(contig_list, samples = c("X_112_BL", "X_112_9M", 
                                                    "X_113_BL", "X_113_9M"))
### addVariable
combined_TCR <- addVariable(combined_TCR, 
                            variable.name = "time_frame", 
                            variables = rep(c("BL", "9M"), 2))
combined_TCR = lapply(combined_TCR, function(x){
  x$cdr3s_pat = paste(x$CTaa, x$sample, sep = "_"); x})

head(combined_TCR[[1]])

# Extract cell barcodes from Seurat object
seurat_barcodes <- colnames(merged_seurat)

# Extract cell barcodes from TCR data (assuming `combined_TCR` contains cell barcodes)
# This depends on the structure of your TCR data, let's assume the barcode column is named "barcode"
# If your TCR data has different naming, replace 'barcode' with the correct column name
tcr_barcodes <- unique(unlist(lapply(combined_TCR, function(x) x$barcode)))

# Subset the Seurat object to retain only cells present in both datasets
merged_seurat_filtered <- merged_seurat[, seurat_barcodes %in% tcr_barcodes]

#Combine TCR data and Seurat Expression data
seurat_sc_combined <- combineExpression(combined_TCR, 
                                        merged_seurat_filtered, 
                                        cloneCall="gene", 
                                        proportion = FALSE, 
                                        cloneSize=c(Single=1, Small=5, Medium=20,
                                        Large=100, Hyperexpanded=500))

seurat_sc_combined@meta.data$cdr3s_pat = paste(seurat_sc_combined@meta.data$CTaa,
                                                seurat_sc_combined@meta.data$sample, 
                                                sep = "_")
seurat_sc_combined@meta.data$cdr3s_pat[is.na(seurat_sc_combined@meta.data$CTaa)] <- NA

#remove files after combining
rm(data)
rm(seurat_objects)
rm(seurat_obj)
rm(S1)
rm(S2)
rm(S3)
rm(S4)
rm(contig_list)

#===========================================ProjectIILS=====================================================
#### Load reference CD8 T map
ref.file <- "CD8T_human_ref_v1.rds"
if (!file.exists(ref.file)) {
  download.file("https://figshare.com/ndownloader/files/38921366", 
                destfile = ref.file)
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
ggsave("plots/CD8T_reference_umap.pdf", height = 6, width = 8)

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



# Subset merged Seurat object based on treatment
merged_seurat_BL <- subset(seurat_sc_combined, subset = treatment == "BL")
merged_seurat_9M <- subset(seurat_sc_combined, subset = treatment == "9M")

#Split by sample
seurat_list_BL = SplitObject(merged_seurat_BL, split.by = "sample")
seurat_list_9M = SplitObject(merged_seurat_9M, split.by = "sample")
seurat_list_all = SplitObject(seurat_sc_combined, split.by = "sample")

seurat_projected_BL = Run.ProjecTILs(seurat_list_BL, ref = ref.cd8, ncores = 6)
seurat_projected_9M = Run.ProjecTILs(seurat_list_9M, ref = ref.cd8, ncores = 6)
seurat_projected_all = Run.ProjecTILs(seurat_sc_combined, ref = ref.cd8, split.by = "sample", ncores = 6)

#### Compare Marker gene expression profiles of query data with reference map

fig.height <- 10
fig.width <- 20
query_sample <- "all"
query_sample <- "X_112_9M"

genes4radar <- c("CD4", "CD8A", "TCF7", "CCR7", "IL7R", "LMNA", "GZMA", "GZMK", "FGFBP2",
                 "XCL1", "CRTAM", "TOX", "PDCD1", "HAVCR2", "PRF1", "GZMB", "TRAV1-2",
                 "KLRB1")
if(query_sample == "all"){
  plot.states.radar(ref.cd8, seurat_projected_all, genes4radar = genes4radar)
}else{
  plot.states.radar(ref.cd8, seurat_projected_all[[query_sample]], genes4radar = genes4radar)
}

plot_filename <- paste("plots/radars_projected_", query_sample, ".pdf", sep = "")

ggsave(plot_filename, height = 11, width = 16)


### Look into the radar projection treatment wise
query_list = SplitObject(seurat_projected_all, split.by = "sample")
plot.states.radar(ref.cd8, query_list[['X_112_9M']], genes4radar = genes4radar)
plot_filename <- paste("plots/radars_projected_treatment_all.pdf", sep = "")

ggsave(plot_filename, height = 11, width = 16)
### Cell Subtype Composition of query data

query_sample <- "X_112_BL"

a <- plot.projection(ref.cd8, query = seurat_projected_all[[query_sample]],
                     pointsize = 0.1, linesize = 0.3, raster=T) + NoLegend() +
  ggtitle(query_sample)

b <- plot.statepred.composition(ref.cd8, query = seurat_projected_all[[query_sample]]) 
a | b
plot_title <- paste("plots/single_patient_composition_", query_sample, ".pdf")
ggsave(plot_title, height=6, width=10)



# ######Combine all samples
# fig.height <- 7
# fig.width <- 14
# 
# #for large enough sample we can compare their composition
# plots = lapply(names(seurat_projected), function(x){
#   plot.projection(ref.cd8, query = seurat_projected[[x]],pointsize = 0.1,
#                   linesize = 0.3, raster=T) + NoLegend() + ggtitle(x) + theme(legend.position = "none")
# })
# plots[[1]] <- plots[[1]] + theme(legend.position = "left")  # Add legend to the first plot
# 
# names(plots) <- names(seurat_projected)
# # Combine and display the plots in a single column
# combined_plot <- wrap_plots(plots, ncol = length(names(seurat_projected)))
# 
# # Display the combined plot
# combined_plot
# ggsave("plots2/all_patient_composition_projection_plot.pdf", combined_plot, height = fig.height, width = fig.width)

### Exclude smal Samples and compare their compositions in terms of cell subtypes
sizes = as.vector(lapply(seurat_projected_all, ncol))
keep = names(sizes)[sizes > 100]
seurat_projected_all = seurat_projected_all[keep]

fig.height <- 7
fig.width <- 14
y_limit_range <- c(0, 3000) 

labels <- ref.cd8[["functional.cluster"]][,1]

states_all <- levels(factor(labels))
tb <- table(factor(seurat_projected_all[["X_112_BL"]][["functional.cluster"]][,1], levels=states_all))


# for large enough sample we can compare their composition
plots = lapply(names(seurat_projected_all), function(x) {
  plot.statepred.composition(ref.cd8, query = seurat_projected_all[[x]],
                             metric = "Count") + ggtitle(x) +
                             theme(legend.position = "none")+ ylim(y_limit_range)
                              })
plots[[1]] <- plots[[1]] + theme(legend.position = "left")  # Add legend to the first plot
names(plots) <- names(seurat_projected_all)
# Combine and display the plots in a single column with shared y-axis
combined_plot <- wrap_plots(plots, ncol = length(names(seurat_projected_all)), sharey = TRUE)
# Display the combined plot
combined_plot
ggsave("plots/all_patient_composition_counts.pdf", combined_plot, height = fig.height, width = fig.width)

### Merge lists of object to obtain a single object
merged_projected = Reduce(merge.Seurat.embeddings, seurat_projected_all)
Idents(merged_projected) = "functional.cluster"


### Merge lists of object to obtain a single object
merged_projected_BL = Reduce(merge.Seurat.embeddings, seurat_projected_BL)
Idents(merged_projected_BL) = "functional.cluster"
### Merge lists of object to obtain a single object
merged_projected_9M = Reduce(merge.Seurat.embeddings, seurat_projected_9M)
Idents(merged_projected_9M) = "functional.cluster"


merged_projected
merged_projected_BL
merged_projected_9M
#============================================Clonal Analysis=============================================
###Identify mmost expanded clones
freqs_BL = lapply(seurat_projected_BL, function(x){
  table(x$cdr3s_pat)/sum(!is.na(x$cdr3s_pat))
})

freqs_9M = lapply(seurat_projected_9M, function(x){
  table(x$cdr3s_pat)/sum(!is.na(x$cdr3s_pat))
})

freqs_all = lapply(seurat_projected_all, function(x){
  table(x$cdr3s_pat)/sum(!is.na(x$cdr3s_pat))
})

freqs_BL = Reduce(c, freqs_BL)
freqs_9M = Reduce(c, freqs_9M)
freqs_all = Reduce(c, freqs_all)


sorted_BL = sort(freqs_BL, decreasing = TRUE)
sorted_9M = sort(freqs_9M, decreasing = TRUE)
sorted_all = sort(freqs_all, decreasing = TRUE)
largest_clones_BL = head(sorted_BL, 6)
largest_clones_9M = head(sorted_9M, 6)
largest_clones_all = head(sorted_all, 6)


fig.height <- 6
fig.width <- 12
plots = list()

largest_clones = largest_clones_all
for(i in 1:length(largest_clones)){
  ctype = names(largest_clones)[i]
  patient <- substr(ctype, nchar(ctype) - 7, nchar(ctype))
  cells = which(merged_projected[["cdr3s_pat"]] == ctype)
  plots[[i]] = plot.projection(ref.cd8, merged_projected[, cells],
                               pointsize = 1, linesize = 0.3, raster=T) + 
      ggtitle(sprintf("%s - Freq %.1f %%", patient, 100*largest_clones[i])) + theme(legend.position = "none")
}
plots[[1]] <- plots[[1]] + theme(legend.position = "left")  # Add legend to the first plot
combined_plot = wrap_plots(plots, ncol = 3)
combined_plot
ggsave("plots/largest_clones_among_all_samples.pdf", combined_plot, height = fig.height, width = fig.width)



#### Clonal Overlay by sample
fig.height <- 6
fig.width <- 10
plot = clonalOverlay(merged_projected, 
              reduction = "umap", 
              cutpoint = 1, 
              bins = 15, 
              facet.by = "sample") + guides(color = "none") +
              ggtitle("Clonal Overlay of Expanded Cells All Sample")
plot
ggsave("plots/clonal_overlay_expanded_cells_all.pdf", plot, height = fig.height, width = fig.width)


#### Shared Clones
shared_clones_BL <- clonalNetwork(merged_projected_BL, reduction = "umap", 
                               group.by = "functional.cluster",
                               cloneCall = "CTaa",  exportClones = TRUE)


#### Shared Clones
shared_clones_9M <- clonalNetwork(merged_projected_9M, reduction = "umap", 
                               group.by = "functional.cluster", 
                               cloneCall = "CTaa", exportClones = TRUE)

#### Shared Clones
shared_clones_all <- clonalNetwork(merged_projected,  reduction = "umap", 
                               group.by = "functional.cluster",
                               cloneCall = "CTaa",  exportClones = TRUE)
head(shared_clones_all)
head(shared_clones_BL)
head(shared_clones_9M)


###### Highlight Clones
merged_projected <- highlightClones(merged_projected, 
                                 cloneCall= "CTaa", 
                                 sequence = c("CAVRGGNTGKLIF;CAVGALGGSARQLTF_CASSNRQGVRLYEQYF", 
                                              "CGTDPTATDKLIF_CASSSTGYTNEKLFF",
                                              "CAVNDIRTYKYIF_CASRHLGVGPYNEQFF",
                                              "CAVPSGGSYIPTF;CVVEGSSNDYKLSF_CSAPAAGGLETQYF"
                                              ))


fig.height <- 3
fig.width <- 10
plot = Seurat::DimPlot(merged_projected, group.by = "highlight") + 
  ggplot2::theme(plot.title = element_blank()) + ggtitle("Highlighted Shared Clones among all Samples")

plot
ggsave("plots/shared_clone_highlight_all.pdf", plot, height = fig.height, width = fig.width)


###### Highlight Clones BL
merged_projected_BL <- highlightClones(merged_projected_BL, 
                                    cloneCall= "CTaa", 
                                    sequence = c("CGTDPTATDKLIF_CASSSTGYTNEKLFF", 
                                                 "CAVPSGGSYIPTF;CVVEGSSNDYKLSF_CSAPAAGGLETQYF",
                                                 "CAMREGRRNTGFQKLVF_CASTRGGSSYEQYF",
                                                 "CAAKRSGNQFYF_CASSQEIPGGIGGNEKLFF"
                                    ))
fig.height <- 3
fig.width <- 10
plot = Seurat::DimPlot(merged_projected_BL, group.by = "highlight") + 
  ggplot2::theme(plot.title = element_blank()) + ggtitle("Highlighted Shared Clones among BL Samples")

plot
ggsave("plots/shared_clone_highlight_BL.pdf", plot, height = fig.height, width = fig.width)

###### Highlight Clones 9M
merged_projected_9M <- highlightClones(merged_projected_9M, 
                                       cloneCall= "CTaa", 
                                       sequence = c("CAVRGGNTGKLIF;CAVGALGGSARQLTF_CASSNRQGVRLYEQYF", 
                                                    "CAVNDIRTYKYIF_CASRHLGVGPYNEQFF",
                                                    "CGTDPTATDKLIF_CASSSTGYTNEKLFF",
                                                    "CAMREGRRNTGFQKLVF_CASTRGGSSYEQYF"
                                       ))
fig.height <- 3
fig.width <- 10
plot = Seurat::DimPlot(merged_projected_9M, group.by = "highlight") + 
  ggplot2::theme(plot.title = element_blank()) + 
  ggtitle("Highlighted Shared Clones among 9M Samples")

plot
ggsave("plots/shared_clone_highlight_9M.pdf", plot, height = fig.height, width = fig.width)


### Clonal Occupy/expansion by T cell subtype

fig.height <- 6
fig.width <- 7

plot = clonalOccupy(merged_projected, 
            x.axis = "functional.cluster")+ 
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
            ggtitle("Occupied Repertoire by All Sample Sample Functional Cluster")
plot
ggsave("plots/occupied_repertorire_by_all_sample_functional_cluster.pdf",
       plot, height = fig.height, width = fig.width)


### Clonotype Proliferation Rate
#Identify cycling cells
##Change the merged projected object here.
BL_merged_projected = merged_projected_9M

BL_merged_projected$is_cycling = ifelse((BL_merged_projected$cycling.score.G1_S > 0.1 |
                                           BL_merged_projected$cycling.score.G2_M > 0.1), yes = "Proliferating",
                                     no = "Resting")
clonotypes = table(BL_merged_projected$cdr3s_pat)
expanded = names(clonotypes)[clonotypes>=2]
frequency_proliferating = sapply(expanded, function(x){
  sub = subset(BL_merged_projected[[]], subset = cdr3s_pat == x)
  sum(sub$is_cycling == "Proliferating")/ncol(sub)
})

sorted <- sort(frequency_proliferating, decreasing = T)
top_proliferating = head(sorted, 6)

#See frequency of top 10 clones
fig.height=3
fig.width=8

df <- data.frame(top_proliferating)
df$clone <- rownames(df)
ggplot(df, aes(x=reorder(clone, top_proliferating), y=top_proliferating)) +
  geom_bar(stat="identity", fill="black") + coord_flip() + theme_light()+ggtitle("Most Proliferative in 9M Sample")
ggsave("plots/most_prolif_clones_9M_barplot.pdf", height=fig.height, width=fig.width)


##### See Most Proliferative clones on the reference space
fig.width=16
fig.height=8

plots <- list()
for (ctype in names(top_proliferating)) {
  
  freq <- top_proliferating[ctype]
  title <- sprintf("%.0f %% proliferating", 100*freq)
  cells_in_clone <- which(BL_merged_projected[["cdr3s_pat"]]==ctype)
  
  plots[[ctype]] <- plot.projection(ref.cd8, BL_merged_projected[,cells_in_clone],
                                    pointsize = 1, linesize = 0.3, raster=T) + 
    ggtitle(title) + NoLegend()
}

wrap_plots(plots, ncol = 3)
ggsave("plots/most_prolif_clones_umap_9M.pdf", height=fig.height, width=fig.width)


##### Clonal Sharing Between T-cell Subtypes
### Clonal overlap
query = merged_projected_9M

plot = clonalOverlap(query, 
              cloneCall = "cdr3s_pat", 
              method = "cosine") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
plot
ggsave("plots/cosine_heatmap_9M.pdf", height=6, width=6)



### Subset connectivity using circos plots
query = merged_projected_9M
circles <- getCirclize(query, 
                       group.by = "functional.cluster",
                       proportion = TRUE
                       )

#Just assigning the normal colors to each cluster
grid.cols <- hue_pal()(length(unique(merged_projected$functional.cluster)))
names(grid.cols) <- unique(merged_projected$functional.cluster)

#Graphing the chord diagram
pdf("plots/circos_subtype_9M.pdf", width = 5, height = 5)
chordDiagram(circles, self.link = 1, grid.col = grid.cols)
dev.off()

### Get only largest clones

large <- subset(merged_projected, subset = cloneSize >= 5)

circles <- getCirclize(large,
                       cloneCall = "cdr3s_pat",
                       group.by = "functional.cluster",
                       proportion = TRUE
                       )

#Just assigning the normal colors to each cluster
grid.cols <- hue_pal()(length(unique(merged_projected$functional.cluster)))
names(grid.cols) <- unique(merged_projected$functional.cluster)

pdf("plots2/circos_subtype_all_largest.pdf", width = 5, height = 5)
circlize::chordDiagram(circles,
                       self.link = 1, 
                       grid.col = grid.cols)
dev.off()

###Quantifying Clonal Bias
##clonalBias
clonalBias(merged_projected, 
           cloneCall = "cdr3s_pat", 
           split.by = "patient",
           group.by = "functional.cluster",
           n.boots = 10, 
           min.expand =0) + 
           scale_color_manual(values = ref.cd8@misc$atlas.palette)


ggsave("plots/clonotype_bias_all.pdf", height=3, width=5)

large <- subset(merged_projected_9M, subset = cloneSize >= 5)

biased <- clonalBias(large, cloneCall = "cdr3s_pat", split.by = "patient",
                        group.by = "functional.cluster", min.expand = 10, exportTable = TRUE)

most.biased <- biased[order(biased$Z.score, decreasing = TRUE),]
head(most.biased)

##### See where the most "biased" clonotypes are found on the reference.
fig.width=16
fig.height=8


plots <- list()
for (i in 1:6) {
  
  ctype <- most.biased[i, "Clone"]
  patient <- most.biased[i, "sample"]
  size <- most.biased[i, "ncells"]
  zscore <- most.biased[i, "Z.score"]
  cells_in_clone <- which(merged_projected[["cdr3s_pat"]]==ctype)
  title <- sprintf("Clone %s, size %s, Z=%.1f", i, size, zscore)
  
  plots[[i]] <- plot.projection(ref.cd8, merged_projected[,cells_in_clone],
                                pointsize = 1, linesize = 0.3, raster=T) + 
    ggtitle(title) + NoLegend()
}

wrap_plots(plots, ncol = 3)

ggsave("plots/biased_clones_umap_9M.pdf", height=7, width=12)

#============================Basic Clonal Visualizations==================================



