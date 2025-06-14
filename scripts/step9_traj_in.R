log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

message("CONFIGURATION STEP")
# A. Parameters: 
# 1. Load libraries. 
suppressMessages(library("GenomeInfoDbData"))
suppressMessages(library("Seurat"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("slingshot"))
suppressMessages(library("SingleCellExperiment"))
suppressMessages(library("rmarkdown"))
suppressMessages(library("gam"))
suppressMessages(library("pheatmap"))
suppressMessages(library("openxlsx"))
if (snakemake@params[["graphics"]]) {suppressMessages(library("rgl"))}
if (snakemake@params[["graphics"]]) {suppressMessages(library("htmltools"))}
message("1. Libraries were loaded.")

# 2. Folder configuration. 
dir.name = snakemake@params[["output_dir"]]
input_data = snakemake@input[["seurat_obj"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_scRepertoire", "5_degs", "6_annotation", "7_gs", "8_traj_in", "9_func_analysis", "10_RNAvelocity")
message("2. Folder paths were set.")

# 3. Get variables from Snakemake.  
selected_res = snakemake@params[["selected_res"]]
start.clus = snakemake@params[["start_clus"]]
end.clus = snakemake@params[["end_clus"]]
n_var_genes = snakemake@params[["n_var_genes"]]
n_plotted_genes = snakemake@params[["n_plotted_genes"]]
random_seed = snakemake@params[["random_seed"]]
pc = snakemake@params[["pc"]]
graphics = snakemake@params[["graphics"]]
ram = snakemake@resources[["mem_mb"]]
threads = snakemake@threads
message("3. Parameters were loaded.")

# 4. Analysis configuration. 
# RAM configuration.
options(future.globals.maxSize = ram*1024^2)
# Set seed.
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}
message(paste0("4. Seed was set at ", random_seed, "."))
message("Configuration finished.")
message("\n")

message("PROCESSING STEP")
# 10. Trajectory inference analysis using slingshot.
# 10.1. Seurat object is converted to SingleCellExperiment object (required by slingshot) & cluster resolution is selected. 
seurat <- readRDS(input_data)
assay_type <- seurat@active.assay
cluster_res <- paste0(assay_type,"_snn_res.",selected_res)
if (!(cluster_res %in% colnames(seurat@meta.data))){
  stop("Specified resolution is not available.")
}
seurat.sim <- as.SingleCellExperiment(seurat)
message("1. Seurat object was loaded.")

if (graphics){
  # Running UMAP to get 3 dimensions.
  seurat3D <- RunUMAP(seurat,dims = 1:pc, n.components = 3, verbose = FALSE)
  seurat.sim <- as.SingleCellExperiment(seurat)
  seurat3D.sim <- as.SingleCellExperiment(seurat3D)
  message("1.5. 3D UMAP was calculated.")
}

# 10.2. Slingshot algorithm (dimensions = UMAP). There are 4 options depending of the start and end cluster in the following trajectory:
if(is.numeric(start.clus) == FALSE && is.numeric(end.clus) == FALSE){
  seurat.sim <- slingshot(seurat.sim, clusterLabels=cluster_res, reducedDim="UMAP")
  if (graphics) {seurat3D.sim <- slingshot(seurat3D.sim, clusterLabels=cluster_res, reducedDim="UMAP")}
  const = FALSE
} else if(is.numeric(start.clus) == TRUE && is.numeric(end.clus) == FALSE){
  seurat.sim <- slingshot(seurat.sim, clusterLabels=cluster_res, reducedDim="UMAP", start.clus = start.clus)
  if (graphics) {seurat3D.sim <- slingshot(seurat3D.sim, clusterLabels=cluster_res, reducedDim="UMAP", start.clus = start.clus)}
  const = TRUE
} else if(is.numeric(start.clus) == FALSE && is.numeric(end.clus) == TRUE){
  seurat.sim <- slingshot(seurat.sim, clusterLabels=cluster_res, reducedDim="UMAP", end.clus = end.clus)
  if (graphics) {seurat3D.sim <- slingshot(seurat3D.sim, clusterLabels=cluster_res, reducedDim="UMAP", end.clus = end.clus)}
  const = TRUE
} else if(is.numeric(start.clus) == TRUE && is.numeric(end.clus) == TRUE){
  seurat.sim <- slingshot(seurat.sim, clusterLabels=cluster_res, reducedDim="UMAP", start.clus = start.clus, end.clus=end.clus)
  if (graphics) {seurat3D.sim <- slingshot(seurat3D.sim, clusterLabels=cluster_res, reducedDim="UMAP", start.clus = start.clus, end.clus=end.clus)}
  const = TRUE
}
message("2. Slingshot trajectories were calculated.")

# 10.3. Plots generation.
# 10.3.1. Set the number of colours from the palette and extend it.
n_col <- length(levels(seurat.sim@colData[, cluster_res]))
getPalette <- colorRampPalette(brewer.pal(9,'Set1'))

if (graphics) {
  # 10.3.2. Set a good window size for the 3D plots.
  r3dDefaults$windowRect <- c(0,50, 1024, 720) 
  
  # 10.3.3. Curves 3D plot with legend --> HTML output.
  plot3d(reducedDims(seurat3D.sim)$UMAP[,1:3], col = getPalette(n_col)[seurat3D.sim@colData[, cluster_res]])
  plot3d.SlingshotDataSet(SlingshotDataSet(seurat3D.sim), type = "curves", add = TRUE, lwd =3)
  legend3d("topright", legend=paste0("Cluster - ", levels(seurat3D.sim@colData[, cluster_res])), pch=16, col=getPalette(n_col), inset=c(0))
  p_curves <- rglwidget(width = 1240, height = 1024)
  htmltools::save_html(htmltools::tagList(p_curves), file = paste0(dir.name, "/", folders[8], "/", paste0("3D_curves_", selected_res, "_res.html")))

  # 10.3.4. Lineage 3D plot with legend --> HTML output.
  plot3d(reducedDims(seurat3D.sim)$UMAP[,1:3], col = getPalette(n_col)[seurat3D.sim@colData[, cluster_res]])
  plot3d.SlingshotDataSet(SlingshotDataSet(seurat3D.sim), type = "lineages", add = TRUE, lwd =3)
  legend3d("topright", legend=paste0("Cluster - ", levels(seurat3D.sim@colData[, cluster_res])), pch=16, col=getPalette(n_col), inset=c(0))
  p_lineages <- rglwidget(width = 1240, height = 1024)
  htmltools::save_html(htmltools::tagList(p_lineages), file = paste0(dir.name, "/", folders[8], "/", paste0("3D_lineages_", selected_res, "_res.html")))
}
# 10.3.5. Curves 2D plot with legend --> pdf output. 
pdf(paste0(dir.name, "/", folders[8], "/2D_curves_", selected_res, "_res.pdf"), width = 11, height = 11)
plot(reducedDims(seurat.sim)$UMAP, col = getPalette(n_col)[seurat.sim@colData[, cluster_res]],
     pch=16, asp = 1, main = "2D curves - Clusters Trajectories")
legend("topright", legend=paste0("Cluster - ", levels(seurat.sim@colData[, cluster_res])), pch=16, col=getPalette(n_col))
lines(SlingshotDataSet(seurat.sim), lwd=2, col = 'black')
dev.off()

# 10.3.6. Lineage 2D plot with legend --> pdf output.
pdf(paste0(dir.name, "/", folders[8], "/2D_lineages_", selected_res, "_res.pdf"), width = 11, height = 11)
plot(reducedDims(seurat.sim)$UMAP, col = getPalette(n_col)[seurat.sim@colData[, cluster_res]],
     pch=16, asp = 1, main = "2D lineages - Clusters Trajectories")
legend("topright", legend=paste0("Cluster - ", levels(seurat.sim@colData[, cluster_res])), pch=16, col=getPalette(n_col))
lines(SlingshotDataSet(seurat.sim), lwd=2, type = 'lineages', col = 'black', show.constraints = const)
dev.off()
message("3. Trajectory lots were done.")

# 10.4. Temporally expressed genes heatmap.
# 10.4.1. Pseudotime is obtained.
for (i in c(1:sum(grepl("slingPseudotime", colnames(seurat.sim@colData))))){
  t <- seurat.sim[[paste0("slingPseudotime_",i)]]

  # 10.4.2. Get the n variable genes.
  Y <- assays(seurat.sim)$logcounts
  var_genes <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:n_var_genes]
  Y <- Y[var_genes,]

  # 10.4.3. Fitting a gam model.
  gam.pval <- apply(Y,1,function(z){
    d <- data.frame(z=z, t=t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
    })
    p <- summary(tmp)[3][[1]][2,3]
    p
  })

  # 10.4.4. Get the top genes selected by p-val and plot the heatmap.
  topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:n_plotted_genes]
  heatdata <- assays(seurat.sim)$logcounts[topgenes, order(t, na.last = NA)]
  heatclus <- seurat.sim@colData[,cluster_res][order(t, na.last = NA)]
  annotation <- data.frame("Cluster" = heatclus, row.names = colnames(heatdata))
  ann_colors <- list("Cluster" = setNames(getPalette(n_col),
                                          levels(heatclus)))
  # 10.4.5. Plot the genes.
  pheatmap(heatdata, cluster_cols = FALSE, 
          color =  colorRampPalette(c("yellow", "red"))(100),
          annotation_col = annotation, annotation_colors = ann_colors,
          show_colnames = FALSE, filename = paste0(dir.name, "/", folders[8], "/temporally_expressed_heatmaps_", selected_res, "_res-trajectory_",i,".pdf"))
  write.xlsx(data.frame("gene" = names(sort(gam.pval, decreasing = FALSE)),
                        "pvalue"=sort(gam.pval, decreasing = FALSE)),
             paste0(dir.name, "/", folders[8], "/temporally_expressed_heatmaps_", selected_res, "_res-trajectory_", i, "_full.xlsx"))
  write.xlsx(data.frame("gene" = names(sort(gam.pval, decreasing = FALSE)[1:n_plotted_genes]),
                        "pvalue"=sort(gam.pval, decreasing = FALSE)[1:n_plotted_genes]),
             paste0(dir.name, "/", folders[8], "/temporally_expressed_heatmaps_", selected_res, "_res-trajectory_", i, "_topN.xlsx"))
}

message("4. Pseudotime heatmap was obtained.")

# Save used objects.
if (graphics) {save(seurat.sim, seurat3D.sim, file = paste0(dir.name, "/", folders[8], "/slingshot_sce_objects.RData"))
} else {save(seurat.sim, file = paste0(dir.name, "/", folders[8], "/slingshot_sce_objects.RData"))}
message("5. R objects were saved.")
