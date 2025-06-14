log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

message("CONFIGURATION STEP")
# A. Parameters: 
# 1. Load libraries.
suppressMessages(library("Seurat"))
suppressMessages(library("velocyto.R"))
suppressMessages(library("SeuratWrappers"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("future"))
message("1. Libraries were loaded.")

# 2. Folder configuration. 
input_data = snakemake@input[["seurat_obj"]]
velocyto_dir = snakemake@params[["velocyto_dir"]]
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_scRepertoire", "5_degs", "6_annotation", "7_gs", "8_traj_in", "9_func_analysis", "10_RNAvelocity")
message("2. Folder paths were set.")

# 3. Get variables from Snakemake.
selected_res = snakemake@params[["selected_res"]]
random_seed = snakemake@params[["random_seed"]]
downsampling = snakemake@params[["downsampling"]]
n_cells = snakemake@params[["n_cells"]]
ram = snakemake@resources[["mem_mb"]]
threads = snakemake@threads
message("3. Parameters were loaded.")

# 4. Analysis configuration. 
# RAM configuration.
options(future.globals.maxSize = ram*1024^2)
#Set parallelization.
plan("multisession", workers = threads)
message(paste0("4. Threads were set at ", threads, "."))
# Set seed.
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}
message(paste0("5. Seed was set at ", random_seed, "."))
message("Configuration finished.")
message("\n")


message("PROCESSING STEP")
# Load seurat object
seurat = readRDS(input_data)

# 12. Compute velocity and velocity plots.
# 12.1. Check cluster resoltution.
assay_type <- seurat@active.assay
cluster_res <- paste0(assay_type, "_snn_res.", selected_res)
if (!(cluster_res %in% colnames(seurat@meta.data))){
  stop("Specified resolution is not available.")
}
message("1. Seurat object was loaded.")

# If Seurat objects are not integrated we need to add the velocyto matrices.
if (!(seurat@active.assay == "integrated" || seurat@project.name == "merged")) {
  # 12.2. Get velocity matrices and place them in a list.
  velo_names = c("spliced", "unspliced", "ambiguous")
  vel_matrices = list()
  for (name in velo_names) {
    vel_matrices[[name]] <- Read10X(data.dir = paste0(velocyto_dir, name))
  }

  # 12.3. Load Velocyto matrices as seurat assays.
  for (name in velo_names) {
    vel_matrices[[name]] <- vel_matrices[[name]][, which(x = colnames(vel_matrices[[name]]) %in% colnames(seurat))] 
    seurat[[name]] <- CreateAssayObject(counts = vel_matrices[[name]])
  }
}
message("2. Seurat object was prepared to perform RNA velocity.")

# 12.4. Downsampling (optional).
if (downsampling == TRUE && n_cells < length(rownames(seurat@meta.data))){
  #n_total_cells = length(rownames(seurat@meta.data))
  random_sample = sample(x = rownames(seurat@meta.data), size = n_cells, replace = FALSE)
  seurat <- subset(x = seurat, cells = random_sample) 
  message("2.5. Downsampling was performed.")
}

# 12.5. Set specific cluster labels as idents.
Idents(seurat) <- seurat@meta.data[[cluster_res]]

# 12.6. Run velocyto from the wrapper.
seurat <- RunVelocity(object = seurat, deltaT = 1, kCells = 25, fit.quantile = 0.02)
message("3. RNA velocity was calculated.")

# 12.7. Obtain palette.
n_col <- length(levels(seurat))
getPalette <- colorRampPalette(brewer.pal(9,'Set1'))

# 12.8. Set colours to clusters. 
ident.colors <- getPalette(n_col) 
names(ident.colors) <- levels(seurat)
cell.colors <- ident.colors[Idents(seurat)]
names(cell.colors) <- colnames(seurat)

# 12.9. Create the RNA velocity plot.
pdf(paste0(dir.name, "/",folders[10], "/RNA_velocity_plot.pdf"))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
show.velocity.on.embedding.cor(emb = Embeddings(object = seurat, reduction = "umap"), vel = Tool(object = seurat, 
                               slot = "RunVelocity"), n = 100, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE,  cell.border.alpha = 0.1, xlab = "UMAP1", ylab = "UMAP2", 
                               main = paste0("RNA velocity plot - Cluster resolution: ", selected_res))
opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
on.exit(par(opar))
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend("topright", inset = c(0.05, 0.115), legend=paste0("Cluster - ", levels(seurat)),
       pch=16, col=getPalette(n_col))
dev.off()
message("4. Velocity plot was produced.")


#12.10. This plot the sample coloring the cell by its origin assay (integrated objects).
if (seurat@active.assay == "integrated") {
  Idents(seurat) <- seurat@meta.data[["assay_name"]]
  ident.colors <- getPalette(n_col)
  names(ident.colors) <- levels(seurat)
  cell.colors <- ident.colors[Idents(seurat)]
  names(cell.colors) <- colnames(seurat)

  # 12.11. Create the RNA velocity plot.
  pdf(paste0(dir.name, "/",folders[10], "/RNA_velocity_plot_by_assay.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  show.velocity.on.embedding.cor(emb = Embeddings(object = seurat, reduction = "umap"), vel = Tool(object = seurat,
                                 slot = "RunVelocity"), n = 100, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5),
                                 cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1,
                                 do.par = FALSE,  cell.border.alpha = 0.1, xlab = "UMAP1", ylab = "UMAP2",
                                 main = paste0("RNA velocity plot - Origin assay"))
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend("topright", inset = c(0.05, 0.115), legend=paste0("Assay - ", levels(seurat)),
         pch=16, col=getPalette(n_col))
  dev.off()
  message("4.5. Velocity plot from integrated object was produced.")

}

#12.12. This plot the sample coloring the cell by its origin assay (mergedd objects).
if (seurat@project.name == "merged") {
  Idents(seurat) <- "assay_name"
  ident.colors <- getPalette(n_col)
  names(ident.colors) <- levels(seurat)
  cell.colors <- ident.colors[Idents(seurat)]
  names(cell.colors) <- colnames(seurat)

  # 12.13. Create the RNA velocity plot.
  pdf(paste0(dir.name, "/",folders[10], "/RNA_velocity_plot_by_assay.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  show.velocity.on.embedding.cor(emb = Embeddings(object = seurat, reduction = "umap"), vel = Tool(object = seurat,
                                 slot = "RunVelocity"), n = 100, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5),
                                 cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1,
                                 do.par = FALSE,  cell.border.alpha = 0.1, xlab = "UMAP1", ylab = "UMAP2",
                                 main = paste0("RNA velocity plot - Origin assay"))
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend("topright", inset = c(0.05, 0.115), legend=paste0("Assay - ", levels(seurat)),
         pch=16, col=getPalette(n_col))
  dev.off()
  message("4.5. Velocity plot from integrated object was produced.")
}

# 12.14. Save seurat object with RNA velocity slots. 
saveRDS(seurat, file = paste0(dir.name, "/",folders[10], "/seurat_velocity.rds"))
message("5. Seurat object was saved.")
