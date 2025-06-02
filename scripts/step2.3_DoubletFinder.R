# https://rpubs.com/kenneditodd/doublet_finder_example
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

message("CONFIGURATION STEP")
# A. Parameters: 
# 1. Load libraries. 
suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("patchwork"))
suppressMessages(library("parallel")) # detectCores()
suppressMessages(library("bitops"))
suppressMessages(library("spam"))
suppressMessages(library("fields"))
require("devtools")
if (!require("DoubletFinder", character.only = TRUE)){
    print("Externally installing...")
    devtools::install_github("maxmeieran/DoubletFinder") # corrected fork version from original remotes::install_github('chris-mcginnis-ucsf/DoubletFinder'), errors in seuratV3
}
suppressMessages(library("DoubletFinder")) # paramSweep_v3()

message("1. Libraries were loaded.")

# 2. Folder configuration. 
input_file = snakemake@input[["seurat_obj"]]
if (is.null(input_file)) {
  cat("Input is empty, skipping processing.\n")
  quit(status = 0)
}
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_scRepertoire", "5_degs", "6_annotation", "7_gs", "8_traj_in", "9_func_analysis", "10_RNAvelocity")
message("2. Folder paths were set.")


# 3. Get variables from Snakemake.  
random_seed = snakemake@params[["random_seed"]]
sample_i = snakemake@params[["sample_i"]]
if (sample_i == "merged"){
    file.create(paste0(dir.name, "/",folders[1], "/seurat_post-qc_singlet.rds"))
    pdf(paste0(dir.name, "/", folders[1], "/5_vlnplot_QC_variables_postDoubletFilter.pdf"))
    dev.off()

    file.create(paste0(dir.name, "/", folders[1], "/6_post_vs_postDoublet_stats.tsv"))
    quit(status = 0)
}
ram = snakemake@resources[["mem_mb"]]
threads = snakemake@threads
message("3. Parameters were loaded.")

# 4. Analysis configuration. 
# RAM configuration.
options(future.globals.maxSize = ram*1024^2)
message(paste0("4. Threads were set at ", threads, "."))
# Set seed.
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}
message(paste0("5. Seed was set at ", random_seed, "."))
message("Configuration finished.")
message("\n")


message("PROCESSING STEP")
# Load seurat object.
seurat <- readRDS(input_file)
assay_type <- seurat@active.assay
message("1. Seurat object was loaded.")

if (length(colnames(seurat)) < 50) {
  message("not enough cells to run a doublet filtering (< 50 cells detected")
  }
if (length(colnames(seurat)) >= 50){ # 50 minimal cells for doublet filtering
	# Pre-process seurat object
	seurat.singlets <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
	seurat.singlets <- FindVariableFeatures(seurat.singlets, selection.method = "vst", nfeatures = 2500)
	seurat.singlets <- ScaleData(seurat.singlets, features = rownames(seurat.singlets))
	seurat.singlets <- RunPCA(seurat.singlets, features = VariableFeatures(object = seurat.singlets), npcs = 50)

	pdf(paste0(dir.name, "/", folders[1], "/7.intermediate_plot-DoubletF.pdf"))
	# Find significant PCs
	stdv <- seurat.singlets[["pca"]]@stdev
	sum.stdv <- sum(seurat.singlets[["pca"]]@stdev)
	percent.stdv <- (stdv / sum.stdv) * 100
	cumulative <- cumsum(percent.stdv)
	co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
	co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
		           percent.stdv[2:length(percent.stdv)]) > 0.1), 
		  decreasing = T)[1] + 1
	min.pc <- min(co1, co2)
	print(paste0("minimal pc - ", min.pc))
	# finish pre-processing
	seurat.singlets <- FindNeighbors(object = seurat.singlets, reduction = "pca", dims = 1:min.pc, future.seed = NULL, k.param = 20, future.seed = NULL)              
	seurat.singlets <- FindClusters(object = seurat.singlets, resolution = 0.4, future.seed = NULL)
	seurat.singlets <- RunUMAP(seurat.singlets, dims = 1:min.pc, n.components = 2, verbose = FALSE, future.seed = NULL)

	message("1.1- Seurat object pre-processed")

	# pK identification (no ground-truth)
	sweep.list <- paramSweep(seurat.singlets, PCs = 1:min.pc, num.cores = min(detectCores(), threads) - 1)
	message("1")
	sweep.stats <- summarizeSweep(sweep.list)
	message("2")
	bcmvn <- find.pK(sweep.stats)
	message("3")
	# Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
	bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
	message("4")
	optimal.pk <- bcmvn.max$pK
	optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
	message("1.2- Optimal pk extracted")

	## Homotypic doublet proportion estimate
	annotations <- seurat.singlets@meta.data$seurat_clusters
	homotypic.prop <- modelHomotypic(annotations) 
	nExp.poi <- round(optimal.pk * nrow(seurat.singlets@meta.data))
	nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

	# run DoubletFinder
	seurat.singlets <- doubletFinder(seu = seurat.singlets, 
		                       PCs = 1:min.pc, 
		                       pK = optimal.pk,
		                       nExp = nExp.poi.adj)
	metadata <- seurat.singlets@meta.data
	dev.off()
	colnames(metadata)[length(metadata)] <- "doublet_finder"
	seurat.singlets@meta.data <- metadata 

	# subset and save
	seurat.singlets <- subset(seurat.singlets, doublet_finder == "Singlet")
	message("2. Doublets were filtered out.")

	# comnpare pre- and post- seurat objects
	print(seurat)
	print(seurat.singlets)
	p1 <- VlnPlot(seurat.singlets, features = c("nFeature_RNA"), pt.size = 0.25, cols = "#9CCCD0", group.by = "orig.ident") + ggtitle("Nº features") + theme(legend.position="bottom") 
	p2 <- VlnPlot(seurat.singlets, features = c("nCount_RNA"), pt.size = 0.25, cols = "#8ADD56", group.by = "orig.ident")  + ggtitle("Nº counts") + theme(legend.position="bottom")
	p3 <- VlnPlot(seurat.singlets, features = c("percent.mt"), pt.size = 0.25, cols = "#F07800", group.by = "orig.ident") + ggtitle("Mitochondrial %") + theme(legend.position="bottom")
	p4 <- VlnPlot(seurat.singlets, features = c("percent.ribo"), pt.size = 0.25, cols = "#E44631", group.by = "orig.ident") + ggtitle("Ribosomal %") + theme(legend.position="bottom")
	p_comp <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4)
	ggsave(paste0(dir.name, "/", folders[1], "/5_vlnplot_QC_variables_postDoubletFilter.pdf"), plot = p_comp, width = 10, height = 8) #, scale = 1.2)

	stats_pre <- c(length(colnames(seurat)), median(seurat@meta.data[["nCount_RNA"]]), median(seurat@meta.data[["nFeature_RNA"]]), median(seurat@meta.data[["percent.mt"]]), median(seurat@meta.data[["percent.ribo"]]))
	stats_post <- c(length(colnames(seurat.singlets)), median(seurat.singlets@meta.data[["nCount_RNA"]]), median(seurat.singlets@meta.data[["nFeature_RNA"]]), median(seurat.singlets@meta.data[["percent.mt"]]), median(seurat.singlets@meta.data[["percent.ribo"]]))
	filtering_df <- data.frame("Number of cells" = c(stats_pre[1],stats_post[1]), "Count median" = c(stats_pre[2],stats_post[2]),"Expressed genes median" = c(stats_pre[3],stats_post[3]), "Mitochondrial percentage median" = c(stats_pre[4],stats_post[4]), "Ribosomal percentage median" = c(stats_pre[5],stats_post[5]))
	row.names(filtering_df) <- c("Post-QC", "Post-DoubletF")
	write.table(filtering_df, file = paste0(dir.name, "/", folders[1], "/6_post_vs_postDoublet_stats.tsv"), sep = "\t", col.names = NA, quote = FALSE)
	message("3. Comparisons were saved.")

	seurat$barcodes <- colnames(seurat)
	singlets.BC <- colnames(seurat.singlets)
	print(length(colnames(seurat)))
	print(length(singlets.BC))
	seurat <- seurat[, singlets.BC]  # filter out singlets from initial raw seurat object
	print(length(colnames(seurat)))
	}
Idents(seurat) <- seurat$orig.ident
saveRDS(seurat, file = paste0(dir.name, "/",folders[1], "/seurat_post-qc_singlet.rds"))
message("4. Seurat object was saved.")

