# edit to cvheck for souporcell results! line 130-

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

message("CONFIGURATION STEP")
# A. Parameters: 
# 1. Load libraries. 
suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
suppressMessages(library("Matrix"))
suppressMessages(library("patchwork"))
message("1. Libraries were loaded.")

# 2. Folder configuration. 
data_dir = snakemake@params[["input_dir"]]
is_cmo = snakemake@params[["is_cmo"]]
dir.name = snakemake@params[["output_dir"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_annotation", "5_scRepertoire", "6_degs", "7_gs", "8_traj_in", "9_func_analysis", "10_RNAvelocity")
message("2. Folder paths were set.")

# 3. Get variables from Snakemake. 
samples_path = snakemake@params[["samples_path"]]
min_cells_per_gene = snakemake@params[["min_cells_per_gene"]]
input_type = snakemake@params[["input_type"]]
technology = snakemake@params[["technology"]]
units_path = snakemake@params[["units_path"]]
sample = snakemake@params[["sample"]]
random_seed = snakemake@params[["random_seed"]]
case =  snakemake@params[["case"]]
ram = snakemake@resources[["mem_mb"]]
message("3. Parameters were loaded.")

# 4. Analysis configuration. 
# RAM configuration.
options(future.globals.maxSize = ram*1024^2)
# Set seed.
if (is.numeric(random_seed)) {
  set.seed(random_seed)
}
message(paste0("4. Seed was set at ", random_seed, "."))

# 5. Set regexp for QC variables calculation.
if (case == "lowercase") {
  mito_grep <- "^mt-"
  ribo_grep <- "^rp[sl][[:digit:]]"
} else if (case == "titlecase") {
  mito_grep <- "^mt-"
  ribo_grep <- "^Rp[sl][[:digit:]]"
} else if (case == "uppercase") {
  mito_grep <- "^MT-"
  ribo_grep <- "^RP[SL][[:digit:]]"
} else {
  message("Please choose a correc case option.")
  quit()
}
message(paste0("5. Case is set as ", case, "."))
message("Configuration finished.")
message("\n")


message("PROCESSING STEP")
units <- read.csv(units_path, header = TRUE, sep = "\t", comment.char = "#")
# 0. Read input and create the expression matrix object.
# If the input file is a fastq file (STARsolo input).
if (input_type == "fastq") {
  if (!(is_cmo)) {
  	file.rename(paste0(data_dir,"/features.tsv"), paste0(data_dir,"/genes.tsv"))
  }
  expression_matrix <- Read10X(data.dir = data_dir)
  if (is.list(expression_matrix)) {
        expression_matrix <- expression_matrix[["Gene Expression"]]
    }
  message("1. Expression matrix was created from alignment files.")
  
# If the input file are matrices (directly read from units.tsv).
} else if (input_type == "matrix") { # units.tsv is loaded
  # If the expression matrix is in 10x like format (matrix, cell barcodes and genes).
  if (technology == "10x") { 
    expression_matrix <- readMM(toString(units[sample,"matrix"]))
    colnames(expression_matrix) <- read.table(toString(units[sample,"cell_names"]))[,1]
    row.names(expression_matrix) <- read.table(toString(units[sample,"gene_names"]))[,1]
    message("1. Expression matrix was created from 10x/CellRanger-like files.")
    
  # If the expression matrix is in TSV format ((genes as row names and cells as column names).
  } else if (technology == "standard") { 
    expression_matrix = read.csv(toString(units[sample,"matrix"]), sep = "\t", header = TRUE, row.names = 1)
    message("1. Expression matrix was created from TSV matrix.")

  } else {
    message("Please specify a correct unit input.")
  }
} else {
  message("Please specify a correct input type.")
}

# 1. Creating a seurat object. 
expression_matrix <- expression_matrix[, colSums(expression_matrix != 0) > 100] # Take into account cells with more than 100 counts, since CreateSeuratObject function breaks.
seurat = CreateSeuratObject(expression_matrix, project = sample, min.features = 1, min.cells = min_cells_per_gene)
message("2. Seurat object was created.")
write.table(colnames(seurat), paste0(dir.name, "/", folders[1], "/barcodes.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# 1.1 Add metadata from samples.tsv file.  
samples_file = read.table(samples_path, sep = "\t", row.names = 1, header = TRUE)
if (is_cmo){ # demultiplexing assigns specific CMO IDs to each demultiplexed cell according to sample
	message("Identifying multiplexed samples.")
	mtpx_file <- file.path(dirname(data_dir), "../multiplexing_analysis/assignment_confidence_table.csv")
  if (!file.exists(mtpx_file)) {
  	stop(paste("Multiplexing annotation file not found at:", mtpx_file))
  }
	assignments <- read.csv(mtpx_file, stringsAsFactors = FALSE)
  rownames(assignments) <- assignments$Barcode
	common_cells <- intersect(colnames(seurat), assignments$Barcode)
  assignments <- assignments[common_cells, ]
  seurat <- AddMetaData(seurat, metadata = assignments)
  
  seurat$CMO.id <- as.character(seurat$Assignment)
	seurat$multiplex_class <- dplyr::case_when(
		is.na(seurat$Assignment)          ~ "No_CMO_Data", # Células del 'raw' sin señal de CMO
		seurat$Assignment == "Multiplet" ~ "Multiplet",
		seurat$Assignment == "Blank"     ~ "Blank",
		seurat$Assignment == "Unassigned"~ "Unassigned",
		TRUE                             ~ "Singlet"
	)
	
	# souporcell results
	# load clusters.tsv
	#multiplexing stats: rescued AND/OR non-rescued methods
	generate_mtpx_stats <- function(seurat_obj, col_name, suffix, output_dir, folders) {
		message(paste0("Generating stats for: ", col_name))
		
		cells_before <- as.data.frame(table(seurat_obj[[col_name]]))
		colnames(cells_before) <- c("Class_ID", "Before_filter")
		
		if (col_name == "CMO.id") {
		  is_singlet <- seurat_obj$multiplex_class == "Singlet"
		} else {
		  invalid_labels <- c("unassigned", "Multiplet", "Blank", "No_Genotype_Data", "Dropped", "Conflict_Discard", "Doublet_Genotype")
		  is_singlet <- !(seurat_obj[[col_name]][,1] %in% invalid_labels)
		}
		
		cells_after <- as.data.frame(table(seurat_obj[[col_name]][is_singlet, ]))
		colnames(cells_after) <- c("Class_ID", "After_filter")
		

		summary_table <- dplyr::left_join(cells_before, cells_after, by = "Class_ID")
		summary_table[is.na(summary_table)] <- 0
		total_row <- data.frame(
		  Class_ID = "TOTAL", 
		  Before_filter = sum(summary_table$Before_filter), 
		  After_filter = sum(summary_table$After_filter)
		)
		summary_table <- rbind(summary_table, total_row)
		
		summary_table$Retention_pct <- round((summary_table$After_filter / summary_table$Before_filter) * 100, 2)
		summary_table[is.na(summary_table)] <- 0 # Por si hay divisiones por cero
		
		file_path <- paste0(output_dir, "/", folders[1], "/0.1_Demultiplexing_stats_", suffix, ".tsv")
		write.table(summary_table, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)
		
		return(summary_table)
	}
	clusters_file <- file.path(dirname(data_dir), "../../../../souporcell/clusters.tsv") 
	if (file.exists(clusters_file)) {
		message("run Souporcell software to rescue Unassigned cells")
		clusters_data <- read.table(clusters_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
		rownames(clusters_data) <- clusters_data$barcode

		common_cells <- intersect(colnames(seurat), clusters_data$barcode)
		message(paste("Common barcodes: ", length(common_cells)))

		# add 'soupAssign' column as metadata
		seurat <- AddMetaData(
				object = seurat, 
				metadata = clusters_data[common_cells, "assignment", drop = FALSE], 
				col.name = "soupAssign"
		)
		seurat$soupAssign <- as.character(seurat$soupAssign)

		### compare and map cellranger w/ soup assignments:
		# Generar una tabla cruzada
		comparison_table <- table(
				CellRanger_CMO = seurat$CMO.id, 
				Soup_Assignment = seurat$soupAssign, 
				useNA = "ifany"
		)
		write.table(comparison_table, 
				        file = paste0(dir.name, "/", folders[1], "/0.3-Comparison_CMO_vs_Soup.tsv"), 
				        sep = "\t", quote = FALSE, row.names=TRUE, col.names = NA)
				        
		# map CMO IDs with soup clustering     
		valid_cmo_rows <- setdiff(unique(seurat$CMO.id), c("Unassigned", "Multiplet", "Blank", "No_CMO_Data", "Dropped", NA))
		unique_soup_vals <- unique(seurat$soupAssign)
		soup_singlets_ids <- unique_soup_vals[!grepl("/", unique_soup_vals) & !is.na(unique_soup_vals)]
		tab <- table(
				seurat$CMO.id[seurat$CMO.id %in% valid_cmo_rows & seurat$soupAssign %in% soup_singlets_ids], 
				seurat$soupAssign[seurat$CMO.id %in% valid_cmo_rows & seurat$soupAssign %in% soup_singlets_ids]
		)
		
		mapping <- apply(tab, 2, function(x) {
		if(sum(x) == 0) return(NA) # Seguridad por si un cluster está vacío en la intersección
			rownames(tab)[which.max(x)]
		})
		seurat$soupAssign_Named <- NA
		can_map <- as.character(seurat$soupAssign) %in% names(mapping)
		seurat$soupAssign_Named <- mapping[as.character(seurat$soupAssign)]
		seurat$soupAssign_Named[can_map] <- mapping[as.character(seurat$soupAssign[can_map])]
				
		# Rescue soup-corrected cells: Cellranger 'unassigned' cells with a propper soup result
		seurat$CMO.id_Rescued <- seurat$CMO.id
		to_rescue <- which(seurat$CMO.id == "Unassigned" & !is.na(seurat$soupAssign_Named))

		seurat$CMO.id_Rescued[to_rescue] <- seurat$soupAssign_Named[to_rescue]

		message(paste("Unassigned cells rescued: ", length(to_rescue)))
		
		if(length(to_rescue) > 0) {
    	seurat$CMO.id_Rescued[to_rescue] <- seurat$soupAssign_Named[to_rescue]
		}
		# Guardamos la tabla de comparación final para verificar el éxito
		final_comparison <- table(Original = seurat$CMO.id, Rescued = seurat$CMO.id_Rescued)
		write.table(final_comparison, 
				        file = paste0(dir.name, "/", folders[1], "/0.2-Final_Rescue_Summary.tsv"), row.names=TRUE, col.names = NA,
				        sep = "\t", quote = FALSE)
			# --- Souporcell correction: Rescued stats ---
		stats_rescued <- generate_mtpx_stats(
				seurat_obj = seurat, 
				col_name = "CMO.id_Rescued", 
				suffix = "rescued_soup", 
				output_dir = dir.name, 
				folders = folders
		)
	}
	
	
	# --- CellRanger stats filtering ---
	stats_original <- generate_mtpx_stats(
		  seurat_obj = seurat, 
		  col_name = "CMO.id", 
		  suffix = "cellranger", 
		  output_dir = dir.name, 
		  folders = folders
	)

	message("Demultiplexing stats (original and rescued) were generated.")
	
	# Final filtering
	if (file.exists(clusters_file)) {
		keep_idx <- !(seurat@meta.data$CMO.id_Rescued %in% c("Unassigned", "Blank", "Multiplet"))
		seurat <- seurat[, keep_idx]
	} else {
		keep_idx <- !(seurat@meta.data$CMO.id %in% c("Unassigned", "Blank", "Multiplet"))
		seurat <- seurat[, keep_idx]
	}

	message("Demultiplexing finished")
} 

if (dim(samples_file)[2] != 0){
	for (i in 5:length(colnames(samples_file))) {
	  seurat <- AddMetaData(seurat, samples_file[sample, i], col.name = colnames(samples_file)[i])
	  }
}

message("3. Metadata from sample.tsv was added to Seurat object.")

# 1.1.1 Add specific cell metadata from metadata.tsv file.
if (input_type == "matrix") {
  if (file.exists(toString(units[sample,"metadata"]))){
    meta_file = read.table(toString(units[sample,"metadata"]), sep = "\t", row.names = 1, header = TRUE)
    meta_file <- subset(meta_file, row.names(meta_file) %in% row.names(seurat@meta.data))
    for (i in 1:length(colnames(meta_file))) {
        seurat <- AddMetaData(seurat, meta_file[,i], col.name = colnames(meta_file)[i])
    }
    message("3.1 Metadata from metadata.tsv was added to Seurat object.")
  } else {
    message("Metadata.tsv was not found.")
  }
}

#1.2 Set idents to avoid new idents based on shared CB names. 
Idents(seurat) <- rep(sample, length(colnames(seurat$RNA@data)))

# 2. Preprocessing: Get QC-related values.
# 2.1. Mitochondrial genes - check levels of expression for mt genes. 
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = mito_grep)
# 2.2. Ribosomal genes - check levels of expression for rb genes.
seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = ribo_grep)
message("4. Ribosomal and mitochondrial percentages were calculated.")

# 2.3. QC: violin plots.
p1 <- VlnPlot(seurat, features = c("nFeature_RNA"), pt.size = 0.25, cols = "#9CCCD0") + ggtitle("Nº features") + theme(legend.position="bottom") 
p2 <- VlnPlot(seurat, features = c("nCount_RNA"), pt.size = 0.25, cols = "#8ADD56")  + ggtitle("Nº counts") + theme(legend.position="bottom")
p3 <- VlnPlot(seurat, features = c("percent.mt"), pt.size = 0.25, cols = "#F07800") + ggtitle("Mitochondrial %") + theme(legend.position="bottom")
p4 <- VlnPlot(seurat, features = c("percent.ribo"), pt.size = 0.25, cols = "#E44631") + ggtitle("Ribosomal %") + theme(legend.position="bottom")
p_comp <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4)
ggsave(paste0(dir.name, "/", folders[1], "/1_vlnplot_QC_variables_prefilt.pdf"), plot = p_comp, scale = 1.2, width = 10, height = 8)
message("5. Combined violin plot was generated.")

# 2.4. QC: GenePlot.
scatter1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.25)+ theme(legend.position="bottom") + labs(title = "Mitochondrial % vs Nº counts", x = "Nº counts", y = "Mitochondrial %")
scatter2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.25) + theme(legend.position="bottom") + labs(title = "Nº features vs Nº counts", x = "Nº counts", y = "Nº features")
p_comb2 <- scatter1 + scatter2 + plot_layout(ncol = 2)
ggsave(paste0(dir.name, "/", folders[1], "/2_geneplot_numi_vs_pctmit_ngene.pdf"), plot = p_comb2, scale = 1.5)
message("6. Scaterplots were generated.")

# 2.5. Save RDS: we can use this object to generate all the rest of the data.
saveRDS(seurat, file = paste0(dir.name, "/" ,folders[1], "/seurat_pre-qc.rds"))
message("7. Seurat object was saved.")
