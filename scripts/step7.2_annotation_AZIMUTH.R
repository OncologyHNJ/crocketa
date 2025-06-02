log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

message("CONFIGURATION STEP")
# A. Parameters: 
# 1. Load libraries. 
suppressMessages(library("ggplot2"))
# if (!require("Seurat", character.only = TRUE)){
#     message("Externally installing...")
#     suppressMessages(devtools::install_github("satijalab/seurat", "seurat5", upgrade_dependencies = FALSE))
#     }
suppressMessages(library("Seurat"))
'
if (!require("SeuratData", character.only = TRUE)){
    message("Externally installing...")
    suppressMessages(devtools::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE))}
    '
suppressMessages(library("SeuratData")) 
'
if (!require("Azimuth", character.only = TRUE)){
    message("Externally installing...")
    suppressMessages(devtools::install_github("satijalab/azimuth", quiet = TRUE))}
    '
suppressMessages(library("Azimuth"))
# suppressMessages(library("patchwork"))
suppressMessages(library("dplyr"))
suppressMessages(library("tibble"))

message("1. Libraries were loaded.")

# 2. Folder configuration. 
cond = snakemake@params[["ident_cond"]]
input_file = snakemake@input[["seurat_obj"]]
dir.name = snakemake@params[["output_dir"]]
# copy folders vble from find-clusters.R
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_scRepertoire", "5_degs", "6_annotation", "7_gs", "8_traj_in", "9_func_analysis", "10_RNAvelocity")
case = snakemake@params[["case"]]
ref_Data = snakemake@params[["reference"]] # reference dataset for Azimuth
message("2. Folder paths were set.")

# 3. Get variables from Snakemake.  
random_seed = snakemake@params[["random_seed"]]
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
# Load seurat object.
seurat <- readRDS(input_file)
if (!grepl("^4.", Version(seurat))){
  seurat <- suppressMessages(UpdateSeuratObject(seurat)) # update if needed
  }
assay_type <- seurat@active.assay
message("1. Seurat object was loaded.")

# define selected cond according to available resolutions
if(grepl("0.", cond)){
    cond <- paste0(assay_type, "_snn_res.", cond)
    if(!(cond %in% colnames(seurat@meta.data))){
      stop("The specified resolution is not available.")
    }
  } else {
    if (!(cond %in% colnames(seurat@meta.data))){
      stop("The specified condition is not available.")
    } 
  }

# Analysis through Azimuth

message("2. starting Azimuth Analysis")

options(timeout = 240) # 60 seconds by default, but not enough
InstallData(ref_Data) # install reference data
seurat <- RunAzimuth(seurat, reference = ref_Data, assay = assay_type)

p1 <- DimPlot(seurat, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) # + NoLegend()
ggsave(paste0(dir.name, "/", folders[6], "/4.1_Azimuth_annotation_perCell.pdf"), plot = p1, scale = 1.5)

azm.df <- seurat@meta.data %>% select(predicted.celltype.l2) %>% rownames_to_column("row")
write.table(azm.df, file =  paste0(dir.name, "/", folders[6], "/4.1_Azimuth_annotation.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

#  Barplots
cols <- c(cond, "seurat_clusters")
for (colName in cols){
  pdf(paste0(dir.name, "/", folders[6], "/4.2_Azimuth_Barplot_annotation_perCluster_", colName, ".pdf"), paper = "a4r", width = 20)
  for (cluster in levels(seurat@meta.data[[colName]])){ # per cluster, estimate population for the different labels identified:
    table_ci <- table(seurat@meta.data %>% filter(get(colName) == cluster) %>% select(predicted.celltype.l2))
    print(barplot(sort(table_ci / sum(table_ci) * 100, decreasing = TRUE), cex.names = 0.3, space = 15, 
                  ylim = c(0,100), col = "#69b3a2", las = 1, ylab = "% annotated cells", xlab = "Cell Label", main = paste0("Cells annotated in cluster ", cluster)))
  }
  dev.off()
}

# absolute barplot for the number of cells observed
pdf(paste0(dir.name, "/", folders[6], "/4.2_Azimuth_Barplot_AbsoluteNumbers.pdf"), paper = "a4r", width = 20)
print(barplot(sort(table_ci / sum(table_ci) * 100, decreasing = TRUE), cex.names = 0.3, space = 15, 
                  ylim = c(0,100), col = "#69b3a2", las = 1, ylab = "% annotated cells", xlab = "Cell Label", main = paste0("Cells annotated in assay")))
dev.off()
message("2.3. Azimuth Results are obtained and plotted")
