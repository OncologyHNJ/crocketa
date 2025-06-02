log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

message("CONFIGURATION STEP")
# A. Parameters: 
# 1. Load libraries. 
suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("data.table"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))
suppressMessages(library("BiocParallel"))
suppressMessages(library("openxlsx"))
suppressMessages(library("future"))
suppressMessages(library("presto"))
message("1. Libraries were loaded.")

# 2. Folder configuration. 
dir.name = snakemake@params[["output_dir"]]
input_data = snakemake@input[["seurat_obj"]]
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_scRepertoire", "5_degs", "6_annotation", "7_gs", "8_traj_in", "9_func_analysis", "10_RNAvelocity")
message("2. Folder paths were set.")

# 3. Get variables from Snakemake.  
selected_cond = snakemake@params[["selected_cond"]]
random_seed = snakemake@params[["random_seed"]]
test = snakemake@params[["test"]]
ranking = snakemake@params[["ranking"]]
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
# Load seurat object and set clustering resolution and assay type.
seurat <- readRDS(input_data)
if (!grepl("^4.", Version(seurat))){
  seurat <- suppressMessages(UpdateSeuratObject(seurat)) # update
}
message(paste0("Actual seurat object version: ", Version(seurat)))
assay_type <- seurat@active.assay
message("1. Seurat object was loaded.")

lapply(selected_cond, function(cond){
  if(grepl("[0|1].", cond)){
    cond <- paste0(assay_type, "_snn_res.", cond)
    if(!(cond %in% colnames(seurat@meta.data))){
      stop("The specified resolution is not available.")
    }
  } else {
    if (!(cond %in% colnames(seurat@meta.data))){
      stop("The specified condition is not available.")
    } 
  }
  
  
  if (length(unique(seurat@meta.data[[cond]])) >1){
    seurat@meta.data[[cond]] <- as.factor(seurat@meta.data[[cond]])
    '
    # Check number of cells for the heatmap to avoid errors.
    if (length(colnames(seurat)) < 10000){
      n_cells = length(colnames(seurat))
    } else {
      n_cells = 10000
    }
    ' # might be required if extremely high number of cells analyzed
    
    #Set styles for xlsx files.
    redStyle <- createStyle(fontColour = "#B60A1C", bgFill = "#FFF06A", textDecoration = c("BOLD"))
    greenStyle <- createStyle(fontColour = "#309143", bgFill = "#FFF06A", textDecoration = c("BOLD"))
    
    # 8. Differentially expressed genes between clusters. 
    seurat@meta.data[[cond]] <- as.factor(seurat@meta.data[[cond]])
    Idents(seurat) <- cond
    
    dir.create(paste0(dir.name, "/", folders[5], "/", cond))
    # If the input seurat object is integrated we only perform the FindConservedMarkers.
    if (seurat@active.assay == "integrated") {
      comparisons <- FetchData(seurat, vars = c("assay_name", cond))
      comparisons$assay_cond<-  apply(comparisons[ ,c("assay_name", cond)] , 1 , paste , collapse = "-" )
      if(length(unique(comparisons$assay_cond))>2){
        for (i in 1:length(unique(Idents(seurat)))){
          # 8.1. Table on differentially expressed genes - using basic filterings.
          clusterX.markers <- FindConservedMarkers(seurat, ident.1 = unique(Idents(seurat))[i], min.pct = 0.25, grouping.var = "assay_name", test.use = test) %>% arrange(desc(avg_log2FC)) %>% rownames_to_column("gene")  #min expressed
          write.table(clusterX.markers, file = paste0(dir.name, "/", folders[5], "/", cond, "/cluster", unique(Idents(seurat))[i],".markers.tsv"), sep = "\t",  row.names=F, quote = FALSE)
        }
      } else {
        message("The specified condition overlaps with the assay level.")
      }
      message("2. DEG analysis was done for integrated assay.")
      
      # If the input seurat object is not integrated.
    } else {
      # 8.1. Table on differentially expressed genes - using basic filterings.
      for (i in 1:length(unique(Idents(seurat)))){
        clusterX.markers <- FindMarkers(seurat, ident.1 = unique(Idents(seurat))[i], min.pct = 0.25, test.use = test) %>% arrange(desc(avg_log2FC)) %>% rownames_to_column("gene")  #min expressed
        write.table(clusterX.markers, file = paste0(dir.name, "/", folders[5], "/", cond, "/cluster", unique(Idents(seurat))[i],".markers.tsv"), sep = "\t",  row.names=F, quote = FALSE)
      }
      message("2. Cluster markers were obtained.")
      
      # 8.2. DE includying all genes - needed for a GSEA analysis. 
      if (ranking){
        for (i in 1:length(unique(Idents(seurat)))){
          clusterX.markers <- FindMarkers(seurat, ident.1 = unique(Idents(seurat))[i], min.pct = 0, logfc.threshold = 0, test.use = test) %>% arrange(desc(avg_log2FC))  #min expressed
          wb <- createWorkbook()
          addWorksheet(wb, "DE analysis")
          writeData(wb, "DE analysis", clusterX.markers, rowNames = TRUE)
          conditionalFormatting(wb, "DE analysis", cols = 1:(ncol(clusterX.markers)+1),
                                rows = 2:(nrow(clusterX.markers) + 1), rule = "AND($C2<0, $F2<0.05)",
                                style = greenStyle)
          conditionalFormatting(wb, "DE analysis", cols = 1:(ncol(clusterX.markers)+1),
                                rows = 2:(nrow(clusterX.markers) + 1), rule = "AND($C2>0, $F2<0.05)",
                                style = redStyle)
          legend <- createComment(comment = c("Red means a positive LogFold\n\n", "Green means a negative LogFold"), style = c(redStyle, greenStyle))
          writeComment(wb, "DE analysis", col = 8, row = 2, comment = legend)
          saveWorkbook(wb, paste0(dir.name, "/", folders[5], "/", cond, "/cluster", unique(Idents(seurat))[i],".DE.xlsx"), overwrite = TRUE)    
          
          # 8.2.1. Create RNK file. 
          rnk = NULL
          rnk = as.matrix(clusterX.markers[,2])
          rownames(rnk)= row.names(clusterX.markers)
          write.table(rnk, file = paste0(dir.name, "/", folders[5], "/", cond, "/cluster", unique(Idents(seurat))[i],".rnk"), sep = "\t", col.names = FALSE, quote = FALSE)
        }
        message("3. All genes DEG and RNK file were finished.")
      }
      
      # 8.3. Find TOP markers.
      seurat.markers <- FindAllMarkers(object = seurat, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25, test.use = test)
      groupedby.clusters.markers = seurat.markers %>% group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% slice_head(n=5)
      message("top genes for heatmap - 5 per cluster -")

      # 8.3.1. HeatMap top10.f
      # setting heatmap text label sizes
      ytext_value <- 70/(length(unique(seurat@meta.data[[cond]]))*1.30)
      if (ytext_value > 10) {ytext_value = 10}
      xtext_value <- 50/(length(unique(seurat@meta.data[[cond]]))*1.25)
      if (xtext_value > 4) {xtext_value = 4}
      
      # setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
      table(seurat@meta.data[[cond]])
      p1 <- DoHeatmap(object = seurat, features = unique(groupedby.clusters.markers$gene), group.by = cond, size = xtext_value, angle = 45, 
                      group.bar = TRUE, draw.lines = F, raster = FALSE) +
        scale_fill_gradientn(colors = c("blue", "white", "red")) + guides(color=FALSE) + 
        theme(axis.text.y = element_text(size = ytext_value)) + 
        theme(legend.position="bottom") 
      ggsave(paste0(dir.name, "/", folders[5], "/", cond, "/1.1_heatmap_top5markers.pdf"), plot = p1, scale = 2)
      write.table(groupedby.clusters.markers, file = paste0(dir.name, "/", folders[5], "/", cond, "/1.1_heatmap_top5markers.tsv"), sep = "\t",  row.names=F, quote = FALSE)
      message("marker heatmap saved")
      pdf(paste0(dir.name, "/", folders[5], "/", cond, "/1.2_FeaturePlot_top5markers.pdf"))
      print(DotPlot(seurat, features = unique(groupedby.clusters.markers$gene), group.by = cond,
      	    cols = "RdBu") + RotatedAxis() +  
      ggplot2::theme(axis.text.x=element_text(angle=90,hjust=1)) + 
        theme(axis.text.y = element_text(size = ytext_value), axis.text.x = element_text(size = xtext_value)))
      stacked <- ifelse(length(groupedby.clusters.markers$gene) ==1, FALSE, TRUE) # if only one gene, non-stacked vlnplot
      print(VlnPlot(seurat, groupedby.clusters.markers$gene, stack = stacked, flip = T) + NoLegend() + 
        theme(axis.text.y = element_text(size = ytext_value), axis.text.x = element_text(size = xtext_value)))
      for (i in seq(1, length(groupedby.clusters.markers$gene), 1)){
    	  print(FeaturePlot(seurat, groupedby.clusters.markers$gene[i], 
    	                  reduction = "umap", label = TRUE, repel = TRUE) &   
    	        scale_colour_gradientn(colours = RColorBrewer::brewer.pal(n = 5, name = "PuRd")))
    	  print(FeaturePlot(seurat, groupedby.clusters.markers$gene[i], 
    	                  reduction = "umap", label = TRUE, repel = TRUE)) # or default color palette
      }
      dev.off()

      message("4. Top marker heatmap was done.")
    }
    
    # 8.4. Save RDS
    saveRDS(seurat, file = paste0(dir.name, "/",folders[5], "/seurat_degs.rds"))
    message("5. Seurat object was saved.")
    
  } else if (length(unique(seurat@meta.data[[cond]])) <=1) {
    if(assay_type == "integrated"){
    	message("The specified condition has only one level.")  
      	if(!file.exists( paste0(dir.name, "/", folders[5], "/seurat_degs.rds"))){
	      saveRDS(seurat, file = paste0(dir.name, "/", folders[5], "/seurat_degs.rds"))
      	}
    } else {
      message("The specified condition has only one level. Nothing to be done.")
    	if(!file.exists( paste0(dir.name, "/", folders[5], "/seurat_degs.rds"))){
      		saveRDS(seurat, file = paste0(dir.name, "/", folders[5], "/seurat_degs.rds"))
    	}	
    }
  }
})
