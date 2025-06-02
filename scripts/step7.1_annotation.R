log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

message("CONFIGURATION STEP")
# A. Parameters: 
# 1. Load libraries. 
suppressMessages(library("Seurat"))
# suppressMessages(library("SeuratDisk"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))
# suppressMessages(library("writexl"))
suppressMessages(library("patchwork"))
suppressMessages(library("ggraph"))
suppressMessages(library("igraph"))
suppressMessages(library("future"))
## packages for scType
suppressMessages(library("HGNChelper"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("data.tree"))
suppressMessages(library("openxlsx"))
suppressMessages(library("RColorBrewer"))
## packages for SingleR
suppressMessages(library("SingleR"))
suppressMessages(library("scRNAseq"))
suppressMessages(library("scuttle"))
suppressMessages(library("scater"))

# load sctype scripts
sctype_gsPrepare = snakemake@params[["sctype_gsPrepare"]]
sctype_score = snakemake@params[["sctype_score"]]
source(sctype_gsPrepare) # "scripts/scType-gene_sets_prepare.R"
source(sctype_score) # "scripts/scType-sctype_score.R"

'
# load gene set preparation function from github
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function from github
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
'
message("1. Libraries were loaded.")
# 2. Folder configuration. 
case = snakemake@params[["case"]]
cond = snakemake@params[["ident_cond"]]
input_file = snakemake@input[["seurat_obj"]]
dir.name = snakemake@params[["output_dir"]]
# copy folders vble from find-clusters.R
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_scRepertoire", "5_degs", "6_annotation", "7_gs", "8_traj_in", "9_func_analysis", "10_RNAvelocity")
interest_feat = snakemake@params[["interest_feat"]] # IF TRUE, check predefined genes, else define genes of interest
gmx_dir = snakemake@params[["gmx_dir"]]
tissue_scType = snakemake@params[["tissue_scType"]] # NULL or e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
restrict_REF_scType = snakemake@params[["restrict_REF_scType"]]
ref_function = snakemake@params[["ref_singleR"]] # reference function for sceM in singleR
label_singleR = snakemake@params[["label_singleR"]] # label where cell type is specified in each of the reference datasets
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
message("1. Seurat object was loaded.")
assay_type <- seurat@active.assay
if (grepl("^3.", Version(seurat))){
	message("Update Seurat Object")
	seurat <- suppressMessages(UpdateSeuratObject(seurat)) # update if needed
 }

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
Idents(seurat) <- "seurat_clusters"
# Manual data explore: markers

if (case == "titlecase"){
  interest_feat <- str_to_title(interest_feat)
}
if (case == "lowercase"){
  interest_feat <- str_to_lower(interest_feat)
}

interest_feat <- interest_feat[interest_feat %in% rownames(seurat)]
if (length(interest_feat) > 1){
	pdf(paste0(dir.name, "/", folders[6], "/1.1_Expression_check_interest.pdf"))
	DotPlot(seurat, features = interest_feat,        
		    cols = "RdBu") + RotatedAxis() +  
	  ggplot2::theme(axis.text.x=element_text(angle=90,hjust=1))
	stacked <- ifelse(length(interest_feat) ==1, FALSE, TRUE) # if only one gene, non-stacked vlnplot
	VlnPlot(seurat, interest_feat, stack = stacked, flip = T)
	for (i in seq(1, length(interest_feat), 1)){
		print(FeaturePlot(seurat, interest_feat[i], 
		                  reduction = "umap", label = FALSE, repel = TRUE, order = TRUE) &   
		        scale_colour_gradientn(colours = RColorBrewer::brewer.pal(n = 5, name = "PuRd")))
		        }
	dev.off()
}


# gmx_dir is the directory were all gmx of interest are stored
# if specified, then perform graphics

if (is.null(gmx_dir) == F){
    sign_ext <- ".*gmx$"
    file_Names <- list.files(path = gmx_dir, pattern = sign_ext, recursive = T, full.names = T)
    files <- lapply(file_Names, read.table, 
                    header = TRUE, sep = "\t", quote = NULL)
    names(files) <- list.files(path = gmx_dir, pattern = sign_ext,
                               recursive = T, full.names = F)
    ## 1- Stacked VlnPlots for each list of genes (each signature file)
    # parsing below created for genes ranked according to log2FC
    # if only gene names provided, just skip the parsing below and pull gene names
    for (i in seq(1,length(file_Names))) {
      pdf(paste0(dir.name, "/", folders[6], "/1.2_", names(files)[i], "_Marker_gene_expression.pdf"), height = 30)  
      print(paste0("Testing gene expression for: ", names(files)[i]))
      df_perFile <- files[[i]]
      for (j in seq(1,length(colnames(df_perFile)))) {
        print(paste0("Testing gene expression for: ", names(files)[i], "--", colnames(df_perFile)[j]))
        gene_sign <- df_perFile[,j]
        if (case == "titlecase"){
	  gene_sign <- str_to_title(gene_sign)
	}
	if (case == "lowercase"){
	  gene_sign <- str_to_lower(gene_sign)
	}
        gene_sign <- gene_sign[gene_sign %in% rownames(seurat)[rowSums(seurat) > 0]]
        if (length(gene_sign) > 0){
          if (length(gene_sign) == 1){
            print(VlnPlot(seurat, gene_sign, stack = F, flip = T) +
                  ggtitle(colnames(df_perFile)[j]) + NoLegend())
          } else{
            seurat <- AddModuleScore(seurat, 
                                      list(gene_sign), 
                                      ctrl = length(gene_sign) * 0.3, 
                                      name = paste0("AddModuleScore_", colnames(df_perFile)[j]), seed = random_seed)
            print(VlnPlot(seurat, gene_sign, stack = T, flip = T) +
              ggtitle(colnames(df_perFile)[j]) + NoLegend())
          }
          if (length(gene_sign) < 5){
            for (i in seq(1, length(gene_sign),1))
            print(FeaturePlot(seurat, gene_sign[i], 
                       reduction = "umap", label = FALSE, repel = TRUE, order = TRUE) &   
             scale_colour_gradientn(colours = RColorBrewer::brewer.pal(n = 5, name = "PuRd")))
    
          print(DoHeatmap(seurat, gene_sign, angle = 90) &
            scale_fill_gradientn(colours =rev(RColorBrewer::brewer.pal(n = 8, 
                                                                name = "BrBG"))))        
          }
        }
      }
      dev.off()
    
   message("check addModuleScores")
    pdf(paste0(dir.name, "/", folders[6], "/1.3_", names(files)[i], "_AddModuleScore.pdf"))
    for (module in colnames(seurat@meta.data)[grepl("AddModuleScore_",
                                          colnames(seurat@meta.data))]){
      if (all(seurat[[module]] != "NaN")){
          print(FeaturePlot(seurat,
                        features = module, 
                        label = FALSE, 
                        repel = TRUE, order = TRUE) + theme(plot.title = element_text(size = 8)) +
                scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, 
                                                                name = "RdBu"))))
      }                                      
    }
    dev.off()
    
    seurat@meta.data[,grepl("AddModuleScore_", colnames(seurat@meta.data))] <- NULL # remove the scoring columns

}
}

message("features of interest saved")

# Analysis: A) scType, B) SingleR

## A) scType
message("2. starting scType Analysis")

# DB file
# db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx" # load original sctype DB
db_ = snakemake@params[["scType_ref"]] # load edited sctype ddbb
# if tissue_scType not known (NULL) then it is inferred computationally:
if (is.null(tissue_scType)){
  sctype_autoDetect = snakemake@params[["sctype_autoDetect"]]
  source(sctype_autoDetect) # "./scripts/scType-auto_detect_tissue_type.R"f
  pdf(paste0(dir.name, "/", folders[6], "/2.0_Inferred_tissue_scType.pdf"), paper = "a4r")
  tissue_guess = auto_detect_tissue_type(path_to_db_file = db_, seuratObject = seurat, scaled = TRUE, assay = assay_type)  
  tissue_scType = pull(tissue_guess, "tissue")[1]
  dev.off()
}

# prepare gene sets
gs_list = gene_sets_prepare(normalizePath(db_), tissue_scType)

# if restrict to any cell pattern:
if (is.null(restrict_REF_scType) == FALSE){ # if restrict_REF is not null, we want to specify restricted labels of my reference data
	gs_list_restrict <- gs_list
	gs_list_restrict$gs_positive <- gs_list_restrict$gs_positive[grepl(restrict_REF_scType, names(gs_list_restrict$gs_positive))]
	gs_list_restrict$gs_negative <- gs_list_restrict$gs_negative[grepl(restrict_REF_scType, names(gs_list_restrict$gs_negative))]
} else {gs_list_restrict <- NULL}

for (x in c("non_restrict", "restrict")){ # analyze both options
    if (x == "restrict"){
        gs_i <- gs_list_restrict
    } else {gs_i <- gs_list}
	if(is.null(gs_i) == FALSE){ # if not null, then restriction is applied and both are analysed (otherwise only gs_list is checked)
		# get cell-type by cell matrix
		es.max = sctype_score(scRNAseqData = seurat[[assay_type]]@scale.data, scaled = TRUE, 
							gs = gs_i$gs_positive, gs2 = gs_i$gs_negative) 
		# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
		# In case Seurat is used, it is either seurat[[assay_type]]@scale.data (default), seurat[["SCT"]]@scale.data, in case sctransform is used for normalization,
		# or seurat[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

		# merge by cluster
		cL_results = do.call("rbind", lapply(unique(seurat@meta.data[["seurat_clusters"]]), function(cl){
		es.max.cl = sort(rowSums(es.max[ ,rownames(seurat@meta.data[seurat@meta.data[["seurat_clusters"]]==cl, ])]), decreasing = !0)
		head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat@meta.data[["seurat_clusters"]]==cl)), 10)
		}))

		suppressMessages({
		# create barplot per cluster to show results
		dir.create(paste0(dir.name, "/", folders[6], "/2.1_scType_annotation_perCluster/"))
		pdf(paste0(dir.name, "/", folders[6], "/2.1_scType_Barplot_annotation_perCluster_", x,".pdf"), paper = "a4r", width = 20)
		for (pop in levels(cL_results$cluster)){ # per cluster, estimate population for the different labels identified:
		table_ci <- cL_results %>% filter(cluster == pop) %>% select(c(scores))
		print(barplot(table_ci$scores, cex.names = 0.3, space = 15,
						col = "#69b3a2", las = 1, names.arg = rownames(table_ci), ylab = "scType Score", xlab = "Cell Label", main = paste0("Cells annotated in cluster ", pop)))
		write.xlsx(table_ci, paste0(dir.name, "/", folders[6], "/2.1_scType_annotation_perCluster/2.1_scType_annotation_Cluster_",pop,"_",x,".xlsx"), row.names = TRUE)
		}
		dev.off()
		})

		message("Annotation performed")
		## 2- FeaturePlots for each list of genes (each signature assessed, positive & negative)
		## ADD MODULE SCORES FOR EACH SIGNATURE
		if (length(colnames(seurat)) < 100000){
			for (i in c(1:length(gs_i))){ # pos & neg
			    for (j in c(1:length(gs_i$gs_positive))){
	    			gene_sign <- unlist(gs_i[[i]][j])
	    			if (case == "titlecase"){
	    			gene_sign <- str_to_title(gene_sign)
	    			}
	    			if (case == "lowercase"){
	    			gene_sign <- str_to_lower(gene_sign)
	    			}
	    			if (sum(gene_sign %in% rownames(seurat)) > 1){ 
	    			seurat <- AddModuleScore(seurat, 
	    			list(gene_sign), ctrl = length(gene_sign) * 0.3, 
	    			name = paste0(names(gs_i)[i], "_",x,"_", names(gs_i[[i]][j])), 
	    			seed = random_seed)
	    			
	    			}
			    }
			}
		
			pdf(paste0(dir.name, "/", folders[6], "/2.2_scType_AddModuleScore.pdf"), width = 8)
			for (module in 
				colnames(seurat@meta.data)[grepl("gs_", # pattern depends on sign_ext
														colnames(seurat@meta.data))]){
			if (length(unique(seurat@meta.data[[module]])) > 1){
				print(FeaturePlot(seurat,
								features = module,
								label = FALSE, 
								repel = TRUE, order = TRUE) +
						scale_colour_gradientn(colours = rev(brewer.pal(n = 11, 
																		name = "RdBu"))))
			}
			}
			dev.off()
				
			# remove these extra columns from metadata  once they are plotted
			seurat@meta.data[,grepl("gs_", # pattern depends on sign_ext
														colnames(seurat@meta.data))] <- NULL
		}
		## and select top cluster
		sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

		# set low-confident (low ScType score) clusters to "unclassifed" (< 20% cells)
		sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/5] = "Unclassified"
		# print(sctype_scores[,1:3])


		# plot cell types on UMAP
		seurat@meta.data$scType_annotation = ""
		for(j in unique(sctype_scores$cluster)){
		cl_type = sctype_scores[sctype_scores$cluster==j,]; 
		seurat@meta.data$scType_annotation[seurat@meta.data[["seurat_clusters"]] == j] = as.character(cl_type$type[1])
		}

		p1 <- DimPlot(seurat, reduction = "umap", label = FALSE, repel = TRUE, group.by = 'scType_annotation') +
				theme(legend.text=element_text(size=rel(0.5)))
		ggsave(paste0(dir.name, "/", folders[6], "/2.3_scType_annotation_",x,".pdf"), plot = p1, scale = 1.5)
		write.xlsx(cL_results, paste0(dir.name, "/", folders[6], "/2.3_scType_annotation_scores_perCluster_",x,".xlsx"), row.names = TRUE)

		message("2.3. scType Results are obtained and plotted")
        
		# prepare edges
		cL_results=cL_results[order(cL_results$cluster),]; edges = cL_results; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

		# prepare nodes
		nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
		ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3", "#00FF66", "#666666", "#FF3399", "#660033", "#6600FF", "#0099FF", "#003333", "#0033FF", "#666600", "#990099", "#CC99FF", "#CCFF99")
		for (i in 1:length(unique(cL_results$cluster))){
			dt_tmp = cL_results[cL_results$cluster == unique(cL_results$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
		}
		nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
		files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
		nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
		if (!(length(unique(nodes[["cluster"]]))>length(ccolss))){
			message("plot graph nodedge")
			mygraph <- graph_from_data_frame(edges, vertices=nodes)
	
			# Make the graph
			gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
            geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + 
            geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
            theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5))) + 
            geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
		    pdf(paste0(dir.name, "/", folders[6], "/2.4_scType_annotation_summary_",x,".pdf"), width = 15)
		    p2 <- gridExtra::grid.arrange(DimPlot(seurat, group.by = "seurat_clusters", reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss), gggr, ncol = 2)
		    p2
		    dev.off()
		    message("2.4. scType Summary of the scType analysis is plotted")
		}
		cond <- "scType_annotation"
		data <- as.data.frame(table(seurat$orig.ident, seurat@meta.data[[cond]]))
		data <- transform(data, rel = Freq / ave(Freq, Var1, FUN = sum))
		pdf(paste0(dir.name, "/", folders[6], "/2.5_scType-perCellType_perOrig_plot.pdf"), width = 10)
		# A) individual barplots per cell type
		print(ggplot(data, aes(x = factor(Var1), y = Freq, fill = factor(Var2))) +
		    geom_bar(stat = "identity") +
		  labs(x = "Condition", y = "Number of Cells", fill = "Cell type") +
		  theme_minimal() +
		  theme(legend.position = "top") +
		  facet_wrap(~Var2, scales = "free_y"))
		# B) stacked varplot: absolute
		print(ggplot(data, aes(x = factor(Var1), y = Freq, fill = factor(Var2))) + 
		    geom_bar(stat = "identity") +
		    geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size = 3) +
		    labs(x = "Condition", y = "Number of Cells", fill = "Cell type") +
		    theme_minimal() +
		    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

		# C) stacked varplot: relative
		print(ggplot(data, aes(x = factor(Var1), y = rel, fill = factor(Var2))) +
		    geom_bar(stat = "identity") +
		    geom_text(aes(label = round(rel,3)), position = position_stack(vjust = 0.5), size = 3) +
		    labs(x = "Condition", y = "% of Cells", fill = "Cell type") +
		    theme_minimal() +
		    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
		dev.off()
		
		message("scType analysis finished")
	}
}


## B) singleR
if (!is.null(ref_function)) {
  if (is.null(label_singleR)){
    stop("Provide a value for label_singleR variable")
  }
  else {
    message("3. starting singleR Analysis")
    for (colName in c(cond, "seurat_clusters")){
		for (i in c(1:length(ref_function))){
		  ref_DS <- gsub("\\()", "", ref_function[i])
		  message(paste0("3.0. Analyzing ", ref_DS))
		  myref <- function() {
		    eval(parse(text = ref_function[i]))
		  }
		  sceM <- myref() # load reference data
		  ref_DS <- paste0(ref_DS, "_", colName)
		  dir.create(paste0(dir.name, "/", folders[6], "/3-", ref_DS))
		  message("3.1. singleR reference is loaded")
		  seurat@meta.data[,grepl(ref_DS, colnames(seurat@meta.data))] <- NULL # remove previous results
		  if(any(!(label_singleR[i] %in% colnames(colData(sceM))),
		          length(label_singleR) != length(ref_function))){
		    stop("Provide the propper labels for singleR reference datasets")
		  }  
		  sceG <- as.SingleCellExperiment(seurat,
		  colData = DataFrame(seurat@meta.data[, colName])) # load experimental data
		  message("3.1.2 singleR expression data is transformed")

		  # IMPORTANT: Check filtering performed in reference paper to assure data quality is optimal
		    
		  # First of all, remove NA values for non-annotated cells
		  sceM <- sceM[!is.na(sceM[[label_singleR[i]]])]
		  # SingleR() expects reference datasets to be normalized and log-transformed.
		  if ("counts" %in%  names(assays(sceM))){
		  sceM <- logNormCounts(sceM)
		  }
		  # and run SingleR
		  pred.grun <- SingleR(test=sceG, ref=sceM, labels=sceM[[label_singleR[i]]], de.method="wilcox", de.n = 15)
		  # saveRDS(pred.grun, paste0(dir.name, "/",folders[6], "/3.0_SingleR_predict_", ref_DS,"_", colName, ".rds"))
		  message("3.2. singleR result is saved")
		  # table(pred.grun$labels)
		  ## Annotation diagnostics
		  pdf(paste0(dir.name, "/", folders[6], "/3-", ref_DS, "/3.1_SingleR_Annotation_diagnostics_",ref_DS, ".pdf"), width = 15)

		  # Heatmap showing scores for all cells at all reference labels (to validate results)
		  print(plotScoreHeatmap(pred.grun))

		  # Difference between the score for the assigned label and rthe 
		  # median across all labels for each cells . 
		  # Low delta = uncertain assignment
		  print(plotDeltaDistribution(pred.grun, ncol = 3)) #, labels.use = c("Astro", "Ependy", "Epith1", "Epith2") # if I wish to remove some cell labels 

		  # pruneScores() function remove potentially poor quality assignments
		  # outlier-based approach 
		  # summary(is.na(pred.grun$pruned.labels))
		  # table(Label=pred.grun$labels,
		  #       Lost=is.na(pred.grun$pruned.labels))
		  dev.off()
		  '
		  ## not very meaningful and really long step creating several heatmaps:
		  # Examine expression of marker genes for each label in
		  # test dataset. We first: extract markers identity
		  all.markers <- metadata(pred.grun)$de.genes
		  sceG$labels <- pred.grun$labels
		  # and second use them in the heatmap. If a single cell label is robust,
		  # we would hope it to be strong expressed in that label`s markers
		  pdf(paste0(dir.name, "SingleR_auto/marker_heatmap_", colName, ".pdf"))
		  for (i in seq(1,length(names(all.markers)))){
		    plotHeatmap(sceG, order_columns_by="labels",
		                features=unique(unlist(all.markers[i])), main = names(all.markers)[i])
		  }
		  dev.off()
		  '
		  ### AND ADD ANNOTATION AS METADATA OF SEURAT
		  # pred.grun@rownames
		  seurat@meta.data <- seurat@meta.data %>% rownames_to_column("row") %>% rowwise()  %>%
		    mutate("SingleR" = pred.grun[row, "pruned.labels"]) %>%
		    rename_with(~ref_DS, "SingleR") %>%
		    as.data.frame() %>% column_to_rownames("row")
		  p3 <- DimPlot(seurat, group.by = ref_DS, label = T, reduction='umap') +
		          theme(legend.text=element_text(size=rel(0.5)))
		  ggsave(paste0(dir.name, "/", folders[6], "/3-", ref_DS, "/3.2_SingleR_annotation_",ref_DS, ".pdf"), plot = p3, width = 7.5, scale = 1.5)

		  ## CREATE BARPLOT TO QUANTIFY ABUNDANCE OF CELL TYPES PER CLUSTER
		  # useful to define very clearly annotated clusters
		  pdf(paste0(dir.name, "/", folders[6],  "/3-", ref_DS, "/3.3_SingleR_Barplot_annotation_",ref_DS, ".pdf"), paper = "a4r", width = 20)
		  for (cluster in levels(seurat@meta.data[[colName]])){ # per cluster, estimate population for the different labels identified:
		    table_ci <- table(seurat@meta.data %>% filter(get(colName) == cluster) %>% select(ref_DS))
		    print(barplot(sort(table_ci / sum(table_ci) * 100, decreasing = TRUE), cex.names = 0.3, space = 15, 
		    ylim = c(0,100), col = "#69b3a2", las = 1, ylab = "% annotated cells", xlab = "Cell Label", main = paste0("Cells annotated in cluster ", cluster)))
		  }
		  dev.off()
		  ## Try to assign per-cluster annotation following commands from scType vignette
		  # change NA SingleR values to "unassigned"
		  seurat@meta.data[[ref_DS]][is.na(seurat@meta.data[[ref_DS]])] <- "unassigned"
		  # create identity matrix to summarize afterwards each label per cluster
		  labels_unique <- seurat@meta.data %>% pull(ref_DS) %>% unique()
		  es.max.seu <- seurat@meta.data %>% select(ref_DS) %>% t()
		  es.max.singleR <- matrix(0, nrow = length(labels_unique), ncol = ncol(es.max.seu), dimnames = list(labels_unique, colnames(es.max.seu)))
		  for (i in 1:length(labels_unique)) {
		    es.max.singleR[i, ] <- as.numeric(es.max.seu == labels_unique[i])
		  }
		  # merge by cluster
		  cL_results_singleR = do.call("rbind", lapply(unique(seurat@meta.data[[colName]]), function(cl){
		    es.max.singleR.cl = sort(rowSums(es.max.singleR[ ,rownames(seurat@meta.data[seurat@meta.data[[colName]]==cl, ])]), decreasing = !0)
		    head(data.frame(cluster = cl, type = names(es.max.singleR.cl), scores = es.max.singleR.cl, ncells = sum(seurat@meta.data[[colName]]==cl)), 10)
		  }))
		  write.xlsx(cL_results_singleR, paste0(dir.name, "/", folders[6],  "/3-", ref_DS, "/3.4_SingleR_annotation_SummPerCluster_",ref_DS, ".xlsx"), row.names = TRUE)
		  ## and select top cluster
		  singleR_scores = cL_results_singleR %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

		  # set low-confident (<25% annotated cells to top celltype), clusters to "unknown"
		  singleR_scores$type[as.numeric(as.character(singleR_scores$scores)) < singleR_scores$ncells/5] = "Unclassified"

		  # plot cell types on UMAP
		  seurat@meta.data$customclassif = ""
		  for(j in unique(singleR_scores$cluster)){
		    cl_type_singleR = singleR_scores[singleR_scores$cluster==j,]; 
		    seurat@meta.data$customclassif[seurat@meta.data[[colName]] == j] = as.character(cl_type_singleR$type[1])
		  }
		  seurat@meta.data <- seurat@meta.data %>% rename_with(~ paste0(ref_DS, "_ClusterSum_", colName), "customclassif")
		  p3.5 <- DimPlot(seurat, reduction = "umap", label = TRUE, repel = TRUE, group.by = paste0(ref_DS, "_ClusterSum_", colName))  +
        		theme(legend.text=element_text(size=rel(0.5)))     
		  ggsave(paste0(dir.name, "/", folders[6], "/3-", ref_DS, "/3.4_SingleR_annotation_SummPerCluster_", ref_DS, ".pdf"), plot = p3.5, , width = 7.5, scale = 1.5)
		  
		  message("3.3. singleR results are plotted")
		}
	  }
	}
}

# D) Load AZIMUTH results (if enabled) to seurat object:
if (sum(grepl("9.4.1_Azimuth_annotation.tsv", list.files(path = paste0(dir.name, "/", folders[6])))) > 0){
	# if Azimuth annotation is present
    message("5. Integrating Azimuth Analysis")
	azimuth_annotation <- read.table(paste0(dir.name, "/", folders[6], "/4.1_Azimuth_annotation.tsv"), header = T, quote = "", sep = "\t")
	seurat@meta.data$Azimuth_annotation <- "Unknown"
	seurat@meta.data <- seurat@meta.data %>% rownames_to_column("row") %>% rowwise() %>%
		      mutate(Azimuth_annotation = if_else(row %in% azimuth_annotation$row,
		                                azimuth_annotation[azimuth_annotation$row == row, "predicted.celltype.l2"],
		                                Azimuth_annotation)) %>% column_to_rownames("row") %>% as.data.frame()
	message(print(table(seurat$Azimuth_annotation)))
}
saveRDS(seurat, file = paste0(dir.name, "/",folders[6], "/seurat_annotated.rds"))
message("6. Seurat object was saved.")
