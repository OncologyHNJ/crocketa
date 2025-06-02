log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

message("CONFIGURATION STEP")
# A. Parameters: 
# 1. Load libraries. 
suppressMessages(library("Seurat"))
packageVersion("Seurat")
# suppressMessages(library("SeuratDisk"))
suppressMessages(library("reticulate"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggalluvial"))
suppressMessages(library("writexl"))
suppressMessages(library("patchwork"))
suppressMessages(library("ggraph"))
suppressMessages(library("future"))
suppressMessages(library("scRepertoire"))
# suppressMessages(install.packages("dplyr", repos = "https://ftp.cixug.es/CRAN/", verbose = F))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
suppressMessages(library("ggvenn"))
suppressMessages(library("ggalluvial"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("readxl"))
suppressMessages(library("UpSetR"))
suppressMessages(library("purrr"))
# suppressMessages(library("ComplexUpset"))
### IMMUNARCH for viral annotation of clones
suppressMessages(library("immunarch"))

message("1. Libraries were loaded.")

# 2. Folder configuration. 
case = snakemake@params[["case"]]
cond = snakemake@params[["ident_cond"]]
input_file = snakemake@input[["seurat_obj"]]
input_type = snakemake@params[["input_type"]]
dir.name = snakemake@params[["output_dir"]]
# copy folders vble from find-clusters.R
folders = c("1_preprocessing", "2_normalization", "3_clustering", "4_scRepertoire", "5_degs", "6_annotation", "7_gs", "8_traj_in", "9_func_analysis", "10_RNAvelocity")
cellranger_enabled = snakemake@params[["cellranger_enabled"]]
csv_dir = snakemake@params[["input_csv"]]
merge = snakemake@params[["merge"]]
cells = snakemake@params[["cells"]]
relative = snakemake@params[["relative"]]
clonetypes = snakemake@params[["cloneTypes"]]
samples_path = snakemake@params[["samples_path"]]
samples_df = read.table(samples_path, header = T, quote = "", sep = "\t")
n = snakemake@params[["n_CTaa"]]

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
message(paste0("Old seurat object version: ", Version(seurat)))
# update seurat object for scRepertoire compatibility
seurat <- suppressMessages(UpdateSeuratObject(seurat))
message(paste0("Actual seurat object version: ", Version(seurat)))
# if seurat not merged, barcodes are yet to be reformatted:
if (!merge){
    colnames(seurat) <- paste0(seurat$orig.ident, "_", colnames(seurat))
}
assay_type <- seurat@active.assay
message("1. Seurat object was loaded.")
color_CTs <- colorRampPalette(c("#8D2063", "#ff0069", "#FCD127", "#99CCD9", "#2AC0E7"))      
color_conds <- colorRampPalette(c("#0D6BA6", "#24576F", "#A091DC", "#18B587", "#A5C360", "#F7E693", "#FCBB8C", "#FD9247", "#F0926F", "#d47170"))
color_clones <- colorRampPalette(c("skyblue2", "darkturquoise", "dodgerblue2","#246B82", "#107D6E","#336B66","#6A3D9A","#CAB2D6","orchid1", "deeppink1","#992160", "#7E0734", "#cf0000", "#FB9A99", "brown",  "#773824", "#D09283", "#FF7F00", "#FDBF6F", "yellow4", "yellow3",  "gold1", "khaki2", "green4", "green1", "palegreen2"))
# load csv files with repetoire data
#### WITHOUT VDJ R PACKAGE
# https://ncborcherding.github.io/vignettes/vignette.html
# https://www.borch.dev/uploads/vignette/vignette (v10)
##########################################################
## 1 annotation file per sample => COMPULSORY perform one cellranger multi run per sample (otherwise, one csv file is obtained for all, but no cell barcoding propperly obtained to match afterwards)
# print(samples_df$contig_csv)
if (cellranger_enabled){
  if (input_type=="fastq"){
    cellranger_dir = snakemake@params[["cellranger_out"]]
    csv_files <- list.files(path = cellranger_dir, pattern = "filtered_contig_annotations.csv", recursive = T, full.names = T)
    csv_files <- sort(csv_files[grepl("outs", csv_files)]) # should be already ordered but just in case
    message("cellranger results employed for clonetype analysis")
  }
  else {
    stop("Cellranger is enabled but input type is not fastq")
  }
} else{
    if(is.null(csv_dir)){
      stop("Cellranger is not enabled and no csv directory provided")
    } else{
      csv_files <- list.files(path = csv_dir, pattern = "filtered_contig_annotations.csv", recursive = T, full.names = T)
      csv_files <- sort(csv_files)
      message("csv_dir files employed for clonetype analysis")

    }
}
if (length(csv_files) == 0){
  stop("No contig annotation csv file loaded")
}
contig_list <- lapply(csv_files, function(file) {
  tryCatch({
    read.csv(file, header = TRUE)
  }, error = function(e) {
    message(sprintf("Error en el archivo '%s': %s", file, e$message))
    return(NULL)  # Devuelve NULL en caso de error
  })
})
names(contig_list) <- csv_files
contig_list <- Filter(Negate(is.null), contig_list)
message(colnames(seurat)[1])
message(length(contig_list)) # all files loaded to contig list?
message(contig_list[[1]][["barcode"]][1]) # format of barcoding

## modify barcodes in both object to match barcodes
if (!all(grepl("-1",colnames(seurat)))){
    for (i in seq_along(contig_list)) {
        contig_list[[i]][["barcode"]] <- gsub("-1", "",  contig_list[[i]][["barcode"]])
    }
}
message(head(contig_list[[1]][["barcode"]][1])) # format of barcoding: -1 suffix removed

# contig_list <- createHTOContigList(contig_list, 
#                                    seurat, 
#                                    group.by = "assay_name")
# message(print(head(contig_list)))
message("contig_list created")

# Merge with seurat object:
## 0.1- generate combined data
# combine(T/B)CR functions require a samples argument
# this is a vector with sample names for each contig in the list

# Combine clonotype data
# We can combine all cell types if wished (t-ab & t-gd & b-cells) in the same list (full.combined)
# There might be some overlapping among barcodes (cells identified in both TCR & BCR). These barcodes are originally ignored.
# Combined single objects are saved in order to check duplicated barcodes and see in final annotation results if the barcode
# is a Bcell or a Tcell in order to keep only the actual barcode and re-run analysis including duplicated cells.

# define samples variables according to contig list (to propperly match barcodes)
samples_vdj <-unlist(lapply(names(contig_list), function(file) {
    # print(file)
    val <- NULL
    for (ident in levels(seurat@meta.data[["orig.ident"]])){
        if (grepl(ident, file)) {
            val <- ident
            break  # Salir del bucle si se encuentra el valor
        }
    }
    if (is.null(val)) {
        message(sprintf("Correspondent value not found in file '%s'", file))
    }
    print(val)  # Imprime el valor encontrado
    return(val)
}))
print(samples_vdj)
if (is.null(cells) | all(!(cells %in% c("T-AB", "T-GD", "B-cells")))){
  stop("Define a propper cell label")
  # cell box is empty (null) or not one of the 3 possible sets of cells
  } else {
    full.combined <- NULL
    for (label in cells){ # might be more than one cell type at a time
      if (label %in% c("T-AB", "T-GD")){
      combined <- combineTCR(contig_list, 
        samples = samples_vdj, 
        ID = NULL,
        cells = label) # for T-cells
      message("combineTCR done")
      } else {
        if (label == "B-cells"){
          combined <- combineBCR(contig_list, 
            samples = samples_vdj,
            ID = NULL) # for B-cells
          message("combineBCR done")
        }
      }
      saveRDS(combined, paste0(dir.name, "/", folders[4], "/", "combined_", label, ".rds")) # save and check BCR/TCR overlapping afterwards
      # if required, remove overlapping barcodes and its line from original .csv file
      # (Check final annotation and see if barcode is actually B or T cell)
      full.combined <- c(full.combined, combined)
    }
  }
message("Clonetype Data has been integrated to seurat")

# check barcoding formats are ok for both datasets
print(colnames(seurat)[1])
print(full.combined[[1]][["barcode"]][1])
print(head(full.combined[[1]]))
for (i in seq_along(contig_list)) {
    if (sum(full.combined[[i]][["barcode"]] %in% colnames(seurat)) == 0){
        stop("No barcode matching seurat & combined repertoire data, check barcode re-formating")
    }
}

## 0.2- Combine with seurat object from previous steps:
seurat <- combineExpression(full.combined, seurat, 
                cloneCall= "aa", 
                group.by = "sample", 
                proportion = relative, 
                cloneTypes=c(Rare=clonetypes[1], Small=clonetypes[2], Medium=clonetypes[3], Large=clonetypes[4], Hyperexpanded=clonetypes[5]))
message("Seurat & Repertoire data were propperly combined.")
for (i in cond){
  if(grepl("^0|1.[0-9]$", i)){
    i <- paste0(assay_type, "_snn_res.", i)
  } 
  seurat@meta.data[[i]] <- as.factor(seurat@meta.data[[i]])
}
# scRepertoire Analysis:
message("1- starting scRepertoire Analysis")

## 0 - Dimplot according to cloneType
seurat@meta.data[["cloneType"]] <- factor(seurat@meta.data[["cloneType"]], 
                levels = c(paste0("Hyperexpanded (", clonetypes[4], " < X <= ", clonetypes[5], ")"), 
                          paste0("Large (", clonetypes[3], " < X <= ", clonetypes[4], ")"), 
                          paste0("Medium (", clonetypes[2], " < X <= ", clonetypes[3], ")"), 
                          paste0("Small (", clonetypes[1], " < X <= ", clonetypes[2], ")"), 
                          paste0("Single (", 0, " < X <= ", clonetypes[1], ")"), 
                          NA))

p0.1 <- DimPlot(seurat, group.by = "cloneType") + 
	theme(legend.text=element_text(size=rel(0.4))) +
  scale_color_manual(values = color_CTs(length(unique(seurat@meta.data[["cloneType"]]))), na.value="grey") +
  theme(plot.title = element_blank())
ggsave(paste0(dir.name, "/", folders[4], "/0.1_Clonetype_freqs_DimPlot.pdf"), plot = p0.1, scale = 1.5)
t0.1 <- table(seurat@meta.data[["cloneType"]])
write_xlsx(as.data.frame(t0.1), paste0(dir.name, "/", folders[4], "/0.1_Clonetype_abundance.xlsx"))
clonotype_filter <- !is.na(seurat$CTaa)
seurat$clonotype.detected <- "no.CT"
seurat$clonotype.detected[clonotype_filter] <- "detected.CT"
p0.2 <- DimPlot(seurat, group.by = "clonotype.detected", label = F, reduction='umap') +
	ggtitle("Identified clonotype") +
	theme(plot.title = element_text(hjust = 0.5)) +
	scale_colour_manual(values = c(color_CTs(2)))
ggsave(paste0(dir.name, "/", folders[4], "/0.2_IdentifiedClones_Dimplot.pdf"), plot = p0.2, scale = 1.5)
t0.2 <- table(seurat@meta.data[["clonotype.detected"]])
write_xlsx(as.data.frame(t0.2), paste0(dir.name, "/", folders[4], "/0.2_IdentifiedClones.xlsx"))

## Immunarch viral annotation (performed at the final step)
# load DDBBs of interest: vdjDB, mc
options(timeout=320)
# ALPHA & EL BETA annotation independently
## VDJDB
message("load vdj")
print(case)
org_vdj <- if_else(case == "uppercase", "HomoSapiens", "MusMusculus")
print(org_vdj)
# vdjdb_ALPHA = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz", "vdjdb", .species = org_vdj, .chain = "TRA")
# vdjdb_BETA = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz", "vdjdb", .species = org_vdj, .chain = "TRB")
vdjdb_ref_ = snakemake@params[["vdjdb_ref"]]
mcpas_ref_ = snakemake@params[["mcpas_ref"]]
db_ = snakemake@params[["scType_ref"]] # load edited sctype ddbb
vdjdb_ALPHA = dbLoad(vdjdb_ref_, "vdjdb", .species = org_vdj, .chain = "TRA")
vdjdb_BETA = dbLoad(vdjdb_ref_, "vdjdb", .species = org_vdj, .chain = "TRB")
## McPAS-TCR
message("load mcpas")
org_mcpas = ifelse(case == "uppercase", "Human", "Mouse")
mcpas = dbLoad(mcpas_ref_, "mcpas", .species = org_mcpas)

cond_i <- cond[1] # condition of interest for viral annotation plots
Idents(seurat) <- seurat@meta.data[[cond_i]]

totaldf <- as.data.frame(table(seurat@meta.data[[cond_i]], !is.na(seurat$CTaa))) %>% filter(Var2 == TRUE)
colnames(totaldf) <- c("condition", "CTaa", "count")
write_xlsx(totaldf, paste0(dir.name, "/", folders[4], "/1_Clonetype_by_", cond_i, ".xlsx"))

## VIRAL ANNOTATION
'
# TBAdb paired ann?
# tbadb = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/TBAdb.xlsx", "tbadb", .species = "Homo Sapiens", .chain = c("TRA-TRB"))
## Cannot access
# if TBAdb error solved, it might be possible to perform analysis
df_chains_paired <- df_chains %>% rowwise() %>%
  mutate("CTaa_paired" = paste0(CTaa_alpha, "_", CTaa_beta))
df_chains_paired <- df_chains_paired %>% 
  distinct(row, CTaa_paired, condition)
df_chains_paired <- split(as.data.frame(df_chains_paired), df_chains_paired$condition)

for (i in seq_along(df_chains_paired)){
  # df_chains[[i]] <- transform(df_chains[[i]], Clones = ave(seq(nrow(df_chains[[i]])), CTaa, FUN=length))
  df_chains_paired[[i]]$Clones <- 1
}
'
df_chains <- seurat@meta.data %>% mutate("row" = rownames(seurat@meta.data)) %>% rowwise()  %>%
        mutate("CTaa_alpha" = ifelse(is.na(CTaa), CTaa, strsplit(CTaa, "_")[[1]][1]),
          "CTaa_beta" = ifelse(is.na(CTaa), CTaa,strsplit(CTaa, "_")[[1]][-1])) %>%
        as.data.frame() # separate paired to single chain clones
df_chains <- separate_rows(df_chains, CTaa_alpha, CTaa_beta, sep = ";") # duplicate for those with 2 clones
df_chains_single <- split(as.data.frame(df_chains), df_chains[[cond_i]]) # list by condition
for (i in seq_along(df_chains_single)){
  df_chains_single[[i]]$Clones <- 1 # to sum clones, as not performed for unique but for single cells separately, set to 1 and it is count downstream
}
# VDJDB - Single
vdjdb_ann_alpha <- dbAnnotate(df_chains_single, vdjdb_ALPHA, "CTaa_alpha", "cdr3")
vdjdb_ann_alpha <- vdjdb_ann_alpha %>%
    mutate(total = rowSums(select(., -CTaa_alpha, -Samples)))
vdjdb_ann_beta <- dbAnnotate(df_chains_single, vdjdb_BETA, "CTaa_beta", "cdr3")
vdjdb_ann_beta <- vdjdb_ann_beta %>%
    mutate(total = rowSums(select(., -CTaa_beta, -Samples)))

# summarizing plot per condition & chain
summarized_vdjdb <- rbind(
  vdjdb_ann_alpha  %>%
  summarize_at(vars(-CTaa_alpha, -Samples, -total),
              funs(sum = sum)),
  vdjdb_ann_beta  %>%
  summarize_at(vars(-CTaa_beta, -Samples, -total),
              funs(sum = sum))
) %>% as.data.frame()

colnames(summarized_vdjdb) <- gsub("_sum", "", colnames(summarized_vdjdb))
summarized_vdjdb$chain <- c("alpha", "beta")
summarized_vdjdb <- melt(setDT(summarized_vdjdb))
p7.1A <- ggplot(summarized_vdjdb, aes(x = chain, y = variable)) + 
  geom_tile(aes(fill = value), height = 0.5) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 5, name = "PuRd")) +
  labs(title = "Viral annotated clones per chain", x = "TCR Chain", y = cond_i, fill = "count") +
  theme(panel.background = element_rect(fill = "white")) +
  geom_text(aes(label = value), color = "black", size = 3)  # Agregar etiquetas de texto con los valores
ggsave(paste0(dir.name, "/", folders[4], "/7.1A_ViralAnn_VDJDB_SummPerChain.pdf"), plot = p7.1A, scale = 1.5, width = 8)

# plot also relative numbers
## totaldf <- as.data.frame(table(seurat@meta.data[[cond_i]], !is.na(seurat$CTaa))) %>% filter(Var2 == TRUE)
## colnames(totaldf) <- c("condition", "CTaa", "count")
summarized_vdjdb <- summarized_vdjdb %>% rowwise() %>% 
  mutate("rel_value" = value/totaldf[totaldf[["condition"]] == variable, "count"])
p7.1B <- ggplot(summarized_vdjdb, aes(x = chain, y = variable)) + 
  geom_tile(aes(fill = rel_value), height = 0.5) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 5, name = "PuRd")) +
  labs(title = "Viral annotated clones per chain (relative to total clones per condition)", x = "TCR Chain", y = cond_i, fill = "count") +
  theme(panel.background = element_rect(fill = "white")) +
  geom_text(aes(label = round(rel_value, 3)), color = "black", size = 3)  # Agregar etiquetas de texto con los valores
ggsave(paste0(dir.name, "/", folders[4], "/7.1B_ViralAnn_VDJDB_SummPerChain_relative.pdf"), plot = p7.1B, scale = 1.5, width = 8)

# plot top 10 clones per chain
for (chain in c("alpha", "beta")){
  if (chain == "alpha"){
    vdjdb_ann <- vdjdb_ann_alpha} else{
      vdjdb_ann <- vdjdb_ann_beta}
  colnames(vdjdb_ann) <- gsub("_alpha|_beta", "", colnames(vdjdb_ann))
  vdjdb_ann <- vdjdb_ann %>%
    arrange(-total) %>% slice_head(n = 10) %>%
    pivot_longer(cols = -c(CTaa, Samples, total), names_to = cond_i, values_to = "value")
  vdjdb_ann$CTaa <- as.factor(vdjdb_ann$CTaa)
  p7.2 <- ggplot(vdjdb_ann, aes(x = reorder(CTaa, -value), y = value, fill = get(cond_i))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste0("Top-10 Abundant viral clones - ", chain, " chain"), x = "CTaa", y = "Number of cells", fill = cond_i) +
    scale_fill_manual(values = color_conds(n_distinct(vdjdb_ann[[cond_i]]))) +
    #scale_fill_manual(values =color_vector(9)[c(5:8)]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(dir.name, "/", folders[4], "/7.2_TOP10_ViralClones_VDJDB-",chain,".pdf"), plot = p7.2, scale = 1.5, width = 8)

}

# plot UMAP highlighting viral clones:
viralClones <- df_chains %>% rowwise() %>%
  mutate("alpha_viral" = ifelse(CTaa_alpha %in% vdjdb_ann_alpha$CTaa_alpha, CTaa, "FALSE"),
          "beta_viral" = ifelse(CTaa_beta %in% vdjdb_ann_beta$CTaa_beta, CTaa, "FALSE")) %>%
  mutate("paired_viral" = ifelse(alpha_viral != "FALSE" & beta_viral != "FALSE", CTaa, "FALSE")) %>% 
  select(CTaa, alpha_viral, beta_viral, paired_viral)

seurat <- highlightClonotypes(seurat, 
                                cloneCall= "aa", 
                                sequence = unique(viralClones$alpha_viral))
seurat$viral_alpha <- seurat@meta.data %>% rowwise() %>%
    mutate("is.viral" = ifelse(is.na(highlight), highlight, "alpha_viral")) %>% pull(is.viral)

seurat <- highlightClonotypes(seurat, 
                                cloneCall= "aa", 
                                sequence = unique(viralClones$beta_viral))

seurat$viral_beta <- seurat@meta.data %>% rowwise() %>%
    mutate("is.viral" = ifelse(is.na(highlight), highlight, "beta_viral")) %>% pull(is.viral)
seurat <- highlightClonotypes(seurat, 
                                cloneCall= "aa", 
                                sequence = unique(viralClones$paired_viral))
seurat$viral_paired <- seurat@meta.data %>% rowwise() %>%
    mutate("is.viral" = ifelse(is.na(highlight), highlight, "paired_viral")) %>% pull(is.viral)

pdf(paste0(dir.name, "/", folders[4], "/7.3_UMAP_viralAnn_VDJDB.pdf"), width = 7)
p7.3.1 <- DimPlot(seurat, group.by = "viral_alpha") + 
    theme(plot.title = element_blank()) +
    scale_color_manual(values = "#ED217C", na.value="grey") +
    theme(legend.text=element_text(size=rel(0.5))) +
    ggtitle("Viral annotated clones - alpha chain - VDJDB")
p7.3.2 <- DimPlot(seurat, group.by = "viral_beta") + 
    theme(plot.title = element_blank()) +
    scale_color_manual(values = "#1B998B", na.value="grey") +
    theme(legend.text=element_text(size=rel(0.5))) +
    ggtitle("Viral annotated clones - beta chain - VDJDB")
p7.3.3 <- DimPlot(seurat, group.by = "viral_paired") + 
    theme(plot.title = element_blank()) +
    scale_color_manual(values = "#0C6B9F", na.value="grey") +
    theme(legend.text=element_text(size=rel(0.5))) +
    ggtitle("Viral annotated clones - paired - VDJDB")
print(p7.3.1)
print(p7.3.2)
print(p7.3.3)
dev.off()

common_clones <- seurat@meta.data %>%
  group_by(CTaa) %>% # grouping AA clones
  summarise("count"=n_distinct(get(cond_i))) %>% # count condition in each group
  filter(count == length(unique(seurat@meta.data[[cond_i]]))) %>% select(CTaa) %>% # filter those groups that don't have alll conditions presents
  filter(!is.na(CTaa))

combined <- expression2List(seurat, 
                            split.by = cond_i)
p7.4.1 <- compareClonotypes(combined, cloneCall="aa", clonotypes = unique(viralClones$alpha_viral), graph = "alluvial") +
  scale_fill_manual(values = rep("#ED217C",length(viralClones$alpha_viral))) + NoLegend()
ggsave(paste0(dir.name, "/", folders[4], "/7.4A_Alluvial_alphaViralClonotypes_",cond_i, ".pdf"), plot = p7.4.1, scale = 1.5, width = 8)
p7.4.2 <- compareClonotypes(combined, cloneCall="aa", clonotypes = viralClones$beta_viral, graph = "alluvial") +
  scale_fill_manual(values = rep("#1B998B",length(viralClones$beta_viral))) + NoLegend()
ggsave(paste0(dir.name, "/", folders[4], "/7.4B_Alluvial_betaViralClonotypes_", cond_i, ".pdf"), plot = p7.4.2, scale = 1.5, width = 8)
p7.4.3 <- compareClonotypes(combined, cloneCall="aa", clonotypes = viralClones$paired_viral, graph = "alluvial") +
  scale_fill_manual(values = rep("#0C6B9F",length(viralClones$paired_viral))) + NoLegend()
ggsave(paste0(dir.name, "/", folders[4], "/7.4C_Alluvial_pairedViralClonotypes_", cond_i, ".pdf"), plot = p7.4.3, scale = 1.5, width = 8)

# MCPAS DDBB: Paired
mcpas_ann <- dbAnnotate(df_chains_single, mcpas, c("CTaa_alpha", "CTaa_beta"), c("CDR3.alpha.aa", "CDR3.beta.aa"))
ann_res <- list(vdjdb_ann_alpha, vdjdb_ann_beta, mcpas_ann)
names(ann_res) <- c("vdjdb_ann_alpha", "vdjdb_ann_beta", "mcpas_ann")
write_xlsx(ann_res, paste0(dir.name, "/", folders[4], "/7.5_ViralAnn.xlsx"))

## and perform now clonotype analysis with & without Viral annotated clones
seurat$viral_alpha[is.na(seurat$viral_alpha)] <- "noViral"
seurat$viral_beta[is.na(seurat$viral_beta)] <- "noViral"
seurat$viral_paired[is.na(seurat$viral_paired)] <- "noViral"
nonviral_Data <- FetchData(object = seurat, vars =c("viral_alpha", "viral_beta", "viral_paired"))
seurat_noViral <- seurat[,rowSums(nonviral_Data == "noViral") == ncol(nonviral_Data)]

message("Viral annotation completed!")

for (x.cond in cond){
  if(grepl("^0|1.[0-9]$", x.cond)){
    x.cond <- paste0(assay_type, "_snn_res.", x.cond)
  } 
  if(!(x.cond %in% colnames(seurat@meta.data))){
    stop("The specified condition is not available.")
  }
  for (set in c("Full_assay", "Non_Viral_assay")){
    message(paste0("1- Performing repertoire Analysis for ", set, " - According to ", x.cond))
    dir.create(paste0(dir.name, "/", folders[4], "/", set))
    dir.create(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond))
    if (set == "Full_assay"){
      seurat_i <- seurat
    } else{
        seurat_i <- seurat_noViral
    }
    Idents(seurat_i) <- x.cond
    ## 1- Dimplot by assay, but with the colors we'll use per each x.cond level
    p1 <- DimPlot(seurat_i, group.by = x.cond) + 
      theme(legend.text=element_text(size=rel(0.4))) +
      scale_color_manual(values = color_conds(length(unique(seurat_i@meta.data[[x.cond]])))) +
      theme(plot.title = element_blank())
    ggsave(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/1.1_Clonetype_DimPlot_by_", x.cond, ".pdf"), plot = p1, scale = 1.5)
		p1.2 <- DimPlot(seurat_i, group.by = x.cond, split.by = "clonotype.detected", label = F, reduction='umap') +
		  ggtitle("Identified Clonotype") +
			scale_color_manual(values = color_conds(length(unique(seurat_i@meta.data[[x.cond]])))) +
		  theme(plot.title = element_text(hjust = 0.5))
	      ggsave(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/1.2_IdentifiedClones_Dimplot_by_", x.cond, ".pdf"), plot = p1.2, scale = 1.5, width = 8)

		t1.2 <- table(seurat_i@meta.data[["clonotype.detected"]], seurat_i@meta.data[[x.cond]])
		write_xlsx(as.data.frame(t1.2), paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/1.2_IdentifiedClones_Dimplot_by_", x.cond, ".xlsx"))

    ' 
    # 2,3 plots from scRepertoire package are not relevant
    ## 2 - clonalOverlay
    p2 <- clonalOverlay(seurat_i, 
                  reduction = "umap", 
                  freq.cutpoint = 30, 
                  bins = 10, 
                  facet = "ident") + 
                    guides(color = "none")
    ggsave(paste0(dir.name, "/", folders[4], "/", x.cond, "/2_ClonalOverlay.pdf"), plot = p2, scale = 1.5)
    # add params for freq.cutpoint & bins (& patient/condition facet?)

    ## 3 - ClonalNetwork
    #No Identity filter
    p3 <- clonalNetwork(seurat_i, 
                  reduction = "umap", 
                  identity = "ident",
                  filter.clones = NULL,
                  filter.identity = NULL,
                  cloneCall = "aa")
    ggsave(paste0(dir.name, "/", folders[4], "/", x.cond, "/3_ClonalNetwork.pdf"), plot = p3, scale = 1.5)
    # and export table
    t3 <- clonalNetwork(seurat_i, 
                      reduction = "umap", 
                      identity = "ident",
                      filter.clones = NULL,
                      filter.identity = NULL,
                      cloneCall = "aa", exportTable = T)
    write_xlsx(t3, paste0(dir.name, "/", folders[4], "/", x.cond, "/3_ClonalNetwork.xlsx"))
    # or select a specific class in filter.identity to only plot network involving this group
    '
    # Clones identified by subgroups
    final_df <- NULL
    seurat_i@meta.data[[x.cond]] <- as.factor(seurat_i@meta.data[[x.cond]])
    for (per_C in levels(seurat_i@meta.data[[x.cond]])){
      dir.create(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/", per_C))
      # create the per_condition subsetting of seurat_i
      seurat_i_sub <- seurat_i[, which(x = FetchData(object = seurat_i, vars = x.cond) == per_C)]
      ## step0 - ADDITIONAL: save .csv files with all sequences split by alpha & beta, to perform viral ddbb enrichment analysis manually (if want to check immunarch results)
      # - VDJDB (viral)
      # - McPAS-TCR (http://friedmanlab.weizmann.ac.il/McPAS-TCR/) (pathologies in general?)
      for (i in c(1,2)){ # 1 for alpha chain; 2 for beta chain
        chain <- NULL
        id <- ifelse(i == 1, "Alpha", "Beta")
        for (pair in strsplit(seurat_i_sub$CTaa, "_")){
          chain_i <- pair[i]
          if (grepl(";", chain_i)){ chain_i <- strsplit(chain_i, ";")[[1]] }
          chain <- c(chain, chain_i)
        }
        chain <- data.frame("clones" = unique(chain[chain != "NA" & !is.na(chain)]))
        write.csv(chain, paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/", per_C, "/7.0-Clones_", per_C, "_", id, "chain.csv"), row.names = F)
      }
      
      ## 4 - highlightClonotypes: highlight the most-frequent clonotypes per condition
      TOPsequence <- seurat_i_sub@meta.data[order(seurat_i_sub@meta.data$Frequency, decreasing = TRUE),]
      write_xlsx(TOPsequence[,c("barcode", "Frequency", "CTaa", "CTnt", "CTgene", "cloneType")], paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/", per_C, "/2_cloneTypes_byFrequency_", x.cond, "_",per_C, ".xlsx"))
      for (i in seq(1,length(n))){ # in case more than 1 top_n specified in n
        TOP_n <- head(unique(TOPsequence[["CTaa"]]), n[i])
        # n value as input parameter if >2 desired clonotypes to plot
        # and plot the n most prominent aa sequences
        seurat_i <- highlightClonotypes(seurat_i, 
                                    cloneCall= "aa", 
                                    sequence = TOP_n)
        p2 <- DimPlot(seurat_i, group.by = "highlight") + 
        theme(plot.title = element_blank()) +
        scale_color_manual(values = rev(color_clones(n[i])), na.value="lightgrey") +
        theme(legend.text=element_text(size=rel(0.35)))
        ggsave(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/", per_C, "/2_highlightClonotypes_", x.cond, "_",per_C, "_TOP",n[i],".pdf"), plot = p2, scale = 1.5)
      }
      ## AND BARPLOT TOP100 CLONES
      TOP_100 <- head(unique(TOPsequence[["CTaa"]]), 100)
      df_100 <- seurat_i_sub@meta.data %>% filter(CTaa %in% TOP_100) %>% arrange(desc(Frequency))
      count_data<- table(df_100$CTaa)
      if (dim(count_data)[1] > 0){
		  count_df100 <- data.frame(CTaa = names(count_data), count = as.vector(count_data)) %>% arrange(desc(count))
		  count_df100$clone <- paste("clone", 1:dim(count_df100)[1], sep = "")
		  count_df100$clone <- factor(count_df100$clone, levels = count_df100$clone)
		  p6.1 <- ggplot(count_df100, aes(x= clone, y = count))+ 
		    geom_bar(stat = "identity", fill = "skyblue",colour='grey') +
		    ggtitle(paste0("TOP 100 Frequent clones per ", x.cond, ", subset ", per_C)) +
		    theme_minimal() +
		    theme( # remove the vertical grid lines
		        panel.grid.minor.x = element_blank() ,
		        panel.grid.major.x = element_blank() ,
		        # explicitly set the horizontal lines (or they will disappear too)
		        panel.grid.major.y = element_line(linewidth=.4, color="grey") ,
		        panel.grid.minor.y = element_line( linewidth=.1, color="grey")) +
		    xlab("Clonetype AA") +
		    ylab("Number of cells") +
		    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) 
		  count_df100$CTaa <- factor(count_df100$CTaa, levels = count_df100$CTaa)
		  p6.2 <- ggplot(count_df100, aes(x= CTaa, y = count))+ 
		    geom_bar(stat = "identity", fill = "skyblue",colour='grey') +
		    ggtitle(paste0("TOP 100 Frequent clones per ", x.cond, ", subset ", per_C)) +
		    theme_minimal() +
		    theme( # remove the vertical grid lines
		        panel.grid.minor.x = element_blank() ,
		        panel.grid.major.x = element_blank() ,
		        # explicitly set the horizontal lines (or they will disappear too)
		        panel.grid.major.y = element_line(linewidth=.4, color="grey") ,
		        panel.grid.minor.y = element_line( linewidth=.1, color="grey")) +
		    xlab("Clonetype AA") +
		    ylab("Number of cells") +
		    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))  
		  pdf(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/", per_C, "/6_TOP100_clonetypes_barplot_", x.cond, "_",per_C, ".pdf"), width = 10)
		  print(p6.1)
		  print(p6.2)
		  dev.off()
		  write_xlsx(count_df100, paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/", per_C, "/6_TOP100_clonetypes_barplot_", x.cond, "_",per_C, ".xlsx"))
		}
      ## 5- Clonetype repertoire deep analysis
      # create a dataframe to store unique values per clonetype range (plotted outside loop)
      t_i <- as.data.frame(table(seurat_i_sub@meta.data[, c("cloneType", x.cond)]))
      df_i <- seurat_i_sub@meta.data %>%
        group_by(cloneType) %>% mutate("unique_CTaa_perCloneType" = n_distinct(c(CTaa))) %>% # length(unique(CTaa)) per cloneType level
        select(x.cond, cloneType, unique_CTaa_perCloneType) %>% # already extracted the individual unique sets
        distinct() %>% inner_join(t_i, by = c("cloneType", x.cond)) %>% as.data.frame() # add the absolute number of cells (to resemble the barplot)
      final_df <- rbind(final_df, df_i)
    }
    
    col <- length(clonetypes)
    p3.1 <- occupiedscRepertoire(seurat_i, x.axis = x.cond) +
        scale_fill_manual(values = c(color_CTs(col))) 
    ggsave(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/3.1A_occupiedscRepertoire_", x.cond, ".pdf"), plot = p3.1, scale = 1.5)
    p3.2 <- occupiedscRepertoire(seurat_i, x.axis = x.cond, proportion = TRUE) +
        scale_fill_manual(values = c(color_CTs(col)))
    ggsave(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/3.1B_occupiedscRepertoire_", x.cond, "_relative.pdf"), plot = p3.2, scale = 1.5)
    t3.1 <- occupiedscRepertoire(seurat_i, x.axis = x.cond, exportTable = T)
    write_xlsx(t3.1, paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/3.1A_occupiedscRepertoire_", x.cond, ".xlsx"))
    ## run the same plot but with unique absolute values per condition instead of number of cells as obtained above (final_df)
    write_xlsx(as.data.frame(final_df), paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/3.1B_occupiedscRepertoire_", x.cond, "_Unique.xlsx"))
    p3.3 <- ggplot(final_df, aes(x = final_df[,x.cond], y = Freq, fill = cloneType, label = unique_CTaa_perCloneType)) + 
        geom_bar(stat = "identity") + 
        geom_text(size = 3, position = position_stack(vjust = 0.5)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_fill_manual(values = c(color_CTs(col))) + 
            ylab("Single Cells") + 
            theme_classic() + 
            theme(axis.title.x = element_blank())
    ggsave(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/3.1C_occupiedscRepertoire_", x.cond, "_Unique.pdf"), plot = p3.3, scale = 1.5)
    ## additional clonal plots
    combined2 <- expression2List(seurat_i, 
                              split.by = x.cond)
    '
    # 10.2 Clonal Homeostasis
    p10_2 <- clonalHomeostasis(combined2, 
                      cloneCall = "aa") +
        scale_fill_manual(values = c(color_vector(col)))
    ggsave(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/3.2_ClonalHomeostasis_", x.cond, ".pdf"), plot = p10_2, scale = 1.5)
    t10_2 <- clonalHomeostasis(combined2, 
                      cloneCall = "aa",
                      exportTable = T)
    write_xlsx(as.data.frame(t10_2), paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/3.2_ClonalHomeostasis_", x.cond, ".xlsx"))
    '
    # 10.3 Clonal Proportion
    p4 <- clonalProportion(combined2, 
                    cloneCall = "aa") +
        scale_fill_manual(values = c(color_CTs(col+1)))
    ggsave(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/4_ClonalProportion_", x.cond, ".pdf"), plot = p4, scale = 1.5)
    t4 <- clonalProportion(combined2, 
                    cloneCall = "aa",
                    exportTable = T)
    write_xlsx(as.data.frame(t4), paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/4_ClonalProportion_", x.cond, ".xlsx"))
    # 10.4 Clonal Overlap
    pdf(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/5.1_ClonalOverlap_", x.cond, ".pdf"), width = 8.5)
    print(clonalOverlap(combined2, 
              cloneCall="aa", 
              method="raw") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))) 
    print(clonalOverlap(combined2, 
                  cloneCall="aa", 
                  method="morisita")+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))) 
    print(clonalOverlap(combined2, 
                  cloneCall="aa", 
                  method="jaccard")+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))) 
    print(clonalOverlap(combined2, 
                  cloneCall="aa", 
                  method="overlap")+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))) 
    dev.off()
    
    # and export table
    t5.1 <- clonalOverlap(combined2, 
                  cloneCall="aa", 
                  method="morisita",
                  exportTable = T)
    write_xlsx(as.data.frame(t5.1), paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/5.1_ClonalOverlap_morisita_", x.cond, ".xlsx"))

    # additional Venn Diagram (expected to be similar to clonal overlap raw) - máx 4
    venn_list <- vector("list", n_distinct(seurat_i@meta.data[[x.cond]]))
    names(venn_list) <- sort(unique(seurat_i@meta.data[[x.cond]]))
    for (i in 1:n_distinct(seurat_i@meta.data[[x.cond]])){ 
      venn_list[[i]] <- seurat_i[, which(x = FetchData(object = seurat_i, vars = x.cond) == names(venn_list)[i])]@meta.data %>% filter(!is.na(CTaa)) %>% pull(CTaa)
    }
    if (length(unique(seurat_i@meta.data[[x.cond]])) <= 4){ # - máx 4 lists to draw
      p5.2A <- ggvenn(venn_list,
            fill_color = color_conds(n_distinct(seurat_i@meta.data[[x.cond]])),
            stroke_size = 0.5, set_name_size = 4) + 
        ggtitle(paste0("Common & unique Clones according to ", x.cond)) +
        theme(plot.title = element_text(hjust = 0.5))
      ggsave(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/5.2A_Venn-ClonalOverlap_", x.cond, ".pdf"), plot = p5.2A, scale = 1.5)
    }

    # plot 2 vs. 2 venn diagrams
    combination2v <- combn(names(venn_list), 2)
    colv2 <- color_conds(n_distinct(seurat_i@meta.data[[x.cond]]))
    names(colv2) <- sort(unique(seurat_i@meta.data[[x.cond]]))
    pdf(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/5.2B_Venn-1vs1-ClonalOverlap_", x.cond, ".pdf"))
    for (i in 1:ncol(combination2v)) {
      set1 <- venn_list[[combination2v[1, i]]]
      set2 <- venn_list[[combination2v[2, i]]]
      list_i <- list(set1, set2)
      names(list_i) <- c(combination2v[1, i], combination2v[2, i])
      pi <- ggvenn(list_i,
            fill_color = as.vector(colv2[names(list_i)]),
            stroke_size = 0.5, set_name_size = 4) + 
        ggtitle(paste0("Common & unique Clones according to ", x.cond)) +
        theme(plot.title = element_text(hjust = 0.5))
      print(pi)
      }
    dev.off()
      
    # venn Diagram for top100 clones per condition
    venn_list_t100 <- vector("list", n_distinct(seurat_i@meta.data[[x.cond]]))
    names(venn_list_t100) <- sort(unique(seurat_i@meta.data[[x.cond]]))
    file_Names <- list.files(path = paste0(dir.name, "/", folders[4], "/", set, "/", x.cond), pattern = "6_TOP100_clonetypes_barplot.*xlsx$", recursive = T, full.names = T)
    files <- lapply(file_Names, read_excel, sheet = 1)
    names(files) <- dirname(list.files(path = paste0(dir.name, "/", folders[4], "/", set, "/", x.cond), pattern = "6_TOP100_clonetypes_barplot.*xlsx$",recursive = T))

    for (i in 1:length(names(files))){ 
      venn_list_t100[[i]] <- seurat_i[, which(x = FetchData(object = seurat_i, vars = x.cond) == names(venn_list_t100)[i])]@meta.data %>% filter(CTaa %in% files[[i]][["CTaa"]]) %>% pull(CTaa)
    }
    if (length(unique(seurat_i@meta.data[[x.cond]])) <= 4){ # - máx 4 lists to draw
      p5.2C <- ggvenn(venn_list_t100,
            fill_color = color_conds(n_distinct(seurat_i@meta.data[[x.cond]])),
            stroke_size = 0.5, set_name_size = 4) + 
        ggtitle(paste0("Top 100 most frequent clones per ", x.cond, " - Common & unique Clones")) +
        theme(plot.title = element_text(hjust = 0.5))
      ggsave(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/5.2C_Venn-Top100_ClonalOverlap_", x.cond, ".pdf"), plot = p5.2C, scale = 1.5)
      t5.2C <- list("common_clones" = Reduce(intersect, venn_list), "common_clones_t100" = Reduce(intersect, venn_list_t100))
      sink(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/5.2_VennCommonClonetypes_", x.cond, ".txt"))
      print(t5.2C)
      sink()
      }
    # add here venn diagrams with number of cells with commmon clones (same but now cells instead of unique clones to check absolute overlapping)
    
    
    # highlight common clones to all conditions
    common_clones <- seurat_i@meta.data %>%
      group_by(CTaa) %>% # grouping AA clones
      summarise("count"=n_distinct(get(x.cond))) %>% # count condition in each group
      filter(count == length(unique(seurat_i@meta.data[[x.cond]]))) %>% select(CTaa) %>% # filter those groups that don't have alll conditions presents
      filter(!is.na(CTaa))
    seurat_i <- highlightClonotypes(seurat_i, 
                                  cloneCall= "aa", 
                                  sequence = common_clones$CTaa)
    seurat_i$common.clones <- seurat_i@meta.data %>% rowwise() %>%
      mutate("is.common" = ifelse(is.na(highlight), highlight, "common_clones")) %>% pull(is.common)
    p5.3A <- DimPlot(seurat_i, group.by = "highlight") + 
      theme(plot.title = element_blank()) +
      scale_color_manual(values = color_clones(length(common_clones$CTaa)), na.value="grey") +
      theme(legend.text=element_text(size=rel(0.5)))
    p5.3B <- DimPlot(seurat_i, group.by = "common.clones") + 
      theme(plot.title = element_blank()) +
      scale_color_manual(values = "red", na.value="grey") 

    ggsave(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/5.3A_highlightCommonClonotypes_", x.cond, ".pdf"), plot = p5.3A, scale = 1.5)
    ggsave(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/5.3B_highlightCommonClonotypes_", x.cond, ".pdf"), plot = p5.3B, scale = 1.5)
    # create sankey diagram with ggplot
    sankey_df <- seurat_i@meta.data %>% filter(CTaa %in% common_clones$CTaa) %>% distinct(CTaa, get(x.cond), cloneType) %>% tibble()
    colnames(sankey_df) <- c("CTaa", x.cond, "cloneType")
    write_xlsx(sankey_df, paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/5.4_CommoncloneTypes_by", x.cond, ".xlsx"))
    
    p5.4 <- compareClonotypes(combined2, cloneCall="aa", clonotypes = pull(common_clones, CTaa), graph = "alluvial") +
      scale_fill_manual(values = color_clones(length(common_clones$CTaa)))
    ggsave(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/5.4A_AlluvialCommonClonotypes_", x.cond, ".pdf"), plot = p5.4, scale = 1.5, width = 8)
    
    ## upsetR plot (similar to venn diagrams)
    message("Creating UpSetR plot")
    all_data <- do.call(rbind.data.frame, combined2) %>% select(CTaa, x.cond)


      data_wide <-
        all_data %>%
        group_by(CTaa) %>%
        count(get(x.cond)) %>% ## count how many clones per condition
        mutate(match = 1,  ## boolean for a count at each condition
            Total = sum(n))
      
      colnames(data_wide) <- c("CTaa", x.cond, "n", "match", "Total")
      data_wide <- data_wide %>% 
        pivot_wider(
          id_cols = c(CTaa, Total),
          names_from = x.cond,
          values_from = match,
          values_fill = list(match = 0)
        )  ## fill empty cells with 0
      level_present <- levels(seurat_i@meta.data[[x.cond]])[levels(seurat_i@meta.data[[x.cond]]) %in% colnames(data_wide)]
      print(level_present)
      ylabel <- "Number of Common Clones"
      xlabel <- "Total Number of Unique Clones"
      # Convert presence-matching columns in factor
      data_wide <- data_wide %>% select(-Total) %>% as.data.frame()
      print(colnames(data_wide))
      data_wide <- data_wide[, c("CTaa", level_present)]
      # Create graphic: UpSetR
      print(dim(common_clones)[1])
      print(common_clones)
      pdf(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/5.2D-upsetR_intersect.pdf"), width = 15, height = 10)
      if (dim(common_clones)[1] > 2){ # if more than 2 clones, highlight common
          message("highlighting")
          print(upset(
              data_wide,
              sets = level_present,
              keep.order=T,
              sets.bar.color = color_conds(length(level_present)),
              order.by="freq",
              text.scale = c(1.3, 1.3, 1.3, 1, 1.5, 2),
              mainbar.y.label = ylabel,
              sets.x.label = xlabel,
              queries = list(
                list(
                query = intersects,
                params = as.list(level_present),
                active = T,
                query.name = paste0("All levels - ", x.cond),
                color = "orange"
                ) #, # can add colors to any group of interactions of interest
                # list(
                # query = intersects,
                # params = list(sort(unique(seurat_i@meta.data[[x.cond]]))[1], sort(unique(seurat_i@meta.data[[x.cond]]))[2]),
                # active = T,
                # query.name = "M0 vs",
                # color = "red"
                # )'
              )
          ))
      } else{
          print(upset(
              data_wide,
              sets = level_present,
              keep.order=T,
              sets.bar.color = color_conds(length(level_present)),
              order.by="freq",
              text.scale = c(1.3, 1.3, 1.3, 1, 1.5, 2),
              mainbar.y.label = ylabel,
              sets.x.label = xlabel,
          ))
      }
      dev.off()
      message("done upsetR")
      write_xlsx(data_wide, paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/5.2D-upsetR_intersect.xlsx"))
    
    
    message("Extract number of cells per clonal overlap - Barplot")
    # extract number of cells per clonal overlap according to x.cond
	  clon_cond <- all_data %>%
			distinct(CTaa, !!sym(x.cond)) %>%
			group_by(CTaa) %>%
			summarise(Conditions = list(sort(unique(!!sym(x.cond))))) %>%
			ungroup()

		all_cond <- sort(unique(all_data[[x.cond]]))

		combinations <- unlist(lapply(1:length(all_cond), 
				                     function(x) combn(all_cond, x, simplify = FALSE)), recursive = FALSE)

		clon_in_comb <- function(comb, clon_cond) {
			sapply(clon_cond$Conditions, function(cnds) all(comb %in% cnds))
		}

		results <- map_df(combinations, function(comb) {
			presence <- clon_cond$CTaa[clon_in_comb(comb, clon_cond)]
			n_cells <- all_data %>% filter(CTaa %in% presence) %>% nrow()
			tibble(Combination = paste(comb, collapse = "&"), 
			 Num_Cells = n_cells,
			 Num_Conditions = length(comb))
		})  %>% arrange(desc(Num_Cells))
		results <- results %>%
			mutate(Combination = as.character(Combination)) %>%
			arrange(Num_Conditions, desc(Num_Cells)) %>%
			mutate(Combination = factor(Combination, levels = unique(Combination)))
		bar_width <- max(0.2,min(0.9, 5/length(results$Combination))
		
		p5.5 <- ggplot(results, aes(x = Combination, y = Num_Cells, fill = Num_Cells)) +
			geom_col(width = bar_width, color = "black", fill = "#006699") +
			labs(title = paste0("Number of cells per Intersection, according to ", x.cond, " intersection"),
			     y = "Intersection level",
			     x = "Number of cells") +
			theme_minimal(base_size = 10) + 
			theme(
			  axis.text.y = element_text(angle = 45, hjust = 1, vjust = 1),
			  legend.position = "none",
			  panel.grid.major.y = element_blank(),
			  panel.grid.minor = element_blank(),
			  plot.title = element_text(face = "bold"),
			  plot.caption = element_text(size = 8, face = "italic")) + 
			geom_text(aes(label = Num_Cells), 
			  vjust = -0.5, 
			  size = 4)
	  ggsave(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/5.2E-Intersect_NCells_per", x.cond, ".pdf"), plot = p5.5, scale = 1.5, width = 8)
	  
    # Draw scatterplots zoomed/scaled
    samples_list <- sort(names(combined2))
    for (tag in c("plot", "plot_zoomed")){
	    pdf(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/8-scatter_compare-1vs1_scaled", tag, ".pdf"))
	    for (i in seq_along(combined2)){
	        print(i)
	        for (j in c(seq_along(combined2)[c(i:length(combined2))])){
	        if (all(i != j, i != samples_list[length(samples_list)])){
	            print(j)
	            xline <- combined2[[i]] %>% 
	                    count(CTaa) %>% mutate(frequency = n / sum(n)) %>% 
	                    arrange(desc(frequency)) %>% slice_head(n=100) %>% 
	                    tail(1) %>% pull(frequency)
	            yline <- combined2[[j]] %>% 
	                    count(CTaa) %>% mutate(frequency = n / sum(n)) %>% 
	                    arrange(desc(frequency)) %>% slice_head(n=100) %>% 
	                    tail(1) %>% pull(frequency)
	            lim <- max(unlist(lapply(combined2, function(x){
	                z <- x %>% count(CTaa) %>% mutate(frequency = n / sum(n)) %>% pull(frequency) %>% max()
	                return(z)})))
	            
	            print(xline)
	            print(yline)
	            if (tag == "plot"){
	                print(scatterClonotype(combined2, 
	                    cloneCall ="aa", 
	                    x.axis = samples_list[i],
	                    y.axis = samples_list[j],
	                    dot.size = "total",
	                    graph = "proportion") + scale_color_manual(values = c("#8B6914", "#336B66","#7EC0EE","#CDCD00", "#90EE90")) +
	                        geom_vline(xintercept = xline,col="red", linetype = "dashed") + 
	                        geom_hline(yintercept= yline, col="red", linetype = "dashed") +
	                        coord_cartesian(xlim = c(0,lim), ylim=c(0,lim)) + scale_size_continuous(limits=c(1,4000)) # adjust coords
	                )
	            } else {
	                print(scatterClonotype(combined2, 
	                    cloneCall ="aa", 
	                    x.axis = samples_list[i],
	                    y.axis = samples_list[j],
	                    dot.size = "total",
	                    graph = "proportion") + scale_color_manual(values = c("#8B6914", "#336B66","#7EC0EE","#CDCD00", "#90EE90")) +
	                        geom_vline(xintercept = xline,col="red", linetype = "dashed") + 
	                        geom_hline(yintercept= yline, col="red", linetype = "dashed")
	                )
	            }
	            scatter_table <- scatterClonotype(combined2, 
	                        cloneCall ="aa", 
	                        x.axis = samples_list[i],
	                        y.axis = samples_list[j],
	                        dot.size = "total",
	                        graph = "proportion", exportTable=T)
	            write.table(scatter_table, paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/8-scatter_compare-", samples_list[i], ".vs.", samples_list[j],".tsv"), row.names = F, quote=F, sep = "\t")
		        }
			} 
		}
		dev.off()    
		}
    message("Scatterplots drawn")
    # plot HEATMAPS: draw top_n clones frequency per subpopulations in cond vector
    if (length(cond) > 1){
        message("Draw Heatmaps")
				for (subset in cond[cond != x.cond]) {
						if(grepl("^0|1.[0-9]$", subset)){
						          subset <- paste0(assay_type, "_snn_res.", subset)
						      } 
						      if(!(subset %in% colnames(seurat@meta.data))){
						          stop("The specified condition is not available.")
						      }
						      if(subset != x.cond){
											seurat_i@meta.data[[x.cond]] <- as.factor(seurat_i@meta.data[[x.cond]])
											seurat_i@meta.data[[subset]] <- as.factor(seurat_i@meta.data[[subset]])
											coul <- brewer.pal(n = 6, "PuRd")
											totaldf <- as.data.frame(table(seurat_i@meta.data[[subset]], !is.na(seurat_i$CTaa))) %>% filter(Var2 == TRUE)
											colnames(totaldf) <- c("subset", "CTaa", "maxcount")
											totaldf$maxcount <- as.numeric(totaldf$maxcount)
											for (i in seq(1,length(n))){ # in case more than 1 top_n specified in n
												TOPsums <- NULL
												for (per_C in levels(seurat_i@meta.data[[x.cond]])){
													print(sum(grepl("6",  list.files(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/", per_C)))))
													if (sum(grepl("6",  list.files(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/", per_C)))) > 1){
															full_ctaa <- read_excel(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/", per_C, "/6_TOP100_clonetypes_barplot_", x.cond, "_", per_C, ".xlsx"))
															top_ctaa <- full_ctaa %>%
															head(n[i]) %>% pull(CTaa)
															# top_ctaa$freq <- top_ctaa$Count/sum(full_ctaa)
															seurat_sub <- seurat_i[, which(x = FetchData(object = seurat_i, vars = x.cond) == per_C)]
															df_sum <- as.data.frame(table(seurat_sub$CTaa, seurat_sub@meta.data[[subset]]))
															colnames(df_sum) <- c("clone", "subset", "Count")
															df_sum <- df_sum %>% filter(clone %in% top_ctaa)
															df_sum$condition <- per_C
															TOPsums <- rbind(TOPsums, df_sum)   
															}
													}
													TOPsums <- TOPsums %>% 
															group_by(subset, condition) %>%
															summarise(total_count = sum(Count), .groups = 'drop') %>% rowwise() %>%
															mutate("rel_value" = total_count/totaldf[totaldf[["subset"]] == subset, "maxcount"]) %>%
															mutate(rel_value = ifelse(is.na(rel_value), 0, rel_value)) %>% as.data.frame()
													write.table(TOPsums, paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/2.B-heatmap_TopClones-per", x.cond, "_by_", subset, "_TOP",n[i], ".tsv"), row.names = F, quote=F, sep = "\t")
													p1 <- ggplot(TOPsums, aes(x = subset, y = condition, fill = total_count)) +
															geom_tile() +
															theme_bw(base_size = 14) +
															ylab("condition") + xlab(subset) +
															ggtitle(paste0("TOP ", n[i] , " Clones per condition across cell subsets")) +
															theme(plot.title = element_text(size=8)) +
															scale_fill_gradientn(limits=c(0, max(TOPsums$total_count)), colours = coul) +
															theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
													p2 <- ggplot(TOPsums, aes(x = subset, y = condition, fill = rel_value)) +
															geom_tile() +
															theme_bw(base_size = 14) +
															ylab("condition") + xlab(subset) +
															ggtitle(paste0("TOP ", n[i] , " Clones per condition across cell subsets")) +
															theme(plot.title = element_text(size=8)) +
															scale_fill_gradientn(limits=c(0, max(TOPsums$rel_value)), colours = coul) +
															theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
													ggsave(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/2.B-heatmap_TopClones-per", x.cond, "_by_", subset, "_TOP",n[i],"_absolute.pdf"), plot = p1, scale = 1.5)
													ggsave(paste0(dir.name, "/", folders[4], "/", set, "/", x.cond, "/2.C-heatmap_TopClones-per", x.cond, "_by_", subset, "_TOP",n[i],"_relative.pdf"), plot = p2, scale = 1.5)
										}
									}
								}
		}
  }
}
message("ANALYSIS FINISHED")
'
# 7.8. Save Seurat objects in AnnData
SaveH5Seurat(seurat, filename = paste0(dir.name, "/", folders[3], "//seurat_scRepertoire.h5Seurat"))
Convert(paste0(dir.name, "/",folders[3], "//seurat_scRepertoire.h5Seurat"), dest = "h5ad")
' # no compatibility for seurat v4.3 & seuratdisk

saveRDS(seurat_noViral, file = paste0(dir.name, "/",folders[4], "/seurat_scRepertoire-noViral.rds"))
saveRDS(seurat, file = paste0(dir.name, "/",folders[4], "/seurat_scRepertoire.rds"))

message("5. Seurat object was saved.")
