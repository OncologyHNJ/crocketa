rule seurat_qc:
    input: 
        seurat_input
    output:
        seurat_obj= f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_pre-qc.rds",
        pre_filt_plot=report(f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/1_vlnplot_QC_variables_prefilt.pdf", caption="../report/conf/pre_filt_plot.rst", category="2_Single-cell QC"),
        gene_umi_plot=report(f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/2_geneplot_numi_vs_pctmit_ngene.pdf", caption="../report/conf/gene_umi_plot.rst", category="2_Single-cell QC")
    log:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.preqc.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.preQC.bmk"
    params:
        input_dir = lambda wc: "{}/star/{}".format(OUTDIR,wc.sample),
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        samples_path = config["samples"],
        units_path = config["units"],
        input_type = config["input_type"],
        technology = config["technology"],
        random_seed = config["random_seed"],
        case = config["case"],
        sample = f"{{sample}}",
        min_cells_per_gene = config["parameters"]["seurat_qc"]["min_cells_per_gene"]
    conda: "../envs/seurat_qc.yaml"
    resources:
        mem_mb=get_resource("seurat_qc","mem_mb"),
        walltime=get_resource("seurat_qc","walltime")
    script: 
        "../scripts/step1_qc.R"

def post_qc_input(wc):
    if wc.sample == "merged":
        return ""
    else:
        return f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_pre-qc.rds"

rule seurat_post_qc:
    input:
        seurat_obj=post_qc_input
    output:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_post-qc.rds",
        post_filt_plot=report(f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/3_vlnplot_QC_variables_postfilt.pdf", caption="../report/conf/post_filt_plot.rst", category="2_Single-cell QC"),
        stats_table=report(f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/4_pre_vs_post_stats.tsv", caption="../report/conf/stats_table.rst", category="2_Single-cell QC")
    log:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.postqc.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.postqc.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        random_seed = config["random_seed"],
        min_feat = config["parameters"]["seurat_postqc"]["min_feat"],
        max_feat = config["parameters"]["seurat_postqc"]["max_feat"],
        min_count = config["parameters"]["seurat_postqc"]["min_count"],
        max_count = config["parameters"]["seurat_postqc"]["max_count"],
        mit = config["parameters"]["seurat_postqc"]["mit_pct"],
        ribo = config["parameters"]["seurat_postqc"]["ribo_pct"],
	write_table=config["write_table"]
    conda: "../envs/seurat_qc.yaml"
    resources:
        mem_mb=get_resource("seurat_postqc","mem_mb"),
        walltime=get_resource("seurat_postqc","walltime")
    script:
        "../scripts/step2.1_postqc.R"

rule seurat_filter:
    input:
        lambda wc: f"{OUTDIR}/seurat/{wc.sample}/1_preprocessing/seurat_post-qc.rds"
    output:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_post-qc-filtered.rds"
    log:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.filter.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.filter.bmk"
    params: 
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        random_seed = config["random_seed"],
        genes = config["parameters"]["seurat_filter"]["genes"],
        filter_out = config["parameters"]["seurat_filter"]["filter_out"],
        threshold = config["parameters"]["seurat_filter"]["threshold"]
    conda: "../envs/seurat_qc.yaml"
    resources:
        mem_mb=get_resource("seurat_filter", "mem_mb"),
        walltime=get_resource("seurat_filter", "walltime")
    script:
        "../scripts/step2.2_filter.R"

def Doublet_input(wc):
    # samples = [s for s in process_samples([u.sample for u in units.itertuples()]) if s != "merged"]
    if config["parameters"]["seurat_filter"]["genes"]:
        return f"{OUTDIR}/seurat/{wc.sample}/1_preprocessing/seurat_post-qc-filtered.rds"
    else:
        return f"{OUTDIR}/seurat/{wc.sample}/1_preprocessing/seurat_post-qc.rds"

rule seurat_DoubletFinder:
    input:
        seurat_obj=Doublet_input
    output:
        f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_post-qc_singlet.rds",
        post_DF_plot=report(f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/5_vlnplot_QC_variables_postDoubletFilter.pdf", caption="../report/conf/post_doublets_plot.rst", category="2_Single-cell QC"),
        stats_table_DF=report(f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/6_post_vs_postDoublet_stats.tsv", caption="../report/conf/stats_table_DF.rst", category="2_Single-cell QC")
    log:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.DoubletFinder.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/1_preprocessing/{{sample}}.DoubletFinder.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        random_seed = config["random_seed"],
        sample_i=f"{{sample}}"
    conda: "../envs/DoubletFinder.yaml"
    resources:
        mem_mb=get_resource("seurat_DoubletFinder","mem_mb"),
        walltime=get_resource("seurat_DoubletFinder","walltime")
    script:
        "../scripts/step2.3_DoubletFinder.R"

        
rule seurat_merge:
    input:
        data=get_merge_input
    output:
        seurat_obj=f"{OUTDIR}/seurat/merged/1_preprocessing/seurat_post-qc.rds"
    log:
        f"{LOGDIR}/seurat/merged/1_preprocessing/merged.normalization.log"
    benchmark:
        f"{LOGDIR}/seurat/merged/1_preprocessing/merged.normalization.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/merged",
        random_seed = config["random_seed"], 
        velocyto = config["parameters"]["velocyto"]["enabled"],
        outdir_config = f"{OUTDIR}",
	write_table=config["write_table"]
    conda: "../envs/seurat_qc.yaml"
    resources:
        mem_mb=get_resource("seurat_merge","mem_mb"),
        walltime=get_resource("seurat_merge","walltime")
    script:
        "../scripts/step2.4_seurat_merge.R"
           
def norm_input(wc):
    if wc.sample == "integrated":
        return ""
    #if wc.sample == "merged":
    #    return ""  
    else:
        if config["parameters"]["seurat_filter"]["genes"]:
            return f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_post-qc-filtered.rds"
        else:
            return f"{OUTDIR}/seurat/{{sample}}/1_preprocessing/seurat_post-qc.rds"
            
rule seurat_normalization:
    input:
        seurat_obj=norm_input
    output:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/2_normalization/seurat_normalized-pcs.rds",
        elbow_plot=report(f"{OUTDIR}/seurat/{{sample}}/2_normalization/3_elbowplot.pdf", caption="../report/conf/elbow_plot.rst", category="3_Normalization"),
        dimplot=report(f"{OUTDIR}/seurat/{{sample}}/2_normalization/2_dimplot.pdf", caption="../report/conf/dimplot.rst", category="3_Normalization"),
        cellcycle_plot=report(f"{OUTDIR}/seurat/{{sample}}/2_normalization/6_cell_cycle_dimplot.pdf", caption="../report/conf/cellcycle_plot.rst", category="3_Normalization")
    log:
        f"{LOGDIR}/seurat/{{sample}}/2_normalization/{{sample}}.normalization.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/2_normalization/{{sample}}.normalization.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        random_seed = config["random_seed"],
        case = config["case"],
        norm_type = config["parameters"]["seurat_normalization"]["norm_type"],
        regress_out = config["parameters"]["seurat_normalization"]["regress_out"]["enabled"],
        vars_to_regress = config["parameters"]["seurat_normalization"]["regress_out"]["vars_to_regress"],
        regress_cell_cycle = config["parameters"]["seurat_normalization"]["regress_cell_cycle"],
        regress_merge_effect = config["parameters"]["seurat_normalization"]["regress_merge_effect"],
        variable_features = config["parameters"]["seurat_normalization"]["variable_features"],
        write_table=config["write_table"]
    conda: "../envs/seurat_norm.yaml"
    resources:
        mem_mb=get_resource("seurat_normalization","mem_mb"),
        walltime=get_resource("seurat_normalization","walltime")
    threads: 
        get_resource("seurat_normalization","threads")
    script:
        "../scripts/step3.1_normalization.R"

rule seurat_integration:
    input:
        data=get_integration_input_sm
    output:
        seurat_obj=f"{OUTDIR}/seurat/integrated/2_normalization/seurat_normalized-pcs.rds"
    log:
        f"{LOGDIR}/seurat/integrated/2_normalization/integrated.normalization.log"
    benchmark:
        f"{LOGDIR}/seurat/integrated/2_normalization/integrated.normalization.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/integrated",
        random_seed = config["random_seed"],
        case = config["case"],
        norm_type = config["parameters"]["seurat_normalization"]["norm_type"],
        variable_features = config["parameters"]["seurat_normalization"]["variable_features"],
        vars_to_regress = config["parameters"]["seurat_normalization"]["regress_out"]["vars_to_regress"],  
        velocyto = config["parameters"]["velocyto"]["enabled"],
        outdir_config = f"{OUTDIR}",
	write_table=config["write_table"]
    conda: "../envs/seurat_norm.yaml"
    resources:
        mem_mb=get_resource("seurat_integration","mem_mb"),
        walltime=get_resource("seurat_integration","walltime")
    threads: 
        get_resource("seurat_integration","threads")
    script:
        "../scripts/step3.2_integration.R"

rule seurat_find_clusters:
    input:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/2_normalization/seurat_normalized-pcs.rds"
    output:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/3_clustering/seurat_find-clusters.rds",
        anndata_obj=f"{OUTDIR}/seurat/{{sample}}/3_clustering/seurat_find-clusters.h5ad",
        UMAP_perAssay=report(f"{OUTDIR}/seurat/{{sample}}/3_clustering/2_umap_by_assay.pdf", caption="../report/conf/dimplot_UMAP.rst", category="4_Clustering"),
        clustree=report(f"{OUTDIR}/seurat/{{sample}}/3_clustering/1_clustree.pdf", caption="../report/conf/clustree.rst", category="4_Clustering")
    log:
        f"{LOGDIR}/seurat/{{sample}}/3_clustering/{{sample}}.find-clusters.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/3_clustering/{{sample}}.find-clusters.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        random_seed = config["random_seed"],
        pc = config["parameters"]["seurat_find_clusters"]["principal_components"],
        res = config["parameters"]["seurat_find_clusters"]["resolutions"],
        k_neighbors = config["parameters"]["seurat_find_clusters"]["k_neighbors"], 
        batch_metadata = config["parameters"]["seurat_find_clusters"]["batch_metadata"]
    conda: "../envs/seurat_clustering.yaml"
    resources:
        mem_mb=get_resource("seurat_find_clusters","mem_mb"),
        walltime=get_resource("seurat_find_clusters","walltime")
    threads: 
        get_resource("seurat_find_clusters","threads")
    script:
        "../scripts/step4_find-clusters.R"

rule scanpy_load_anndata:
    input:
        anndata_obj=f"{OUTDIR}/seurat/{{sample}}/3_clustering/seurat_find-clusters.h5ad"
    output:
        check=f"{OUTDIR}/scanpy/{{sample}}/check.txt"
    log:
        f"{LOGDIR}/scanpy/{{sample}}/check.log"
    benchmark:
        f"{LOGDIR}/scanpy/{{sample}}/check.bmk"
    params:
        output_dir = f"{OUTDIR}/scanpy/{{sample}}",
    conda: "../envs/scanpy.yaml"
    resources:
        mem_mb=get_resource("seurat_find_clusters","mem_mb"),
        walltime=get_resource("seurat_find_clusters","walltime")
    threads: 
        get_resource("seurat_find_clusters","threads")
    script:
        "../scripts/readAnnData.py"

rule seurat_scRepertoire:
    input:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/3_clustering/seurat_find-clusters.rds"
    output:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/4_scRepertoire/seurat_scRepertoire.rds",
        Identified_clones=report(f"{OUTDIR}/seurat/{{sample}}/4_scRepertoire/0.2_IdentifiedClones_Dimplot.pdf", caption="../report/conf/vdj_clones.rst", category="5_ImmuneRepertoire"),
        Clonetype_Freqs=report(f"{OUTDIR}/seurat/{{sample}}/4_scRepertoire/0.1_Clonetype_freqs_DimPlot.pdf", caption="../report/conf/CTaa_freqs.rst", category="5_ImmuneRepertoire")
    log:
        f"{LOGDIR}/seurat/{{sample}}/4_scRepertoire/{{sample}}.scRepertoire.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/4_scRepertoire/{{sample}}.scRepertoire.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        random_seed = config["random_seed"],
        ident_cond = config["parameters"]["seurat_scRepertoire"]["ident_cond"],
        input_csv = config["parameters"]["seurat_scRepertoire"]["input_csv"],
        merge = config["parameters"]["seurat_merge"]["enabled"],
        input_type = config["input_type"],
        cellranger_enabled = config["parameters"]["cellranger"]["enabled"],
        cells = config["parameters"]["seurat_scRepertoire"]["cells"],
        relative = config["parameters"]["seurat_scRepertoire"]["relative"],
        cloneTypes = config["parameters"]["seurat_scRepertoire"]["cloneTypes"],
        n_CTaa = config["parameters"]["seurat_scRepertoire"]["n_CTaa"],
        vdjdb_ref = config["parameters"]["seurat_scRepertoire"]["vdjdb_ref"],
        mcpas_ref = config["parameters"]["seurat_scRepertoire"]["mcpas_ref"],
        samples_path = config["samples"],
        cellranger_out = f"{OUTDIR}/cellranger",
        case = config["case"]
    conda: "../envs/seurat_scRepertoire.yaml"
    resources:
        mem_mb=get_resource("seurat_scRepertoire","mem_mb"),
        walltime=get_resource("seurat_scRepertoire","walltime")
    threads: 
        get_resource("seurat_scRepertoire","threads")
    script:
        "../scripts/step5_scRepertoire.R"

rule seurat_degs:
    input:
        seurat_obj=get_optional_repertoire_input
    output:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/5_degs/seurat_degs.rds"
    log:
        f"{LOGDIR}/seurat/{{sample}}/5_degs/{{sample}}.seurat_degs.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/5_degs/{{sample}}.seurat_degs.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        random_seed = config["random_seed"],
        selected_cond = config["parameters"]["seurat_degs"]["selected_cond"],
        ranking = config["parameters"]["seurat_degs"]["ranking"],
        test = config["parameters"]["seurat_degs"]["test"],
    conda: "../envs/seurat_degs.yaml"
    resources:
        mem_mb=get_resource("seurat_degs","mem_mb"),
        walltime=get_resource("seurat_degs","walltime")
    threads: 
        get_resource("seurat_degs","threads")
    script:
        "../scripts/step6_degs.R"

rule seurat_annotation_AZIMUTH:
    input:
        seurat_obj=get_optional_repertoire_input
    output:
        f"{OUTDIR}/seurat/{{sample}}/6_annotation/9.4.1_Azimuth_annotation.tsv"
    log:
        f"{LOGDIR}/seurat/{{sample}}/6_annotation/{{sample}}.Azimuth_annotation.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/6_annotation/{{sample}}.Azimuth_annotation.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        random_seed = config["random_seed"],
        ident_cond = config["parameters"]["seurat_annotation_AZIMUTH"]["ident_cond"],
        reference = config["parameters"]["seurat_annotation_AZIMUTH"]["reference"]
    conda: "../envs/seurat_annotation_AZIMUTH.yaml"
    resources:
        mem_mb=get_resource("seurat_annotation_AZIMUTH","mem_mb"),
        walltime=get_resource("seurat_annotation_AZIMUTH","walltime")
    threads: 
        get_resource("seurat_annotation_AZIMUTH","threads")
    script:
        "../scripts/step7.2_annotation_AZIMUTH.R"

rule seurat_annotation:
    input:
        seurat_obj=get_optional_repertoire_input
    output:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/6_annotation/seurat_annotated.rds",
        scType_ann=report(f"{OUTDIR}/seurat/{{sample}}/6_annotation/2.3_scType_annotation_non_restrict.pdf", caption="../report/conf/scType_ann.rst", category="6_Annotation")
    log:
        f"{LOGDIR}/seurat/{{sample}}/6_annotation/{{sample}}.annotation.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/6_annotation/{{sample}}.annotation.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        random_seed = config["random_seed"],
        sctype_gsPrepare = f"./scripts/scType-gene_sets_prepare.R",
        sctype_score = f"./scripts/scType-sctype_score.R",
        sctype_autoDetect = f"./scripts/scType-auto_detect_tissue_type.R",
        interest_feat = config["parameters"]["seurat_annotation"]["interest_feat"],
        gmx_dir = config["parameters"]["seurat_annotation"]["gmx_dir"],
        ident_cond = config["parameters"]["seurat_annotation"]["ident_cond"],
        tissue_scType = config["parameters"]["seurat_annotation"]["tissue_scType"],
        scType_ref = config["parameters"]["seurat_annotation"]["ref_scType"],
        restrict_REF_scType = config["parameters"]["seurat_annotation"]["restrict_REF_scType"],
        ref_singleR = config["parameters"]["seurat_annotation"]["ref_singleR"],
        label_singleR = config["parameters"]["seurat_annotation"]["label_singleR"],
        case = config["case"]
    conda: "../envs/seurat_annotation.yaml"
    resources:
        mem_mb=get_resource("seurat_annotation","mem_mb"),
        walltime=get_resource("seurat_annotation","walltime")
    threads: 
        get_resource("seurat_annotation","threads")
    script:
        "../scripts/step7.1_annotation.R"


def get_rnk_path_gspreranked(wc):
    if config["parameters"]["gs_preranked"]["rnk_path"] == True:
        samples = process_samples([u.sample for u in units.itertuples()])
        return expand("{OUTDIR}/seurat/{sample}/5_degs/",sample=samples,OUTDIR=OUTDIR)
    else:
        return config["parameters"]["gs_preranked"]["rnk_path"]
    
rule gs_preranked:
    input: 
        seurat_obj=get_optional_repertoire_input
    output:
        f"{OUTDIR}/seurat/{{sample}}/7_gs/gspreranked/gspreranked_touchCheck.txt"
    log:
        f"{LOGDIR}/seurat/{{sample}}/7_gs/gspreranked/{{sample}}.gspreranked.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/7_gs/gspreranked/{{sample}}.gspreranked.bmk"
    params:
        output_dir=f"{OUTDIR}/seurat/{{sample}}/7_gs/gspreranked/",
        selected_cond = config["parameters"]["seurat_degs"]["selected_cond"],
        rnk_path=get_rnk_path_gspreranked,
        ref_values=config["parameters"]["gs_preranked"]["ref_path"],
        norm_type=config["parameters"]["seurat_normalization"]["norm_type"]
        #rnk_path=f"{OUTDIR}/seurat/{{sample}}/5_degs/",
    resources:
        mem_mb=get_resource("gs_preranked","mem_mb"),
        walltime=get_resource("gs_preranked","walltime")
    threads: 
        get_resource("gs_preranked","threads")
    shell:
        """
        bash ./scripts/step8.2-gs_preranked.sh {params.output_dir} {params.rnk_path} "{params.selected_cond}" "{params.ref_values}" {params.norm_type} &> {log}
        touch {params.output_dir}/gspreranked_touchCheck.txt 
        """

rule seurat_gs:
    input:
        seurat_obj=get_optional_repertoire_input
    output:
        seurat_obj=f"{OUTDIR}/seurat/{{sample}}/7_gs/seurat_complete.rds"
    log:
        f"{LOGDIR}/seurat/{{sample}}/7_gs/{{sample}}.seurat_complete.log"
    benchmark:
        f"{LOGDIR}/seurat/{{sample}}/7_gs/{{sample}}.seurat_complete.bmk"
    params:
        output_dir = f"{OUTDIR}/seurat/{{sample}}",
        random_seed = config["random_seed"],
        resolutions = config["parameters"]["seurat_find_clusters"]["resolutions"],
        norm_type = config["parameters"]["seurat_normalization"]["norm_type"],
        gs_collection = config["parameters"]["seurat_gs"]["geneset_collection"],
        geneset_percentage = config["parameters"]["seurat_gs"]["geneset_percentage"]
    conda: "../envs/seurat_gs.yaml"
    resources:
        mem_mb=get_resource("seurat_gs","mem_mb"),
        walltime=get_resource("seurat_gs","walltime")
    script:
        "../scripts/step8.1_gs-scoring.R"

rule slingshot:
    input:
        seurat_obj=get_optional_repertoire_input
    output:
        seurat_obj=f"{OUTDIR}/slingshot/{{sample}}/8_traj_in/slingshot_sce_objects.RData"
    log:
        f"{LOGDIR}/slingshot/{{sample}}/8_traj_in/{{sample}}.slingshot.log"
    benchmark:
        f"{LOGDIR}/slingshot/{{sample}}/8_traj_in/{{sample}}.slingshot.bmk"
    params:
        output_dir = f"{OUTDIR}/slingshot/{{sample}}",
        random_seed = config["random_seed"],
        selected_res = config["parameters"]["slingshot"]["selected_res"],
        start_clus = config["parameters"]["slingshot"]["start_clus"],
        end_clus = config["parameters"]["slingshot"]["end_clus"],
        n_var_genes = config["parameters"]["slingshot"]["n_var_genes"],
        n_plotted_genes = config["parameters"]["slingshot"]["n_plotted_genes"],
        pc = config["parameters"]["seurat_find_clusters"]["principal_components"],
	graphics = config["graphics"]
    conda: "../envs/slingshot.yaml"
    resources:
        mem_mb=get_resource("slingshot","mem_mb"),
        walltime=get_resource("slingshot","walltime")
    script:
        "../scripts/step9_traj_in.R"

rule vision:
    input:
        seurat_obj=get_optional_repertoire_input
    output:
        seurat_obj=f"{OUTDIR}/vision/{{sample}}/9_func_analysis/vision_object.rds"
    log:
        f"{LOGDIR}/vision/{{sample}}/9_func_analysis/{{sample}}.vision.log"
    benchmark:
        f"{LOGDIR}/vision/{{sample}}/9_func_analysis/{{sample}}.vision.bmk"
    params:
        output_dir = f"{OUTDIR}/vision/{{sample}}",
        random_seed = config["random_seed"],
        selected_res = config["parameters"]["vision"]["selected_res"],
        geneset_collection = config["parameters"]["vision"]["geneset_collection"],
        meta_columns = config["parameters"]["vision"]["meta_columns"],
        regress_out = config["parameters"]["seurat_normalization"]["regress_out"]["enabled"],
        vars_to_regress = config["parameters"]["seurat_normalization"]["regress_out"]["vars_to_regress"],
        regress_cell_cycle = config["parameters"]["seurat_normalization"]["regress_cell_cycle"]
    conda: "../envs/vision.yaml"
    resources:
        mem_mb=get_resource("vision","mem_mb"),
        walltime=get_resource("vision","walltime")
    threads: 
        get_resource("vision","threads")
    script:
        "../scripts/step10_func_analysis.R"
