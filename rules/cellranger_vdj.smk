rule cellranger_input:
    input:
        samples_path=config["samples"]
    output:
        f"{OUTDIR}/cellranger/{{sample}}/csvRule_touchCheck.txt"
    log:
        f"{LOGDIR}/cellranger/{{sample}}/{{sample}}.writecsv.log"
    benchmark:
        f"{LOGDIR}/cellranger/{{sample}}/{{sample}}.writecsv.bmk"
    params:
        units_path=config["units"],
        ref_GEX=config["ref"]["cellranger_ann"],
        ref_VDJ=config["ref"]["cellranger_vdj"],
        output_dir = f"{OUTDIR}/cellranger/{{sample}}/",
        sample_i=f"{{sample}}",
        cmo_ids=config["ref"]["cmo_ids"],
        cmo_refSet=config["ref"]["cmo_reference"],
        is_cmo = lambda wildcards: is_cmo_run()
    shell:"""
    file_count=$( find "./scripts/" -name  "cellranger-7.2.0*" -type d | wc -l )
    if [[ $file_count -le 0 ]]; then
    	echo "required cellranger-7.2.0 installation"
    	exit 1
    fi
    bash ./scripts/cellranger_writeCSV.sh {params.ref_GEX} {params.ref_VDJ} {input.samples_path} {params.units_path} {params.output_dir} {params.sample_i} {params.cmo_ids} {params.cmo_refSet}  {params.is_cmo}&> {log}
    """

rule cellranger:
    input:
        f"{OUTDIR}/cellranger/{{sample}}/csvRule_touchCheck.txt"
    output:
        f"{OUTDIR}/cellranger/{{sample}}/cellranger_touchCheck.txt",
        f"{OUTDIR}/cellranger/{{sample}}/{{sample}}_cellR/outs/config.csv" if is_cmo_run() else []
    resources:
        mem_mb=get_resource("cellranger","mem_mb"),
        walltime=get_resource("cellranger","walltime")
    threads: 
        get_resource("cellranger","threads")  
    log:
        f"{LOGDIR}/cellranger/{{sample}}/{{sample}}.cellranger.log"
    benchmark:
        f"{LOGDIR}/cellranger/{{sample}}/{{sample}}.cellranger.bmk"
    params:
        samples_path=config["samples"],
        cellranger_csv=f"{OUTDIR}/cellranger/{{sample}}/multi_vdj_{{sample}}.csv",
        sample_i=f"{{sample}}",
        output_dir=f"{OUTDIR}/cellranger/{{sample}}/",
	      ref_GEX=config["ref"]["cellranger_ann"],
        units_path=config["units"],
        mem_gb = lambda wildcards, threads, resources: resources.mem_mb // 1024
    shell:"""
    outdir_final="{params.output_dir}{params.sample_i}_cellR"

    if [[ -d "$outdir_final" ]]; then
        echo "output dir detected, removing to create new pipeinstance"
        rm -rf "$outdir_final"
    fi
    file_count=$( find {params.output_dir} -name  "multi_vdj_*" -type f | wc -l )
    if [[ $file_count -gt 0 ]]; then
        echo "run cellranger multi for {wildcards.sample}"
    	./scripts/cellranger-7.2.0/cellranger multi --id out_cellranger_{params.sample_i} --csv {params.cellranger_csv} --output-dir {params.output_dir}{params.sample_i}_cellR --disable-ui --localcores {threads} --localmem {params.mem_gb} &> {log}
    else
        echo "run cellranger count"
        bash ./scripts/cellranger_GEX.sh {params.ref_GEX} {params.samples_path} {params.units_path} {params.output_dir} {params.sample_i} {threads} {params.mem_gb} &> {log}
    fi
    touch {params.output_dir}/cellranger_touchCheck.txt
    """
