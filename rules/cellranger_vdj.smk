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
        sample_i=f"{{sample}}"
    shell:"""
    file_count=$( find "./scripts/" -name  "cellranger-7.2.0*" -type d | wc -l )
    # echo $( find "./scripts/" -name  "cellranger-7.2.0*" -type d )
    # echo $file_count
    if [[ $file_count -le 0 ]]; then
    	echo "installing cellranger-7.2.0 for first time"
        wget -O ./scripts/cellranger-7.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.2.0.tar.gz?Expires=1738980243&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=TiB3~mGzxSIomToV7ef5SXuT2isc6vg1BF0zSWoNMyQVut9SN9et8wm~~JZIVrJGe~shgPKwahprNrDIEpdceVITeoqrC3WsbtZUuJL~fvVOcd3-NZdTpj~ZSOBucxrCGoxtZA9kmxNka1TRIGO1wcmaAm~8x2emu7x5Jlq-loko-Ex~BX1ziyuBGMrs4E-jZYzfSulWaFCVo1DCoSoS04lrprGg0O3h6c5w5dLp5HGooJmPyuUOjZunri-gb7lZ2iHR6iJc9kHFv3b-WeN4R83IYLFt-nH1Tb0uMaJd~FiQMVkiVYb9ZtTtUO95U~0CVrCjrHBoZ~~qLxVmrrQs-w__"
        tar -xzvf scripts/cellranger-7.2.0.tar.gz --directory=scripts/
    fi
    bash ./scripts/cellranger_writeCSV.sh {params.ref_GEX} {params.ref_VDJ} {input.samples_path} {params.units_path} {params.output_dir} {params.sample_i} &> {log}
    """

rule cellranger:
    input:
        f"{OUTDIR}/cellranger/{{sample}}/csvRule_touchCheck.txt"
    output:
        f"{OUTDIR}/cellranger/{{sample}}/cellranger_touchCheck.txt"
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
    resources:
        mem_mb=get_resource("cellranger","mem_mb"),
        walltime=get_resource("cellranger","walltime")
    threads:
        get_resource("cellranger","threads")
    shell:"""
    file_count=$( find {params.output_dir} -name  "multi_vdj_*" -type f | wc -l )
    echo $( find {params.output_dir} -name  "multi_vdj_*" -type f )
    echo $file_count
    if [[ $file_count -gt 0 ]]; then
        echo "run cellranger multi VDJ + GEX"
    	./scripts/cellranger-7.2.0/cellranger multi --id out_cellranger_{params.sample_i} --csv {params.cellranger_csv} --output-dir {params.output_dir}{params.sample_i}_cellR --disable-ui --localcores {threads} --localmem {resources.mem_mb // 1024} &> {log}
    else
        echo "run cellranger GEX"
        bash ./scripts/cellranger_GEX.sh {params.ref_GEX} {params.samples_path} {params.units_path} {params.output_dir} {params.sample_i} {threads} {resources.mem_mb} &> {log}
    fi
    touch {params.output_dir}/cellranger_touchCheck.txt
    """
