# RSEQC
#
rule fastqc:
    input:
        lambda wc: units.loc[(wc.sample,wc.unit)]['fq' + wc.read]
    output:
        html="{}/qc/fastqc/{{sample}}.{{unit}}.r{{read}}_fastqc.html".format(OUTDIR),
        zip="{}/qc/fastqc/{{sample}}.{{unit}}.r{{read}}_fastqc.zip".format(OUTDIR)
    threads: get_resource("fastqc","threads")
    resources:
        mem_mb=get_resource("fastqc","mem_mb"),
        walltime=get_resource("fastqc","walltime")
    params: 
        lambda wc: "-t {}".format(get_resource("fastqc","threads"))
    log:
        "{}/fastqc/{{sample}}.{{unit}}.r{{read}}.log".format(LOGDIR)
    benchmark:
        "{}/fastqc/{{sample}}.{{unit}}.r{{read}}.bmk".format(LOGDIR)
    wrapper:
        "0.35.0/bio/fastqc"

rule fastq_screen_indexes:
    output:
        "res/FastQ_Screen_Genomes/fastq_screen.conf"
    log:
        "{}/fastq_screen_indexes.log".format(LOGDIR)
    benchmark:
        "{}/fastq_screen_indexes.bmk".format(LOGDIR)
    conda:
        "../envs/fastq_screen.yaml"
    threads: get_resource("fastq_screen_indexes","threads")
    params:
        outdir=config["parameters"]["fastq_screen_indexes"]["outdir"]
    resources:
        mem_mb=get_resource("fastq_screen_indexes","mem_mb"),
        walltime=get_resource("fastq_screen_indexes","walltime")
    shell:"""
        fastq_screen --threads {threads} --get_genomes --outdir {params.outdir}/ &> {log}
    """

rule fastq_screen:
    input:
        lambda wc: units.loc[(wc.sample,wc.unit)]['fq' + wc.read],
        conf="{}/FastQ_Screen_Genomes/fastq_screen.conf".format(config["parameters"]["fastq_screen_indexes"]["outdir"])
    output:
        txt="{}/qc/fastq_screen/{{sample}}.{{unit}}.r{{read}}.fastq_screen.txt".format(OUTDIR),
        png="{}/qc/fastq_screen/{{sample}}.{{unit}}.r{{read}}.fastq_screen.png".format(OUTDIR)
    log:
        "{}/fastq_screen/{{sample}}.{{unit}}.r{{read}}.log".format(LOGDIR)
    benchmark:
        "{}/fastq_screen/{{sample}}.{{unit}}.r{{read}}.bmk".format(LOGDIR)
    threads: get_resource("fastq_screen","threads")
    resources:
        mem_mb=get_resource("fastq_screen","mem_mb"),
        walltime=get_resource("fastq_screen","walltime")
    params:
        fastq_screen_config="{}/FastQ_Screen_Genomes/fastq_screen.conf".format(config["parameters"]["fastq_screen_indexes"]["outdir"]),
        subset=100000,
        aligner='bowtie2'
    wrapper:
        "0.35.0/bio/fastq_screen"

rule rseqc_gtf2bed:
    input:
        config["ref"]["annotation"]
    output:
        bed="{}/qc/rseqc/annotation.bed".format(OUTDIR),
        db=temp("{}/qc/rseqc/annotation.db".format(OUTDIR))
    log:
        "{}/rseqc/rseqc_gtf2bed.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_gtf2bed.bmk".format(LOGDIR)
    conda:
        "../envs/gffutils.yaml"
    threads: get_resource("rseqc_gtf2bed","threads")
    resources:
        mem_mb=get_resource("rseqc_gtf2bed","mem_mb"),
        walltime=get_resource("rseqc_gtf2bed","walltime")
    script:
        "../scripts/gtf2bed.py"
       
rule rseqc_junction_annotation:
    input:
        bam="{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR),
        bed="{}/qc/rseqc/annotation.bed".format(OUTDIR)
    output:
        "{}/qc/rseqc/{{sample}}.junctionanno.junction.bed".format(OUTDIR)
    priority: 1
    log:
        "{}/rseqc/rseqc_junction_annotation/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_junction_annotation/{{sample}}.bmk".format(LOGDIR)
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix="{}/qc/rseqc/{{sample}}.junctionanno".format(OUTDIR)
    conda:
        "../envs/rseqc.yaml"
    threads: get_resource("rseqc_junction_annotation","threads")
    resources:
        mem_mb=get_resource("rseqc_junction_annotation","mem_mb"),
        walltime=get_resource("rseqc_junction_annotation","walltime")
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"
        
rule rseqc_junction_saturation:
    input:
        bam="{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR),
        bed="{}/qc/rseqc/annotation.bed".format(OUTDIR)
    output:
        "{}/qc/rseqc/{{sample}}.junctionsat.junctionSaturation_plot.pdf".format(OUTDIR)
    priority: 1
    log:
        "{}/rseqc/rseqc_junction_saturation/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_junction_saturation/{{sample}}.bmk".format(LOGDIR)
    params:
        extra=r"-q 255", 
        prefix="{}/qc/rseqc/{{sample}}.junctionsat".format(OUTDIR)
    conda:
        "../envs/rseqc.yaml"
    threads: get_resource("rseqc_junction_saturation","threads")
    resources:
        mem_mb=get_resource("rseqc_junction_saturation","mem_mb"),
        walltime=get_resource("rseqc_junction_saturation","walltime")
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"

rule rseqc_stat:
    input:
        "{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR),
    output:
        "{}/qc/rseqc/{{sample}}.stats.txt".format(OUTDIR)
    priority: 1
    log:
        "{}/rseqc/rseqc_stat/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_stat/{{sample}}.bmk".format(LOGDIR)
    conda:
        "../envs/rseqc.yaml"
    threads: get_resource("rseqc_stat","threads")
    resources:
        mem_mb=get_resource("rseqc_stat","mem_mb"),
        walltime=get_resource("rseqc_stat","walltime")
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"

        
rule rseqc_infer:
    input:
        bam="{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR),
        bed="{}/qc/rseqc/annotation.bed".format(OUTDIR)
    output:
        "{}/qc/rseqc/{{sample}}.infer_experiment.txt".format(OUTDIR)
    priority: 1
    log:
        "{}/rseqc/rseqc_infer/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_infer/{{sample}}.bmk".format(LOGDIR)
    conda:
        "../envs/rseqc.yaml"
    threads: get_resource("rseqc_infer","threads")
    resources:
        mem_mb=get_resource("rseqc_infer","mem_mb"),
        walltime=get_resource("rseqc_infer","walltime")
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"
        
rule rseqc_innerdis:
    input:
        bam="{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR),
        bed="{}/qc/rseqc/annotation.bed".format(OUTDIR)
    output:
        "{}/qc/rseqc/{{sample}}.inner_distance_freq.inner_distance.txt".format(OUTDIR)
    priority: 1
    log:
        "{}/rseqc/rseqc_innerdis/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_innerdis/{{sample}}.bmk".format(LOGDIR)
    params:
        prefix="{}/qc/rseqc/{{sample}}.inner_distance_freq".format(OUTDIR)
    conda:
        "../envs/rseqc.yaml"
    threads: get_resource("rseqc_innerdis","threads")
    resources:
        mem_mb=get_resource("rseqc_innerdis","mem_mb"),
        walltime=get_resource("rseqc_innerdis","walltime")
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"

rule rseqc_readdis:
    input:
        bam="{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR),
        bed="{}/qc/rseqc/annotation.bed".format(OUTDIR)
    output:
        "{}/qc/rseqc/{{sample}}.readdistribution.txt".format(OUTDIR)
    priority: 1
    log:
        "{}/rseqc/rseqc_readdis/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_readdis/{{sample}}.bmk".format(LOGDIR)
    conda:
        "../envs/rseqc.yaml"
    threads: get_resource("rseqc_readdis","threads")
    resources:
        mem_mb=get_resource("rseqc_readdis","mem_mb"),
        walltime=get_resource("rseqc_readdis","walltime")
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"

rule rseqc_readdup:
    input:
        "{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR)
    output:
        "{}/qc/rseqc/{{sample}}.readdup.DupRate_plot.pdf".format(OUTDIR)
    priority: 1
    log:
        "{}/rseqc/rseqc_readdup/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_readdup/{{sample}}.bmk".format(LOGDIR)
    params:
        prefix="{}/qc/rseqc/{{sample}}.readdup".format(OUTDIR)
    conda:
        "../envs/rseqc.yaml"
    threads: get_resource("rseqc_readdup","threads")
    resources:
        mem_mb=get_resource("rseqc_readdup","mem_mb"),
        walltime=get_resource("rseqc_readdup","walltime")
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"
        
rule rseqc_readgc:
    input:
        "{}/star/{{sample}}/Aligned.sortedByCoord.out.bam".format(OUTDIR)
    output:
        "{}/qc/rseqc/{{sample}}.readgc.GC_plot.pdf".format(OUTDIR)
    priority: 1
    log:
        "{}/rseqc/rseqc_readgc/{{sample}}.log".format(LOGDIR)
    benchmark:
        "{}/rseqc/rseqc_readgc/{{sample}}.bmk".format(LOGDIR)
    params:
        prefix="{}/qc/rseqc/{{sample}}.readgc".format(OUTDIR)
    conda:
        "../envs/rseqc.yaml"
    threads: get_resource("rseqc_readgc","threads")
    resources:
        mem_mb=get_resource("rseqc_readgc","mem_mb"),
        walltime=get_resource("rseqc_readgc","walltime")
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"

def multiqc_input(wc):
    f = expand("{OUTDIR}/qc/fastqc/{unit.sample}.{unit.unit}.r{read}_fastqc.zip", unit=units.itertuples(), read=('1','2'), OUTDIR=OUTDIR)
    f += expand("{OUTDIR}/star/{unit.sample}/Aligned.sortedByCoord.out.bam", unit=units.itertuples(), OUTDIR=OUTDIR)
    f += expand("{OUTDIR}/qc/rseqc/{unit.sample}.junctionanno.junction.bed", unit=units.itertuples(), OUTDIR=OUTDIR)
    f += expand("{OUTDIR}/qc/rseqc/{unit.sample}.junctionsat.junctionSaturation_plot.pdf", unit=units.itertuples(), OUTDIR=OUTDIR)
    f += expand("{OUTDIR}/qc/rseqc/{unit.sample}.infer_experiment.txt", unit=units.itertuples(), OUTDIR=OUTDIR)
    f += expand("{OUTDIR}/qc/rseqc/{unit.sample}.stats.txt", unit=units.itertuples(), OUTDIR=OUTDIR)
    f += expand("{OUTDIR}/qc/rseqc/{unit.sample}.inner_distance_freq.inner_distance.txt", unit=units.itertuples(), OUTDIR=OUTDIR)
    f += expand("{OUTDIR}/qc/rseqc/{unit.sample}.readdistribution.txt", unit=units.itertuples(), OUTDIR=OUTDIR)
    f += expand("{OUTDIR}/qc/rseqc/{unit.sample}.readdup.DupRate_plot.pdf", unit=units.itertuples(), OUTDIR=OUTDIR)
    f += expand("{OUTDIR}/qc/rseqc/{unit.sample}.readgc.GC_plot.pdf", unit=units.itertuples(), OUTDIR=OUTDIR)
    f += expand("{LOGDIR}/rseqc/rseqc_junction_annotation/{unit.sample}.log", unit=units.itertuples(), LOGDIR=LOGDIR)
    try:
        if config["parameters"]["fastq_screen"]["enabled"]:
            f += expand("{OUTDIR}/qc/fastq_screen/{unit.sample}.{unit.unit}.r2.fastq_screen.txt", unit=units.itertuples(), OUTDIR=OUTDIR)
    except KeyError:
        print("FASTQ_SCREEN disabled by config file. Skipping...")
    return f

rule multiqc:
    input:
        multiqc_input
    output:
        report(f"{OUTDIR}/qc/multiqc_report.html", caption="../report/conf/multiqc.rst", category="1_QC")
    params:
        config["parameters"]["multiqc"]
    benchmark:
        "{}/multiqc.bmk".format(LOGDIR)
    log:
        "{}/multiqc.log".format(LOGDIR)
    threads: get_resource("multiqc","threads")
    resources:
        mem_mb=get_resource("multiqc","mem_mb"),
        walltime=get_resource("multiqc","walltime")
    wrapper:
        "v1.22.0/bio/multiqc"
