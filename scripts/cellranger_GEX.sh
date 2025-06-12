#!/bin/bash

## command line: bash cellranger_writeCSV.sh -gex "/Users/gsoria/Desktop/bolledit_v2/ref/gex" -vdj "/Users/gsoria/Desktop/bolledit_v2/ref/VDJ" -units "/Users/gsoria/Desktop/bolledit_v2/cellranger_draft/edited_units_forCellRnger.tsv" -o "/Users/gsoria/Desktop/bolledit_v2/cellranger_draft/" -i "s3_GeneX_C"

# read arguments from command line
ref_GEX="$1"
samples_path="$2"
units_path="$3"
out_dir="$4"
sample_i="$5"
cores="$6"
mem_mb="$7"

echo "cores: $cores"
echo "mem: $mem_mb"
# Verify all variables are loaded
if [[ -z $out_dir || -z $ref_GEX || -z $samples_path || -z $units_path || -z $sample_i ]]; then
    echo "Must provide all required arguments"
    exit 1
fi

# Create outdir if not present
if [[ ! -d "$out_dir" ]]; then
    echo "Creating output directory: $out_dir"
    mkdir -p "$out_dir"
fi

# field for gex fastq files is $3 (UNITS.TSV)
gex_fq=$(awk -v sample="$sample_i" -F '\t' '$1 == sample {print $3}' "$units_path" | sort | uniq | tail -1)
gex_fqPATH=$(dirname "$gex_fq")
echo ${gex_fqPATH}

# lanes are obtained from the very same parameter $3, as they should be the same for gex & tcr
max_lane=$(echo "$gex_fq" | grep -o 'L00[0-9]\+' | sed 's/L00//')
lanes=$(seq 1 "$max_lane" | tr '\n' '|')
lanes=${lanes%|} # remove last '|'
echo ${lanes}

# Verify all new variables are generated
if [[ -z $gex_fqPATH || -z $lanes ]]; then
    echo "Error while computing any of the variables"
    exit 1
fi

echo "All parameters OK, starting cellranger GEX..."
./scripts/cellranger-7.2.0/cellranger count --id=out_cellranger_${sample_i} --transcriptome=${ref_GEX} --fastqs=${gex_fqPATH} --sample=${sample_i} --output-dir=${out_dir}${sample_i}_cellR --localcores=${cores} --localmem=${mem_mb}
