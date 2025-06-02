#!/bin/bash

## command line: bash cellranger_writeCSV.sh -gex "/Users/gsoria/Desktop/bolledit_v2/ref/gex" -vdj "/Users/gsoria/Desktop/bolledit_v2/ref/VDJ" -units "/Users/gsoria/Desktop/bolledit_v2/cellranger_draft/edited_units_forCellRnger.tsv" -o "/Users/gsoria/Desktop/bolledit_v2/cellranger_draft/" -i "s3_GeneX_C"

# read arguments from command line
ref_GEX="$1"
ref_VDJ="$2"
samples_path="$3"
units_path="$4"
out_dir="$5"
sample_i="$6"

# Verify all variables are loaded
if [[ -z $out_dir || -z $ref_GEX || -z $ref_VDJ || -z $samples_path || -z $units_path || -z $sample_i ]]; then
    echo "Must provide all required arguments"
    exit 1
fi

# Create outdir if not present
if [[ ! -d "$out_dir" ]]; then
    echo "Creating output directory: $out_dir"
    mkdir -p "$out_dir"
fi

### VDJ-TCR
# field for vdj-T fastq files is $4 (SAMPLES.TSV)
vdj_fqPATH_T=$(awk -v sample="$sample_i" -F '\t' '$1 == sample {print $4}' "$samples_path" | sort | uniq)
echo ${vdj_fqPATH_T}
# field for vdj-T sample name is $2 (SAMPLES.TSV)
sample_vdj_T=$(awk -v sample="$sample_i" -F '\t' '$1 == sample {print $2}' "$samples_path" | sort | uniq)
echo ${sample_vdj_T}

### VDJ-BCR
# field for vdj-B fastq files is $5 (SAMPLES.TSV)
vdj_fqPATH_B=$(awk -v sample="$sample_i" -F '\t' '$1 == sample {print $5}' "$samples_path" | sort | uniq)
echo ${vdj_fqPATH_B}
# field for vdj-B sample name is $3 (SAMPLES.TSV)
sample_vdj_B=$(awk -v sample="$sample_i" -F '\t' '$1 == sample {print $3}' "$samples_path" | sort | uniq)
echo ${sample_vdj_B}

# if any TCR|BCR is not computed, column must be = NULL, then sample is commented and skipped - add '#' beforehand
prefix_T=""
if [[ -z "$sample_vdj_T" || "$sample_vdj_T" == "NULL" || -z "$vdj_fqPATH_T" || "$vdj_fqPATH_T" == "NULL" ]]; then
	prefix_T="# "
fi

prefix_B=""
if [[ -z "$sample_vdj_B" || "$sample_vdj_B" == "NULL" || -z "$vdj_fqPATH_B" || "$vdj_fqPATH_B" == "NULL" ]]; then
	prefix_B="# "
fi
echo ${prefix_B}
echo ${prefix_T}
## if both prefix are "#", then no BCR/TCR data so cellranger GEX must be executed instead of GEX+VDJ
if [[ "$prefix_T" == "# " && "$prefix_B" == "# " ]]; then
	# touch output file
	touch "${out_dir}/csvRule_touchCheck.txt"
	echo "stop multi cellranger run..."
	exit 0
fi
echo "multi vdj cellranger will be applied"
### GEX
# field for gex fastq files is $3 (UNITS.TSV)
gex_fq=$(awk -v sample="$sample_i" -F '\t' '$1 == sample {print $3}' "$units_path" | sort | uniq | tail -1)
gex_fqPATH=$(dirname "$gex_fq")
echo ${gex_fqPATH}

# lanes are obtained from the very same parameter $3, as they should be the same for gex & tcr
max_lane=$(echo "$gex_fq" | grep -o 'L00[0-9]\+' | sed 's/L00//')
lanes=$(seq 1 "$max_lane" | tr '\n' '|')
lanes=${lanes%|} # remove last '|'
echo ${lanes}

# Verify all variables are generated
if [[ -z "$gex_fqPATH" || -z "$lanes" || -z "$vdj_fqPATH_T" || -z "$vdj_fqPATH_B" || -z "$sample_vdj_T" || -z "$sample_vdj_B" ]]; then
    echo "Error while computing some variables"
    exit 1
fi

# save csv for next step
# if else -- solves no vdj problem when performing cellranger
cat <<EOL > "${out_dir}/multi_vdj_${sample_i}.csv"
# This template shows the possible cellranger multi config CSV options for analyzing Single Cell Gene Expression with Feature Barcode Technology (Antibody Capture, CRISPR Guide Capture, Cell Multiplexing, Antigen Capture), Fixed RNA Profiling, or Single Cell Immune Profiling data. 
# These options cannot be used all together - see section descriptions for detail.
# Use 'cellranger multi-template --parameters' to see descriptions of all parameters.
# Please see cellranger multi documentation for details and experimental design-specific examples at https://www.10xgenomics.com/support.

[gene-expression]
reference,${ref_GEX}
# probe-set,# Required, Fixed RNA Profiling only. 
# filter-probes,<true|false>, # Optional, Fixed RNA Profiling only. 
# r1-length,<int>
# r2-length,<int>
# chemistry,SC5P-PE
# expect-cells,<int>
# force-cells,<int>
# no-secondary,<true|false>
# no-bam,<true|false>
# check-library-compatibility,<true|false> # default to true
# target-panel,/path/to/target/panel, # Required, Targeted GEX only.
# no-target-umi-filter,<true|false>, # Optional, Targeted GEX only.
# include-introns,<true|false>
# min-assignment-confidence,<0.9>, # Optional, Cell Multiplexing only.
# cmo-set,/path/to/CMO/reference, # Optional, Cell Multiplexing only.
# barcode-sample-assignment,/path/to/barcode-sample-assignment/csv, # Optional, Cell Multiplexing only.

# [feature] # For Feature Barcode libraries only
# reference,/storage/scratch01/users/gsoria/cellranger-7.1.0/ref/
# r1-length,<int>
# r2-length,<int>

[vdj] # For TCR and BCR libraries only
reference,${ref_VDJ}
# inner-enrichment-primers,/path/to/primers
# r1-length,<int>
# r2-length,<int>

[libraries]
fastq_id,fastqs,lanes,feature_types
${sample_i},${gex_fqPATH},${lanes},Gene Expression
# Antibody1,/path/to/fastqs,Antibody Capture
# CRISPR1,path/to/CRISPR_fastqs,CRISPR Guide Capture
# CMO1,/path/to/fastqs,Multiplexing Capture, # Cell Multiplexing only
${prefix_T}${sample_vdj_T},${vdj_fqPATH_T},${lanes},VDJ-T
${prefix_B}${sample_vdj_B},${vdj_fqPATH_B},${lanes},VDJ-B
# VDJ_T1,path/to/vdj_T_fastqs,VDJ-T, # 5' Immune Profiling only
# VDJ_T_GD1,path/to/vdj_T_GD_fastqs,VDJ-T-GD, # 5' Immune Profiling only for gamma-delta TCR
# Antigen1,path/to/antigen_capture_fastqs,Antigen Capture #5' Antigen Capture only

# [antigen-specificity] # for 5' BCR/TCR Antigen Capture only
# control_id,mhc_allele
# Antigen1,AG001
# Antigen2,AG002

# [samples] # for Cell Multiplexing libraries only
# sample_id,cmo_ids
# sample1,CMO301
# sample2,CMO303

# [samples] # for Fixed RNA Profiling multiplexed libraries only
# sample_id,probe_barcode_ids,description
# sample1,BC001,Control
# sample2,BC003,Treated
EOL

# touch output file
touch "${out_dir}/csvRule_touchCheck.txt"
