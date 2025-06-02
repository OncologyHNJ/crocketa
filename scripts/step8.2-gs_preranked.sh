#!/bin/bash

# read arguments from command line
out_dir="$1"
rnk_path="$2"
selected_cond="$3"
ref_values="$4"
norm_type="$5"

echo "1: $out_dir ; 2: $rnk_path ; 3: $selected_cond ; 4: $ref_values ;"
# Convert values to array
IFS=' ' read -r -a ref_array <<< "$ref_values"
IFS=' ' read -r -a cond_array <<< "$selected_cond"

# Verify all variables are loaded
if [[ -z $out_dir || -z $ref_values || -z $selected_cond || -z $rnk_path ]]; then
    echo "Must provide all required arguments: -o output dir -i max index subpopulation -ref reference file -rnk .rnk files path"
    exit 1
fi

echo "--------------------------------------------------"
# Print first argument: output dir
echo "Selected output directory: $out_dir"

# Print second argument: rnk files dir
echo ".rnk files path: $rnk_path"

# Print third argument: ref files
for sub_ref in "${ref_array[@]}"
do
	echo "Selected reference file: $sub_ref"
done

echo "--------------------------------------------------"

# get list of rnk files from gmx_dir
for cond in "${cond_array[@]}"
    do
	echo "selected condition: $cond"
	if [[ "$norm_type" == "SCT" ]]; then
            assay_type="SCT"
        elif [[ "$norm_type" == "standard" ]]; then
            assay_type="RNA"
        fi
        echo "assay type: $assay_type"
        if [[ "$cond" =~ ^0\.[0-9]+$|^1\.0+$ ]]; then
            cond="${assay_type}_snn_res.${cond}"
            echo "modified clustering condition: $cond"
        fi
	rnk_files=($rnk_path$cond/*.rnk) # ultimo cambio, antes era: "$rnk_path""$cond"/*.rnk
	echo "rnk files are stored at: $rnk_files"
	for rnk_file in "${rnk_files[@]}"
	do
	    echo "rnk_file:  $rnk_file"
	    rnk_file_name="${rnk_file##*/}"  # extract file name (no path)
	    echo "rnk_file_name:  $rnk_file_name"
	    rnk_file_name_no_ext="${rnk_file_name%.rnk}"  # Remove .rnk extension
	    echo "rnk_file_name_no_ext:  $rnk_file_name_no_ext"
	    echo "-- run analysis for rnk file: $rnk_file_name_no_ext"
	    for sub_ref in "${ref_array[@]}"
	    do
		ref_file_name="${sub_ref##*/}"
		echo "---- run analysis for reference: $ref_file_name"
		echo "---- save results at: ${out_dir}${cond}"
		    bash scripts/GSEA_4.3.3/gsea-cli.sh GSEAPreranked -rnk ${rnk_file} -gmx $sub_ref -collapse false -mode Max_probe -norm meandiv -nperm 1000 -scoring_scheme weighted -rpt_label ${rnk_file_name_no_ext}-${ref_file_name} -include_only_symbols true -make_sets true -plot_top_x 50 -rnd_seed 123 -set_max 500 -set_min 10 -zip_report false -out ${out_dir}${cond}
	    done
	done
done
# remove daymonth folder from gsea (created by default in crocketa's main folder)
current_date=$(date +%b%d | tr '[:upper:]' '[:lower:]')
folders_to_delete=$(find . -maxdepth 1 -type d -name "*$current_date*")
for folder in $folders_to_delete; do
  rm -rf "$folder"
  echo "Folder removed: $folder"
done

touch ${out_dir}/gspreranked_ok.txt
