#!/bin/bash
# To run: Replace the souporcell directory path and conda env name BEFORE running this script
# Usage: sh replace_paths.sh

# Function to replace the old path and conda env name in a script with the current directory
replace_paths() {
    workflow_script="$1"
    souporcell_scripts="$2"
    conda_env_name="$3"
    utils_script="$4"
    ribo_gene_path="$5"

    # Get the current directory
    new_path=$(pwd)

    # Read the script content
    script_content=$(<"$workflow_script")

    # Replace script path
    script_content=${script_content//CURRENT_DIRECTORY/$new_path}
    
    # Replace souporcell path
    script_content=${script_content//SOUPORCELL_DIRECTORY/$souporcell_scripts}

    # Replace CONDA_ENV_NAME
    script_content=${script_content//CONDA_ENV_NAME/$conda_env_name}

    # Write the updated content back to the script
    echo "$script_content" > "$workflow_script"

    # Read the utils.py content
    script_content=$(<"$utils_script")

    # Replace script path
    script_content=${script_content//RIBO_GENE_PATH_REPLACE/$ribo_gene_path}

    # Write the updated content back to the script
    echo "$script_content" > "$utils_script"
    
}

# Example script file
workflow_script="multiome_qc_workflow.sh"
utils_script="utils.py"

# Change these values accordingly
souporcell_scripts=SOUPORCELL_DIRECTORY     # example: "/data/niecr/ramaiah/multiome/analysis/s0uporcell"
conda_env_name=CONDA_ENV_NAME               # example: "multiome_winner"
ribo_gene_path=RIBO_GENE_PATH_REPLACE       # example: "/data/niecr/ramaiah/multiome/analysis/mymodule/RB_genes_human"

# Call the function to replace paths
replace_paths "$workflow_script" "$souporcell_scripts" "$conda_env_name" "$utils_script" "$ribo_gene_path"
