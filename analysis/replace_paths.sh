#!/bin/bash
# To run: Replace the souporcell directory path and conda env name BEFORE running this script
# Usage: sh replace_paths.sh

# Function to replace the old path and conda env name in a script with the current directory
replace_paths() {
    script_path="$1"
    souporcell_scripts="$2"
    conda_env_name="$3"

    # Get the current directory
    new_path=$(pwd)

    # Read the script content
    script_content=$(<"$script_path")

    # Replace script path
    script_content=${script_content//CURRENT_DIRECTORY/$new_path}
    
    # Replace souporcell path
    script_content=${script_content//SOUPORCELL_DIRECTORY/$souporcell_scripts}

    # Replace CONDA_ENV_NAME
    script_content=${script_content//CONDA_ENV_NAME/$conda_env_name}

    # Write the updated content back to the script
    echo "$script_content" > "$script_path"
}

# Example script file
script_path="multiome_qc_workflow.sh"

# Change these values accordingly
souporcell_scripts=SOUPORCELL_DIRECTORY     # example: "/data/niecr/ramaiah/multiome/analysis/suporcell"
conda_env_name=CONDA_ENV_NAME               # example: "multiome_winner"

# Call the function to replace paths
replace_paths "$script_path" "$souporcell_scripts" "$conda_env_name"
