#!/bin/bash

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --all_samples)
            all_samples=true
            shift
            ;;
        --include_samples)
            include_samples="$2"
            shift 2
            ;;
        --skip_samples)
            skip_samples="$2"
            shift 2
            ;;
        --rna)
            rna=true
            shift
            ;;
        --atac)
            atac=true
            shift
            ;;
        --souporcell)
            souporcell=true
            shift
            ;;
        --overwrite_rna)
            overwrite_rna=true
            shift
            ;;
        -n|--n_genotypes)
            n_genotypes="$2"
            shift 2
            ;;
        -i|--input)
            input_directory="$2"
            shift 2
            ;;
        -o|--output)
            output_directory="$2"
            shift 2
            ;;
        -d|--dd_file)
            dd_file="$2"
            shift 2
            ;;
        -r|--ref_genome)
            ref_genome="$2"
            shift 2
            ;;
        -h|--help)
            help=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Display help message
if [[ $help == true ]]; then
    echo "Usage: sh multiome_qc_workflow.sh [OPTIONS]"
    echo "Example1: sh multiome_qc_workflow.sh --all_samples --rna --atac -i 'input_directory' -o 'output_directory'"
    echo "Example2: sh multiome_qc_workflow.sh --include_samples sample1,sample2,sample3 --rna -i 'input_directory' -o 'output_directory'"
    echo "Options:"
    echo "  --all_samples      Process all samples"
    echo "  --include_samples  Include specific samples (comma-separated)"
    echo "  --skip_samples     Skip specific samples (comma-separated)"
    echo "  --rna              Process RNA data"
    echo "  --atac             Process ATAC data"
    echo "  --souporcell       Run souporcell doublet detection on RNA data"
    echo "  --overwrite_rna    Overwrite RNA preprocessing results."
    echo "  -i, --input        Input directory (required)"
    echo "  -o, --output       Output directory (required)"
    echo "  -d, --dd_file      Input file for souporcell double detection (required to run souporcell)"
    echo "  -r, --ref_genome   Refernce genome file for doublet detection (genome.fa)"
    exit 0
fi

# Check for mandatory arguments
if [[ -z $input_directory || -z $output_directory ]]; then
    echo "Error: Input and output directories are required. Refer --help function"
    exit 1
fi

# Check if RNA or ATAC workflow is specified
if [[ $rna != true && $atac != true && $souporcell != true ]]; then
    echo "Error: RNA or ATAC or souporcell workflow must be specified. Refer --help function"
    exit 1
fi

# Function to submit LSF job
submit_job() {
    sample_name="$1"
    script_type="$2"
    input_directory="$3"
    output_directory="$4"
    
    # Find cellranger output directory
    if [[ -d "${input_directory}/${sample_name}/raw_feature_bc_matrix" ]]; then
        cr_outs="${input_directory}/${sample_name}"
    elif [[ -d "${input_directory}/${sample_name}/outs/raw_feature_bc_matrix" ]]; then
        cr_outs="${input_directory}/${sample_name}/outs"
    else
        echo "Error: raw_feature_bc_matrix directory not found for sample $sample_name"
        return 1
    fi
    
    #-------------------------------------------------#
                 # PRE-PROCESS RNA DATA
    #-------------------------------------------------#
    if [[ $script_type == "RNA" ]]; then
        echo "Submitting job for $sample_name (type: $script_type)"
        
        overwrite_rna="$5"
        if [[ "$overwrite_rna" == true ]]; then
            overwrite_flag="--overwrite"
        else
            overwrite_flag=""
        fi
    
        # Create a temporary LSF script file for RNA
        temp_script=$(mktemp)
        cat > $temp_script <<EOF
#!/bin/bash
#BSUB -J "rna_qc_$sample_name"
#BSUB -q gpuqueue
#BSUB -W 3:00
#BSUB -n 2
#BSUB -R "rusage[mem=25]"
#BSUB -e "$output_directory/$sample_name/%J.err"
#BSUB -o "$output_directory/$sample_name/%J.out"

source ~/.bashrc
cd /data/niecr/ramaiar/multiome_lucyrepo/analysis/
conda run -n multiome_winner_test3 python rna_init.py -s "$sample_name" -i "$cr_outs" -o "$output_directory" $overwrite_flag
EOF

        # Submit the job using the temporary script file
        bsub < $temp_script

        # Remove the temporary script file
        rm $temp_script

    #-------------------------------------------------#
                # PRE-PROCESS ATAC DATA
    #-------------------------------------------------#
    elif [[ $script_type == "ATAC" ]]; then
        echo "Submitting job for $sample_name (type: $script_type)"
    
        # Create a temporary LSF script file for ATAC
        temp_script=$(mktemp)
        cat > $temp_script <<EOF
#!/bin/bash
#BSUB -J "atac_qc_$sample_name"
#BSUB -q gpuqueue
#BSUB -W 3:00
#BSUB -n 2
#BSUB -R "rusage[mem=25]"
#BSUB -e "$output_directory/$sample_name/%J.err"
#BSUB -o "$output_directory/$sample_name/%J.out"

source ~/.bashrc
cd /data/niecr/ramaiar/multiome_lucyrepo/analysis/
conda run -n multiome_winner_test3 Rscript atac_init.R "$sample_name" "$cr_outs" "$output_directory"
EOF

        # Submit the job using the temporary script file
        bsub < $temp_script

        # Remove the temporary script file
        rm $temp_script

    #-------------------------------------------------#
             # SOUPORCELL DOUBLET DETECTION
    #-------------------------------------------------#
    elif [[ $script_type == "souporcell" ]]; then
        echo "Submitting job for $sample_name (type: $script_type)"
        
        n_genotypes="$5"
        vcf_files="$6"
    
        # Create a temporary LSF script file for ATAC
        temp_script=$(mktemp)
        cat > $temp_script <<EOF
#!/bin/bash
#BSUB -J "souporcell_$sample_name"
#BSUB -q gpuqueue
#BSUB -W 24:00
#BSUB -n 2
#BSUB -R "rusage[mem=50]"
#BSUB -e "$output_directory/$sample_name/%J.err"
#BSUB -o "$output_directory/$sample_name/%J.out"

source ~/.bashrc

module load singularity/3.7.1
cd /data/niecr/ramaiar/multiome_lucyrepo/analysis/souporcell

# ilc pre
singularity exec souporcell_latest.sif souporcell_pipeline.py \
    -i "$cr_outs"/gex_possorted_bam.bam \
    -b "$cr_outs"/filtered_feature_bc_matrix/barcodes.tsv.gz \
    -f "$ref_genome" \
    -t 8 \
    -o "$output_directory"/"$sample_name"/RNAsoup_out \
    -k "$n_genotypes"
    
singularity exec Demuxafy.sif bash souporcell_summary.sh "$output_directory"/"$sample_name"/RNAsoup_out/clusters.tsv > "$output_directory"/"$sample_name"/RNAsoup_out/souporcell_summary.tsv

singularity exec Demuxafy.sif Assign_Indiv_by_Geno.R \
    -r "$vcf_files" \
    -c "$output_directory"/"$sample_name"/RNAsoup_out/cluster_genotypes.vcf \
    -o "$output_directory"/"$sample_name"/RNAsoup_out
EOF

        # Submit the job using the temporary script file
        bsub < $temp_script

        # Remove the temporary script file
        rm $temp_script
    
    fi

}

# Read all sample names in the input directory
sample_names=($(ls "$input_directory"))

# Process samples based on the command-line arguments
if [[ $all_samples == true ]]; then
    for sample_name in "${sample_names[@]}"; do
        if [[ $rna == true ]]; then
            submit_job "$sample_name" "RNA" "$input_directory" "$output_directory" "$overwrite_rna"
        fi
        if [[ $atac == true ]]; then
            submit_job "$sample_name" "ATAC" "$input_directory" "$output_directory"
        fi
    done
elif [[ -n $include_samples ]]; then
    include_samples_array=($(echo "$include_samples" | tr ',' '\n'))
    for sample_name in "${include_samples_array[@]}"; do
        if [[ " ${sample_names[@]} " =~ " $sample_name " ]]; then
            if [[ $rna == true ]]; then
                submit_job "$sample_name" "RNA" "$input_directory" "$output_directory" "$overwrite_rna"
            fi
            if [[ $atac == true ]]; then
                submit_job "$sample_name" "ATAC" "$input_directory" "$output_directory"
            fi
        else
            echo "Sample $sample_name not found in the input directory."
        fi
    done
elif [[ -n $skip_samples ]]; then
    skip_samples_array=($(echo "$skip_samples" | tr ',' '\n'))
    for sample_name in "${sample_names[@]}"; do
        if [[ ! " ${skip_samples_array[@]} " =~ " $sample_name " ]]; then
            if [[ $rna == true ]]; then
                submit_job "$sample_name" "RNA" "$input_directory" "$output_directory" "$overwrite_rna"
            fi
            if [[ $atac == true ]]; then
                submit_job "$sample_name" "ATAC" "$input_directory" "$output_directory"
            fi
        fi
    done
fi

if [[ $souporcell == true ]]; then
    # Path to your CSV file
    csv_file=$dd_file

    # Read the CSV file line by line
    while IFS=, read -r sample_name n_genotypes vcf_files; do
        # Skip the header line
        if [[ $sample_name == "sample_name" ]]; then
            continue
        fi
    
        echo "Processing sample: $sample_name"
        echo "Number of genotypes: $n_genotypes"
        echo "VCF files: $vcf_files"
    
        submit_job "$sample_name" "souporcell" "$input_directory" "$output_directory" "$n_genotypes" "$vcf_files" "$ref_genome"
   
    done < "$csv_file"
fi
