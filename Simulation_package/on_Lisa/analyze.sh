#!/bin/bash
#SBATCH -t 00:10:00 -N 1
#SBATCH --job-name=an

module load 2021
module load Julia/1.6.1-linux-x86_64

#module load mpicopy/4.2-gompi-2020a   # define mpicopy environment
#module load OpenMPI/4.0.3-GCC-9.3.0  # define your environment

mkdir "$TMPDIR"/function_environment
mkdir -p "$TMPDIR"/text_files


cp -r "$HOME"/JULIA/Leaves/all_functions_5d.jl "$TMPDIR"/function_environment  #upload all functions to function environment
cp -r "$HOME"/JULIA/Leaves/single_run_analysis_functions.jl "$TMPDIR"/function_environment  #upload all functions to function environment
cp -r "$HOME"/JULIA/Leaves/automated_comparison_functions.jl "$TMPDIR"/function_environment  #upload all functions to function environment
cp -r "$HOME"/JULIA/Leaves/analyze_and_store_metrics.jl "$TMPDIR"/function_environment  #upload all functions to function environment



#for no_leaves in `seq 0 1`; do

#touch filename.txt
julia "$HOME"/JULIA/Leaves/analyze_and_store_metrics.jl  "$TMPDIR"/function_environment/ "$HOME"/JULIA/Leaves/outputs/output_dir_par/STORED_b_3.0_beta_0.95_leaves_original_trim_population_states.jld2 "$TMPDIR"/text_files/ "LEAVES_"
julia "$HOME"/JULIA/Leaves/analyze_and_store_metrics.jl  "$TMPDIR"/function_environment/ "$HOME"/JULIA/Leaves/outputs/output_dir_par/STORED_b_3.0_beta_0.95_no_leaves_original_trim_population_states.jld2 "$TMPDIR"/text_files/ "NO_LEAVES_"


#done
wait
cp -r -u  "$TMPDIR"/text_files "$HOME"/JULIA/Leaves/outputs
echo "copying done"
