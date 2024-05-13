#!/bin/bash
#SBATCH -t 00:10:00 -N 1
#SBATCH --job-name=n00n

module load 2021
module load Julia/1.6.1-linux-x86_64

#module load mpicopy/4.2-gompi-2020a   # define mpicopy environment
#module load OpenMPI/4.0.3-GCC-9.3.0  # define your environment

mkdir "$TMPDIR"/function_environment
mkdir -p "$TMPDIR"/output_dir_par


cp -r "$HOME"/JULIA/Leaves/all_functions_5d.jl "$TMPDIR"/function_environment  #upload all functions to function environment
cp -r "$HOME"/JULIA/Leaves/iterators_and_parameters.jl "$TMPDIR"/function_environment  #upload all functions to function environment


### arguments are: bi, betai, no_leaves, trimming, function_environment, output_...
### bi=5 --> b=3
### betai =5 --? beta = 0.95
for no_leaves in `seq 0 1`; do
                                julia "$HOME"/JULIA/Leaves/one_run_and_save_full.jl  5 5 $no_leaves  -1  "$TMPDIR"/function_environment/ "$TMPDIR"/output_dir_par/ &
done
wait
cp -r -u  "$TMPDIR"/output_dir_par "$HOME"/JULIA/Leaves/outputs
echo "copying done"
