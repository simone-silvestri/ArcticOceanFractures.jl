#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -p pi_raffaele # sched_mit_raffaele_gpu # 
#SBATCH --gres=gpu:1
#SBATCH --time=10:00:00
#SBATCH --mem=50GB

module load nvhpc
folder="EVP-rheology"

julia --project=$folder --check-bounds=no $folder/evp_simulation.jl 
