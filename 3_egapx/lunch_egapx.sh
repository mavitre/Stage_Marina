#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH -A dedicated-cpu@cirad-normal
#SBATCH -J egapx_s
#SBATCH -N 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 128G
#SBATCH --time=24:00:00

#SBATCH --output=egapx_searsii.out
#SBATCH --error=egapx_searsii.err


START=$(date +%s)

INPUT_YAML=$1
OUTPUT_DIR=$2
WORKDIR=$3  #mettre dans le scratch


# Lancer EGAPx
pixi run python3 egapx/ui/egapx.py $INPUT_YAML -e slurm -o $OUTPUT_DIR -w $WORKDIR 

END=$(date +%s)
echo "Temps : $((END - START)) secondes"
