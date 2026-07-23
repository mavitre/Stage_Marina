#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH -A dedicated-cpu@cirad-normal
#SBATCH -J orthofinder
#SBATCH --cpus-per-task 8
#SBATCH --mem-per-cpu 20G
#SBATCH --time=10:00:00
#SBATCH -o orthofinder_%j.out
#SBATCH -e orthofinder_%j.err


INPUT_DIR= $1  # dossier contenant les .faa / .fa
OUTPUT_DIR= $2
THREADS=$SLURM_CPUS_PER_TASK



echo "=== Début : $(date) ==="
echo "Input  : $INPUT_DIR"
echo "Output : $OUTPUT_DIR"
echo "Threads: $THREADS"
 
pixi run orthofinder -f "$INPUT_DIR" -o "$OUTPUT_DIR" -t "$THREADS" -a "$THREADS"
 
echo "=== Fin : $(date) ==="