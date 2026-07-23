#!/bin/bash
#SBATCH -A dedicated-gpu@cirad
#SBATCH -J 2helixer
#SBATCH --cpus-per-task 8
#SBATCH --time=48:00:00
#SBATCH -o helixer_%j.out
#SBATCH -e helixer_%j.err
#SBATCH --nodes=1                                    # Généralement 1 noeud pour commencer
#SBATCH --ntasks=1
#SBATCH --mem=128G                                    # La RAM CPU est aussi importante
#SBATCH --partition=gpu-dedicated                    # ADAPTER la partition !
#SBATCH --gres gpu:nvidia_h100_80gb_hbm3_1g.10gb:1

echo "=== Début : $(date) ==="

module load bioinfo-ifb
module load helixer/0.3.3

FASTA=$1
OUTPUT=$2

singularity run helixer-docker_latest.sif Helixer.py --fasta-path $FASTA --lineage land_plant --gff-output-path $OUTPUT 



echo "=== Fin : $(date) ==="