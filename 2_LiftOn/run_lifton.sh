#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH -A dedicated-cpu@cirad-normal
#SBATCH -J lifton_speltoides
#SBATCH -N 1
#SBATCH --cpus-per-task 6
#SBATCH --mem-per-cpu 20G
#SBATCH --time=24:00:00

#SBATCH --output=lifton_%J.out
#SBATCH --error=lifton_%J.err

gff=$1
ref_fasta=$2
target_fasta=$3
output=$4
nb_cpus=6

START=$(date +%s)

singularity exec /storage/replicated/cirad/projects/GE2POP/APPTAINER_IMAGES/lifton/LiftOn.sif lifton -g $gff -o $output $target_fasta $ref_fasta -t $nb_cpus

END=$(date +%s)
echo "Temps : $((END - START)) secondes"
