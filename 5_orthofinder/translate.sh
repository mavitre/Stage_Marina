#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH -A dedicated-cpu@cirad-normal
#SBATCH --mem 8G
#SBATCH --time=1:30:00

module load bioinfo-ifb
module load agat/1.4.2

GFF=$1
FASTA=$2
OUTPUT=$3

agat_sp_extract_sequences.pl -g GFF -f FASTA -t cds -o OUTPUT # pour cds seulement

#agat_sp_extract_sequences.pl -g GFF -f FASTA -p -o OUTPUT # pour protéines

