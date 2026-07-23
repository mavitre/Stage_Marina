#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH -A dedicated-cpu@cirad-normal
#SBATCH -J remove
#SBATCH -N 1
#SBATCH --mem 16G
#SBATCH --time=1:00:00


#python3 standardize.py /storage/replicated/cirad/projects/GE2POP/REFERENCES/TRITICUM/CHINESE_SPRING/ANNOTATION/annotation_ncbi/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.gff /scratch/users/vitrem/annot_2N_WWR/1_parse_ncbi_gff/test_standardize/Ae_tauschii_standardize2.gff -c /scratch/users/vitrem/annot_2N_WWR/1_parse_ncbi_gff/test_standardize/test_tauschii.txt

python3 standardize.py /storage/replicated/cirad/projects/GE2POP/REFERENCES/TRITICUM/CHINESE_SPRING/ANNOTATION/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3 /scratch/users/vitrem/annot_2N_WWR/1_parse_ncbi_gff/test_standardize/T_aestivum_standardize.gff -c /scratch/users/vitrem/annot_2N_WWR/1_parse_ncbi_gff/test_standardize/aestivum_URGI.txt