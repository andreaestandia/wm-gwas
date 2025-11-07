#!/bin/bash
#SBATCH --nodes=1
#SBATCH --clusters=htc
#SBATCH --ntasks-per-node=1
#SBATCH --time=8-10:00:00 
#SBATCH --job-name=ld_pruning
#SBATCH --partition=long

#########################################################################################################

#Load samtools and bcftools module
ml PLINK/2.00a2.3_x86_64
ml BCFtools/1.14-GCC-11.2.0

# Set path to genotype_calling path
PATH_GT=/data/biol-gt-genomics/sjoh4959/wm-gwas/data/genotype_calling/

cd $PATH_GT

# generate pruned SNP list
plink2 --bfile wm-gwas.qc \
    --indep-pairwise 50 5 0.2 \
    --out wm-gwas.prune

# keep pruned SNPs
plink2 --bfile wm-gwas.qc \
    --extract wm-gwas.prune.prune.in \
    --make-bed \
    --out wm-gwas.pruned