#!/bin/bash
#SBATCH --nodes=1
#SBATCH --clusters=htc
#SBATCH --ntasks-per-node=1
#SBATCH --time=8-10:00:00 
#SBATCH --job-name=relatedness
#SBATCH --partition=long

#########################################################################################################

#Load samtools and bcftools module
ml PLINK/2.00a2.3_x86_64
ml BCFtools/1.14-GCC-11.2.0

# Set path to genotype_calling path
PATH_GT=/data/biol-gt-genomics/sjoh4959/wm-gwas/data/genotype_calling/

cd $PATH_GT

# KING-robust (plink2)
plink2 --bfile wm-gwas.pruned \
    --make-king-table \
    --not-chr chrZ \
    --allow-extra-chr \
    --out wm-gwas.king

# IBD (PI_HAT) 
/data/biol-gt-genomics/sjoh4959/0.0_wytham_great_tit/src/plink --bfile wm-gwas.pruned \
    --genome \
    --not-chr chrZ \
    --allow-extra-chr \
    --out wm-gwas.genome

# PCA
plink2 --bfile wm-gwas.pruned \
    --pca \
    --not-chr chrZ \
    --allow-extra-chr \
    --out wm-gwas.pca 

plink2 --bfile wm-gwas.pruned  \
     --make-rel square   \
     --allow-extra-chr  \   
     --not-chr chrZ  \   
     --out wm-gwas.pruned