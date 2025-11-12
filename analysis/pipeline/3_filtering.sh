#!/bin/bash
#SBATCH --nodes=1
#SBATCH --clusters=htc
#SBATCH --ntasks-per-node=1
#SBATCH --time=8-10:00:00 
#SBATCH --job-name=filtering
#SBATCH --partition=long

#########################################################################################################

#Load plink and bcftools module
ml PLINK/2.00a2.3_x86_64
ml BCFtools/1.14-GCC-11.2.0

# Set path to reference assembly and list of bam files (bam.list)
# Note: bam files need to be indexed (using samtools index) 
REF=/data/biol-wm/sjoh4959/0.0_winter-moth/data/ref_genome/GCA_932527175.1/GCA_932527175.1_ilOpeBrum1.1_genomic_renamed.fna
BAMs=/data/biol-gt-genomics/sjoh4959/wm-gwas/data/resources/list_bam
PATH_OUT=/data/biol-gt-genomics/sjoh4959/wm-gwas/data/genotype_calling/

cd $PATH_OUT

# Keep biallelic SNPs only and apply site-level filters between 2x and 24x average depth and at least QUAL 10
bcftools view -m2 -M2 -v snps -Ou wm-gwas.bcf | bcftools view -i 'QUAL>10 && INFO/DP>=1000 && INFO/DP<=12000' -Ob -o wm.sitefiltered.bcf

# convert to vcf.gz
bcftools view wm.sitefiltered.bcf -Oz -o wm-gwas_tmp.vcf.gz

# add unique IDs in the format CHROM_POS so when doing LD pruning we can track positions
bcftools annotate -x ID -I +'%CHROM\_%POS' wm-gwas_tmp.vcf.gz -Oz -o wm-gwas.vcf.gz

# convert to plink2 pgen (or make-bed)
plink2 --vcf wm-gwas.vcf.gz \
    --make-bed \
    --out wm-gwas \
    --allow-extra-chr

# apply MAF & missingness filters. MAF 0.05 because sample size is small
plink2 --allow-extra-chr \
       --bfile wm-gwas \
       --geno 0.3 \
       --maf 0.05 \
       --make-bed \
       --out wm-gwas.qc
