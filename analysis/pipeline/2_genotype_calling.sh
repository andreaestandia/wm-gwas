#!/bin/bash
#SBATCH --nodes=1
#SBATCH --clusters=htc
#SBATCH --ntasks-per-node=1
#SBATCH --time=8-10:00:00 
#SBATCH --job-name=gt_call
#SBATCH --partition=long

#########################################################################################################

#Load samtools and bcftools module
ml SAMtools/1.16.1-GCC-11.3.0
ml BCFtools/1.14-GCC-11.2.0

# Set path to reference assembly and list of bam files (bam.list)
# Note: bam files need to be indexed (using samtools index) 
REF=/data/biol-wm/sjoh4959/0.0_winter-moth/data/ref_genome/GCA_932527175.1/GCA_932527175.1_ilOpeBrum1.1_genomic_renamed.fna
BAMs=/data/biol-gt-genomics/sjoh4959/wm-gwas/data/resources/list_bam
PATH_OUT=/data/biol-gt-genomics/sjoh4959/wm-gwas/data/genotype_calling/

cd $PATH_OUT

#snp calling based on a list of paths to the sorted BAM files
bcftools mpileup -Ou -f $REF -b $BAMs | \
bcftools call -mv -Ob -o "${PATH_OUT}wm-gwas.bcf"
