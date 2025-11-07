#!/bin/bash
#SBATCH --nodes=1
#SBATCH --clusters=htc
#SBATCH --ntasks-per-node=1
#SBATCH --time=10-10:00:00 
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
PATH_OUT=/data/biol-gt-genomics/sjoh4959/wm-gwas/data/genotype_calling/GLs/

#snp calling with genotype likelihoods - keep PL (phred-scaled likelihoods) and GP (genotype probabilities)
bcftools mpileup -Ou -f $REF -b $BAMs \
  --min-MQ 10 \
  --min-BQ 10 \
  --annotate FORMAT/DP,FORMAT/AD,INFO/AD | \
bcftools call -mv \
  -f GQ,GP \
  -Ou | \
bcftools +fill-tags -Ob -o "${PATH_OUT}wm-gwas.bcf" -- -t AF,MAF

# Keep biallelic SNPs only and apply site-level filters
# Adjust depth thresholds based on sample size (N_samples * 2 to N_samples * 10)
# For 500 samples at 4.5x avg: 1000 to 5000x total depth is reasonable
bcftools view -m2 -M2 -v snps -Ou "${PATH_OUT}wm-gwas.bcf" | \
  bcftools view -i 'QUAL>20 && INFO/DP>=1000 && INFO/DP<=12000' \
  -Ou -o "${PATH_OUT}wm.sitefiltered.bcf"