#!/bin/bash
#SBATCH --nodes=1
#SBATCH --clusters=htc
#SBATCH --ntasks-per-node=1
#SBATCH --time=8-10:00:00 
#SBATCH --job-name=ld_pruning
#SBATCH --partition=long

#########################################################################################################

# =========================
# Define paths
# =========================
PHENO_PATH="/data/biol-gt-genomics/sjoh4959/wm-gwas/data/resources/"
QCTOOL_PATH="/data/biol-gt-genomics/sjoh4959/wm-gwas/src/qctool/"
GEMMA_PATH="/data/biol-gt-genomics/sjoh4959/wm-gwas/src/"

# =========================
# Step 1: Filter variants
# =========================
# Filter SNPs by MAF and missingness
bcftools +fill-tags wm-gwas.vcf.gz -- -t MAF,F_MISSING | \
  bcftools view -i 'MAF>=0.05 && F_MISSING<=0.3' -Oz -o wm-gwas.maf05.geno03.vcf.gz

# =========================
# Step 2: Filter samples
# =========================
# Identify samples with >10% missing genotypes
bcftools query -f '%SAMPLE\t[%GT]\n' wm-gwas.maf05.geno03.vcf.gz | \
  awk '{
    miss=0; total=0;
    for(i=2;i<=NF;i++){
      total++;
      if($i=="." || $i=="./.") miss++;
    }
    if(miss/total>0.1) print $1
  }' > samples.to.remove

# Remove those samples from the dataset
bcftools view -S ^samples.to.remove -Oz -o wm-gwas.filtered.vcf.gz wm-gwas.maf05.geno03.vcf.gz

# =========================
# Step 3: Convert to BIMBAM format
# =========================
"${QCTOOL_PATH}qctool" \
  -g wm-gwas.filtered.vcf.gz \
  -ofiletype bimbam_dosage \
  -og wm-gwas.bimbam

# =========================
# Step 4: Compute relatedness matrix
# =========================
"${GEMMA_PATH}gemma" \
  -g wm-gwas.bimbam \
  -gk 1 \
  -p "${PHENO_PATH}pheno_gemma.txt" \
  -o wm-gwas.relatedness 

plink2 --bfile wm-gwas.pruned \
  --make-rel square \
  --out output/relatedness \
  --allow-extra-chr \
  --not-chr chrZ

# =========================
# Step 5: Run LMM analysis
# =========================
"${GEMMA_PATH}gemma" \
  -g wm-gwas.bimbam \
  -k output/wm-gwas.relatedness.cXX.txt \
  -lmm 4 \
  -p "${PHENO_PATH}pheno_gemma.txt" \
  -o wm-gwas.gemma

"${GEMMA_PATH}gemma" \
  -g wm-gwas.bimbam \
  -k output/wm-gwas.relatedness.cXX.txt \
  -bslmm 1 \
  -p "${PHENO_PATH}pheno_gemma.txt" \
  -o wm-gwas_bslmm.gemma