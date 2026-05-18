#!/bin/bash
# =============================================================================
# SV depth-based genotyping – chr10:9,102,316-9,108,300 structural variant
# SLURM array job – one task per sample
#
# Strategy (reads mapped to primary GCA_932527175.1):
#   The ~5.9 kb primary-specific gap (chr10:9,102,316-9,108,300) is absent from
#   the alt haplotype.  Reads originating from the alt allele cannot map there,
#   so depth in this window is proportional to the number of primary allele
#   copies in each individual:
#       depth_ratio ≈ 1.0  →  primary/primary  (0|0)
#       depth_ratio ≈ 0.5  →  heterozygous      (0|1)
#       depth_ratio ≈ 0.0  →  alt/alt           (1|1)
#
# Split reads at the two breakpoints are also counted as a complementary metric.
#
# Submit:
#   N=$(wc -l < /data/biol-gt-genomics/sjoh4959/wm-gwas/data/raw/list_dir)
#   sbatch --array=1-${N} 8_sv_depth_genotype.sh
#
# Merge per-sample results afterwards:
#   head -1 ${OUTDIR}/sample_1.tsv > sv_depth.tsv
#   tail -n +2 -q ${OUTDIR}/*.tsv  >> sv_depth.tsv
#
# Sample list format (two columns, no header):
#   SAMPLE_NAME   /path/to/bam/directory
#
# Requires: samtools >= 1.10
# =============================================================================

#SBATCH --job-name=sv_depth
#SBATCH --output=/data/biol-gt-genomics/sjoh4959/wm-gwas/logs/sv_depth_%A_%a.out
#SBATCH --error=/data/biol-gt-genomics/sjoh4959/wm-gwas/logs/sv_depth_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --array=1-501:1
#SBATCH --cpus-per-task=1

set -euo pipefail

ml SAMtools/1.14-GCC-11.2.0

# ── Sample list and output directory ─────────────────────────────────────────
Sample_List=/data/biol-gt-genomics/sjoh4959/wm-gwas/data/raw/list_dir
OUTDIR=/data/biol-gt-genomics/sjoh4959/wm-gwas/reports/structural_variant/chr10_9064949_9205154/sv_depth_per_sample
mkdir -p "${OUTDIR}"
mkdir -p "$(dirname "${OUTDIR}")"/../../../logs

# ── STEP 2: resolve sample name and BAM path from array task ID ──────────────
SAMPLE_NAME=$(cat "${Sample_List}" | head -n "${SLURM_ARRAY_TASK_ID}" | tail -1 | awk '{print $1}')
SAMPLE_DIRECTORY=$(cat "${Sample_List}" | head -n "${SLURM_ARRAY_TASK_ID}" | tail -1 | awk '{print $2}')

BAM="${SAMPLE_DIRECTORY}/${SAMPLE_NAME}.bam"

echo "Array task : ${SLURM_ARRAY_TASK_ID}"
echo "Sample     : ${SAMPLE_NAME}"
echo "BAM        : ${BAM}"
echo ""

if [[ ! -f "${BAM}" ]]; then
    echo "ERROR: BAM not found: ${BAM}" >&2
    exit 1
fi

if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
    echo "ERROR: No index found for ${BAM}" >&2
    exit 1
fi

# ── Region definitions (primary chr10 coordinates) ───────────────────────────
CHR="chr10"

# The SV inner region – absent from alt haplotype
SV_INNER="${CHR}:9102316-9108300"    # ~5.9 kb

# Flanking anchors – conserved in both haplotypes; used to normalise depth
ANCHOR_L="${CHR}:9080000-9100000"    # 20 kb left of breakpoint 1
ANCHOR_R="${CHR}:9110000-9130000"    # 20 kb right of breakpoint 2

# Breakpoints for split-read counting (±50 bp windows)
BP1="${CHR}:9102266-9102366"
BP2="${CHR}:9108250-9108350"

# ── Helper: mean depth over a region ─────────────────────────────────────────
mean_depth() {
    samtools depth -a -r "${2}" "${1}" 2>/dev/null \
        | awk '{sum+=$3; n++} END {if(n>0) printf "%.4f", sum/n; else print "NA"}'
}

# ── Helper: count soft-clipped reads near a breakpoint ───────────────────────
count_clipped() {
    samtools view -F 2308 "${1}" "${2}" 2>/dev/null \
        | awk '$6 ~ /S/ {n++} END {print (n>0)?n:0}'
}

# ── Compute metrics ───────────────────────────────────────────────────────────
D_SV=$(mean_depth "${BAM}" "${SV_INNER}")
D_AL=$(mean_depth "${BAM}" "${ANCHOR_L}")
D_AR=$(mean_depth "${BAM}" "${ANCHOR_R}")

D_ANC=$(awk -v l="${D_AL}" -v r="${D_AR}" 'BEGIN {
    if (l=="NA" || r=="NA") print "NA"
    else printf "%.4f", (l+r)/2
}')

SV_RATIO=$(awk -v sv="${D_SV}" -v anc="${D_ANC}" 'BEGIN {
    if (sv=="NA" || anc=="NA" || anc+0==0) print "NA"
    else printf "%.4f", sv/anc
}')

C_BP1=$(count_clipped "${BAM}" "${BP1}")
C_BP2=$(count_clipped "${BAM}" "${BP2}")

echo "depth_sv_inner : ${D_SV}"
echo "depth_anchor   : ${D_ANC}  (L=${D_AL}, R=${D_AR})"
echo "sv_ratio       : ${SV_RATIO}"
echo "clipped BP1/2  : ${C_BP1} / ${C_BP2}"

# ── Write per-sample TSV ──────────────────────────────────────────────────────
OUTFILE="${OUTDIR}/${SAMPLE_NAME}.tsv"

printf "sample\tdepth_sv_inner\tdepth_anchor_L\tdepth_anchor_R\t" >  "${OUTFILE}"
printf "depth_anchor_mean\tsv_ratio\tclipped_BP1\tclipped_BP2\n"  >> "${OUTFILE}"

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${SAMPLE_NAME}" "${D_SV}" "${D_AL}" "${D_AR}" \
    "${D_ANC}" "${SV_RATIO}" "${C_BP1}" "${C_BP2}" \
    >> "${OUTFILE}"

echo ""
echo "Done. Written to: ${OUTFILE}"
