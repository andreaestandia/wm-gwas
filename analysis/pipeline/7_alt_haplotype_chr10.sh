#!/bin/bash
# =============================================================================
# Alt-haplotype vs primary genome comparison – chr10 target region
#
# Target region : chr10:9064949-9205154 (potential structural variant / supergene)
# Primary ref   : GCA_932527175.1_ilOpeBrum1.1_genomic_renamed.fna
# Alt haplotype : GCA_932527245.1_ilOpeBrum1.1_alternate_haplotype_genomic.fna
# Bombyx db     : GCF_000151625.1_Bmor_1.0_protein.faa  (for gene annotation)
#
# Strategy:
#   1. Extract target region from primary reference
#   2. Align entire alt haplotype assembly to primary with minimap2 (asm5)
#   3. Parse PAF to find alt-haplotype contigs overlapping the target region
#   4. Extract those alt-haplotype contig sequences
#   5. Build Bombyx mori protein BLAST DB (downloaded if absent)
#   6. blastx: primary target + alt contig(s) vs Bombyx proteins
# =============================================================================

set -euo pipefail

# ── Paths ────────────────────────────────────────────────────────────────────
BASE=/Volumes/LaCie/wm-gwas
REF=${BASE}/data/ref_genome/GCA_932527175.1/GCA_932527175.1_ilOpeBrum1.1_genomic_renamed.fna
ALT=${BASE}/data/ref_genome/GCA_932527245.1/GCA_932527245.1_ilOpeBrum1.1_alternate_haplotype_genomic.fna
BOMBYX_PROT=${BASE}/data/ref_genome/bombyx_mori/GCF_000151625.1_ASM15162v1_protein.faa
BOMBYX_DB=${BASE}/data/ref_genome/bombyx_mori/bombyx_protein_db

MINIMAP2=${BASE}/src/minimap2/minimap2
SEQTK=${BASE}/src/seqtk/seqtk

OUTDIR=${BASE}/reports/structural_variant/chr10_9064949_9205154
mkdir -p "${OUTDIR}"
mkdir -p "${BASE}/data/ref_genome/bombyx_mori"

# Target region
CHR=chr10
START=9064949
END=9205154
REGION="${CHR}:${START}-${END}"
REGION_LABEL="${CHR}_${START}_${END}"

echo "=== chr10 alt-haplotype analysis ==="
echo "Target: ${REGION}"
echo "Output: ${OUTDIR}"
echo ""

# ── Step 1: Extract target region from primary reference ─────────────────────
echo "[1/6] Extracting ${REGION} from primary reference …"
PRIMARY_FA=${OUTDIR}/primary_${REGION_LABEL}.fa

${SEQTK} subseq <(echo "${REGION}") "${REF}" > "${PRIMARY_FA}" 2>/dev/null || true

# seqtk subseq expects a name-list for whole sequences; use faidx-style region
# Fall back to awk extraction (no samtools required)
python3 - "${REF}" "${CHR}" "${START}" "${END}" "${PRIMARY_FA}" <<'PYEOF'
import sys
ref, chrom, start, end, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
in_target = False
seq_parts = []
with open(ref) as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith('>'):
            in_target = (line[1:].split()[0] == chrom)
            continue
        if in_target:
            seq_parts.append(line)
full_seq = ''.join(seq_parts)
region_seq = full_seq[start-1:end]   # convert to 0-based
with open(out, 'w') as fh:
    fh.write(f'>{chrom}:{start}-{end}\n')
    # wrap at 60 chars
    for i in range(0, len(region_seq), 60):
        fh.write(region_seq[i:i+60] + '\n')
print(f"  Primary region length: {len(region_seq):,} bp")
PYEOF

# ── Step 2: Align alt haplotype to primary reference (asm5 = same-species) ──
echo "[2/6] Aligning alt haplotype to primary reference with minimap2 (asm5) …"
PAF=${OUTDIR}/alt_vs_primary.paf

if [[ ! -s "${PAF}" ]]; then
    ${MINIMAP2} \
        -x asm5 \
        --secondary=no \
        -t 4 \
        "${REF}" \
        "${ALT}" \
        > "${PAF}"
    echo "  PAF written: ${PAF}"
else
    echo "  PAF already exists, skipping alignment."
fi

echo "  Total alignments: $(wc -l < "${PAF}")"

# ── Step 3: Parse PAF – find alt contigs overlapping target region ───────────
echo "[3/6] Parsing PAF for contigs overlapping ${REGION} …"

OVERLAP_TSV=${OUTDIR}/alt_contigs_overlapping_target.tsv
CONTIG_LIST=${OUTDIR}/alt_contig_names.txt

python3 - "${PAF}" "${CHR}" "${START}" "${END}" "${OVERLAP_TSV}" "${CONTIG_LIST}" <<'PYEOF'
import sys
paf, chrom, start, end, out_tsv, out_names = \
    sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6]

hits = []
seen_contigs = set()
with open(paf) as fh:
    for line in fh:
        cols = line.rstrip().split('\t')
        if len(cols) < 12:
            continue
        qname    = cols[0]
        qlen     = int(cols[1])
        qstart   = int(cols[2])
        qend     = int(cols[3])
        strand   = cols[4]
        tname    = cols[5]
        tstart   = int(cols[7])
        tend     = int(cols[8])
        mapq     = int(cols[11])
        # Overlap with target window
        if tname == chrom and tstart < end and tend > start:
            overlap_start = max(tstart, start)
            overlap_end   = min(tend, end)
            overlap_len   = overlap_end - overlap_start
            hits.append((qname, qlen, qstart, qend, strand, tstart, tend, mapq, overlap_len))
            seen_contigs.add(qname)

print(f"  Found {len(hits)} PAF records overlapping target")
print(f"  Unique alt contigs: {len(seen_contigs)}")

header = "alt_contig\talt_len\talt_start\talt_end\tstrand\tprimary_start\tprimary_end\tmapq\toverlap_bp"
with open(out_tsv, 'w') as fh:
    fh.write(header + '\n')
    for h in sorted(hits, key=lambda x: -x[8]):
        fh.write('\t'.join(str(x) for x in h) + '\n')

with open(out_names, 'w') as fh:
    for c in seen_contigs:
        fh.write(c + '\n')

# Print summary table
print(f"\n  {'Contig':<30} {'Len':>10} {'Strand':>6} {'Primary aln':>25} {'Overlap bp':>10} {'MAPQ':>5}")
print("  " + "-"*90)
for h in sorted(hits, key=lambda x: -x[8]):
    qname, qlen, qstart, qend, strand, tstart, tend, mapq, ov = h
    print(f"  {qname:<30} {qlen:>10,} {strand:>6} {tstart:>12,}–{tend:>12,} {ov:>10,} {mapq:>5}")
PYEOF

# ── Step 4: Extract alt haplotype contig sequences ───────────────────────────
echo ""
echo "[4/6] Extracting alt haplotype contig sequences …"
ALT_FA=${OUTDIR}/alt_contigs_${REGION_LABEL}.fa

if [[ ! -s "${CONTIG_LIST}" ]]; then
    echo "  WARNING: No overlapping contigs found — alt haplotype may not cover this region."
else
    python3 - "${ALT}" "${CONTIG_LIST}" "${ALT_FA}" <<'PYEOF'
import sys
genome_fa, name_list, out_fa = sys.argv[1], sys.argv[2], sys.argv[3]
with open(name_list) as fh:
    keep = set(line.strip() for line in fh if line.strip())
print(f"  Extracting {len(keep)} contig(s): {', '.join(sorted(keep))}")
writing = False
count = 0
with open(genome_fa) as fh, open(out_fa, 'w') as oh:
    for line in fh:
        if line.startswith('>'):
            name = line[1:].split()[0]
            writing = name in keep
            if writing:
                count += 1
        if writing:
            oh.write(line)
print(f"  Written {count} sequence(s) to {out_fa}")
PYEOF
fi

# ── Step 5: Download Bombyx mori proteins ─────────────────────────
echo ""
echo "[5/6] Checking Bombyx mori protein database …"
BOMBYX_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/151/625/GCF_000151625.1_ASM15162v1/GCF_000151625.1_ASM15162v1_protein.faa.gz"

if [[ ! -s "${BOMBYX_PROT}" ]]; then
    echo "  Downloading Bombyx mori proteins from NCBI …"
    curl -L --retry 3 -o "${BOMBYX_PROT}.gz" "${BOMBYX_URL}"
    gunzip "${BOMBYX_PROT}.gz"
    echo "  Download complete: ${BOMBYX_PROT}"
else
    echo "  Bombyx mori proteins already present."
fi

if [[ ! -f "${BOMBYX_DB}.phr" ]]; then
    echo "  Building BLAST protein database …"
    makeblastdb \
        -in "${BOMBYX_PROT}" \
        -dbtype prot \
        -out "${BOMBYX_DB}" \
        -title "Bombyx_mori_proteins"
    echo "  BLAST DB created: ${BOMBYX_DB}"
else
    echo "  BLAST DB already exists."
fi

# ── Step 6: blastx – primary target region vs Bombyx proteins ────────────────
echo ""
echo "[6/6] Running blastx: primary target region + alt contig(s) vs Bombyx mori …"

BLAST_PRIMARY=${OUTDIR}/blastx_primary_${REGION_LABEL}.tsv
BLAST_ALT=${OUTDIR}/blastx_alt_contigs_${REGION_LABEL}.tsv

BLAST_FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"

echo "  Primary region …"
blastx \
    -query "${PRIMARY_FA}" \
    -db "${BOMBYX_DB}" \
    -out "${BLAST_PRIMARY}" \
    -outfmt "${BLAST_FMT}" \
    -evalue 1e-5 \
    -num_threads 4 \
    -max_target_seqs 20

echo "  Primary BLAST hits: $(wc -l < "${BLAST_PRIMARY}")"

if [[ -s "${ALT_FA}" ]]; then
    echo "  Alt haplotype contigs …"
    blastx \
        -query "${ALT_FA}" \
        -db "${BOMBYX_DB}" \
        -out "${BLAST_ALT}" \
        -outfmt "${BLAST_FMT}" \
        -evalue 1e-5 \
        -num_threads 4 \
        -max_target_seqs 20

    echo "  Alt BLAST hits: $(wc -l < "${BLAST_ALT}")"
fi

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "=== Analysis complete ==="
echo "Output directory : ${OUTDIR}"
echo ""
echo "Key files:"
echo "  Primary region FASTA : $(basename "${PRIMARY_FA}")"
echo "  Alt-vs-primary PAF   : $(basename "${PAF}")"
echo "  Overlapping contigs  : $(basename "${OVERLAP_TSV}")"
echo "  Alt contig FASTA     : $(basename "${ALT_FA}")"
echo "  BLAST (primary)      : $(basename "${BLAST_PRIMARY}")"
[[ -s "${ALT_FA}" ]] && echo "  BLAST (alt)          : $(basename "${BLAST_ALT}")"
echo ""
echo "Top Bombyx mori hits in primary target region:"
if [[ -s "${BLAST_PRIMARY}" ]]; then
    echo ""
    printf "  %-20s %-12s %-8s %-8s %-12s %s\n" "Hit" "Identity%" "AlnLen" "E-value" "Bitscore" "Description"
    printf "  %s\n" "$(python3 -c "print('-'*100)")"
    awk 'BEGIN{OFS="\t"} {printf "  %-20s %-12s %-8s %-8s %-12s %s\n", $2, $3, $4, $11, $12, $13}' \
        "${BLAST_PRIMARY}" | sort -t$'\t' -k12 -rn | head -20
fi
