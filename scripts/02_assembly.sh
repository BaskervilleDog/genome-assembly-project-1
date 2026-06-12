#!/usr/bin/env bash
#==============================================================================
# STAGES 2–6 — FILTER, DECONTAMINATE, ASSEMBLE & POLISH
# Genome Assembly Pipeline: Trichoderma harzianum TW11
#==============================================================================
# Author : Gianlucca de Urzêda Alves
# Date   : 11-02-2026
#
# Purpose:
#   - Filter PacBio reads (fastplong) and Illumina reads (fastp)
#   - Decontaminate reads with Kraken2
#   - Assemble with three independent assemblers (Flye, Raven, wtdbg2)
#   - Polish: Racon (×2) → Polypolish → NextPolish (×2)
#
# Tools (mamba environments):
#   fastplong    | fastplong
#   fastp        | fastp
#   seqkit       | seqkit
#   fastqc       | fastqc
#   multiqc      | multiqc
#   kraken2      | kraken2
#   flye292      | flye
#   raven        | raven
#   wtdbg2       | wtdbg2, wtpoa-cns
#   minimap2     | minimap2
#   racon1420    | racon
#   bwa          | bwa
#   polypolish   | polypolish
#   nextpolish39 | nextPolish
#
# Outputs: 02_filtered/, 03_qc_filtered/, 04_contaminants/,
#          05_decontaminated/, 06_assemblies/, 08_polishing/
#==============================================================================

set -euo pipefail

#==============================================================================
# CONFIGURATION
#==============================================================================

THREADS=18
KRAKEN_DB="${HOME}/databases/kraken2_db"
PROJECT_DIR="$(pwd)"

# Directories
DOWNLOAD_DIR="${PROJECT_DIR}/00_downloads"
FILTERED_DIR="${PROJECT_DIR}/02_filtered"
QC_FILTERED_DIR="${PROJECT_DIR}/03_qc_filtered"
CONTAMINANTS_DIR="${PROJECT_DIR}/04_contaminants"
DECON_DIR="${PROJECT_DIR}/05_decontaminated"
ASM_DIR="${PROJECT_DIR}/06_assemblies"
POLISH_DIR="${PROJECT_DIR}/08_polishing/longreads"
LOG_DIR="${PROJECT_DIR}/logs"

#==============================================================================
# DIRECTORY SETUP
#==============================================================================

mkdir -p \
    "${FILTERED_DIR}/pacbio" \
    "${FILTERED_DIR}/illumina" \
    "${QC_FILTERED_DIR}" \
    "${CONTAMINANTS_DIR}" \
    "${DECON_DIR}" \
    "${ASM_DIR}/flye" \
    "${ASM_DIR}/raven" \
    "${ASM_DIR}/wtdbg2" \
    "${POLISH_DIR}" \
    "${LOG_DIR}"

#==============================================================================
# UTILITY FUNCTIONS
#==============================================================================

log_step() {
    echo ""
    echo "========================================================================"
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP: $*"
    echo "========================================================================"
}

log_info() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] [INFO] $*"
}

#==============================================================================
# STEP 1 — FILTER PACBIO READS (fastplong)
#==============================================================================
# fastplong is long-read-aware; standard short-read trimmers (fastp,
# Trimmomatic) are not appropriate for PacBio reads.
#
# Quality filtering is explicitly DISABLED (-Q) because quality scores are
# absent in this SRA dataset (Q=0 placeholders). Only length filtering is
# applied: reads <1 kb are too short to contribute to assembly contiguity
# and may introduce false overlaps in the assembler.

log_step "Filtering PacBio long reads (fastplong)"

mamba run -n fastplong fastplong \
    -w "${THREADS}" \
    -Q \
    --length_required 1000 \
    -i "${DOWNLOAD_DIR}/SRR10848482_1.fastq.gz" \
    -o "${FILTERED_DIR}/pacbio/SRR10848482_1_filtered.fastq.gz" \
    --html "${QC_FILTERED_DIR}/SRR10848482_1.fastplong.html" \
    --json "${QC_FILTERED_DIR}/SRR10848482_1.fastplong.json" \
    > "${LOG_DIR}/SRR10848482_1.fastplong.stdout.log" \
    2> "${LOG_DIR}/SRR10848482_1.fastplong.stderr.log"

log_info "PacBio filtering complete → ${FILTERED_DIR}/pacbio/"

#==============================================================================
# STEP 2 — FILTER ILLUMINA READS (fastp)
#==============================================================================
# fastp combines adapter trimming, quality filtering, and overlap-based error
# correction in a single pass.
#
# --detect_adapter_for_pe : auto-detects PE adapter sequences (no hardcoding)
# --qualified_quality_phred 30 : Q30 minimum base quality
# --length_required 150 : discard reads too short after trimming
# --correction : overlap-based error correction using read-pair overlaps;
#                corrects sequencing errors before short-read polishing

log_step "Filtering Illumina short reads (fastp)"

for srr in SRR10848483 SRR10848484; do
    log_info "Filtering ${srr}..."
    mamba run -n fastp fastp \
        -w "${THREADS}" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 30 \
        --length_required 150 \
        --correction \
        -i "${DOWNLOAD_DIR}/${srr}_1.fastq.gz" \
        -I "${DOWNLOAD_DIR}/${srr}_2.fastq.gz" \
        -o "${FILTERED_DIR}/illumina/${srr}_1_filtered.fastq.gz" \
        -O "${FILTERED_DIR}/illumina/${srr}_2_filtered.fastq.gz" \
        --html "${QC_FILTERED_DIR}/${srr}.fastp.html" \
        --json "${QC_FILTERED_DIR}/${srr}.fastp.json" \
        > "${LOG_DIR}/${srr}.fastp.stdout.log" \
        2> "${LOG_DIR}/${srr}.fastp.stderr.log"
done

log_info "Illumina filtering complete → ${FILTERED_DIR}/illumina/"

#==============================================================================
# STEP 3 — QC: FILTERED READS (SeqKit + FastQC + MultiQC)
#==============================================================================

log_step "QC: Filtered reads — SeqKit stats"

# PacBio filtered reads
mamba run -n seqkit seqkit \
    stats -a \
    "${FILTERED_DIR}/pacbio/SRR10848482_1_filtered.fastq.gz" \
    > "${QC_FILTERED_DIR}/SRR10848482_1_filtered_stats.txt" \
    2>> "${LOG_DIR}/seqkit_filtered.stderr.log"

# Illumina filtered reads
for srr in SRR10848483 SRR10848484; do
    for read in 1 2; do
        mamba run -n seqkit seqkit \
            stats -a \
            "${FILTERED_DIR}/illumina/${srr}_${read}_filtered.fastq.gz" \
            > "${QC_FILTERED_DIR}/${srr}_${read}_filtered_stats.txt" \
            2>> "${LOG_DIR}/seqkit_filtered.stderr.log"
    done
done

log_step "QC: Filtered Illumina reads — FastQC + MultiQC"

mamba run -n fastqc fastqc \
    -t "${THREADS}" \
    -o "${QC_FILTERED_DIR}/" \
    "${FILTERED_DIR}/illumina/SRR10848483_1_filtered.fastq.gz" \
    "${FILTERED_DIR}/illumina/SRR10848483_2_filtered.fastq.gz" \
    "${FILTERED_DIR}/illumina/SRR10848484_1_filtered.fastq.gz" \
    "${FILTERED_DIR}/illumina/SRR10848484_2_filtered.fastq.gz" \
    > "${LOG_DIR}/fastqc_filtered.stdout.log" \
    2> "${LOG_DIR}/fastqc_filtered.stderr.log"

mamba run -n multiqc multiqc \
    "${QC_FILTERED_DIR}/" \
    -o "${QC_FILTERED_DIR}/multiqc_report" \
    > "${LOG_DIR}/multiqc_filtered.stdout.log" \
    2> "${LOG_DIR}/multiqc_filtered.stderr.log"

log_info "Filtered QC reports written to ${QC_FILTERED_DIR}/"

#==============================================================================
# STEP 4 — DECONTAMINATION (Kraken2)
#==============================================================================
# Fungal cultures (especially soil-derived T. harzianum) frequently carry
# bacterial and human DNA from culture media, researcher handling, or reagents.
# Assembling without decontamination causes chimeric contigs and inflated
# genome size estimates.
#
# Kraken2 uses exact k-mer matching against a reference database (standard
# 16 Gb build: bacteria, archaea, viruses, human).
#
# Reads marked UNCLASSIFIED (no database hit) are retained as the clean
# dataset. Classified (contaminant) reads are saved separately for inspection.

log_step "Decontamination: Kraken2 taxonomic classification"

# PacBio long reads (single-end)
log_info "Classifying PacBio reads..."
mamba run -n kraken2 kraken2 \
    --db "${KRAKEN_DB}" \
    --threads "${THREADS}" \
    --report "${CONTAMINANTS_DIR}/pacbio_report.k2report" \
    --unclassified-out "${DECON_DIR}/SRR10848482_unclassified.fastq" \
    --classified-out "${CONTAMINANTS_DIR}/SRR10848482_classified.fastq" \
    "${FILTERED_DIR}/pacbio/SRR10848482_1_filtered.fastq.gz" \
    > "${CONTAMINANTS_DIR}/pacbio_output.kraken2" \
    2> "${LOG_DIR}/kraken2_pacbio.stderr.log"

# Illumina paired-end reads (both runs)
for srr in SRR10848483 SRR10848484; do
    log_info "Classifying ${srr}..."
    mamba run -n kraken2 kraken2 \
        --db "${KRAKEN_DB}" \
        --threads "${THREADS}" \
        --paired \
        --report "${CONTAMINANTS_DIR}/${srr}_report.k2report" \
        --unclassified-out "${DECON_DIR}/${srr}_unclassified#.fastq" \
        --classified-out "${CONTAMINANTS_DIR}/${srr}_classified#.fastq" \
        "${FILTERED_DIR}/illumina/${srr}_1_filtered.fastq.gz" \
        "${FILTERED_DIR}/illumina/${srr}_2_filtered.fastq.gz" \
        > "${CONTAMINANTS_DIR}/${srr}_output.kraken2" \
        2> "${LOG_DIR}/kraken2_${srr}.stderr.log"
done

# Compress all decontaminated FASTQ files
log_info "Compressing decontaminated reads..."
pigz -p "${THREADS}" "${DECON_DIR}"/*.fastq

log_info "Decontaminated reads → ${DECON_DIR}/"
log_info "Contaminant reports → ${CONTAMINANTS_DIR}/"

#==============================================================================
# STEP 5 — GENOME ASSEMBLY (Flye + Raven + wtdbg2)
#==============================================================================
# Three assemblers with fundamentally different graph approaches are run on
# the same decontaminated PacBio reads. BUSCO and QUAST are applied to all
# three outputs (in 03_evaluate.sh) before selecting the best for polishing.
#
# Running multiple assemblers mitigates algorithm-specific biases. Empirical
# evaluation — not theoretical assumptions — determines which performs best
# on this genome and coverage depth.

PACBIO_CLEAN="${DECON_DIR}/SRR10848482_unclassified.fastq.gz"

# — — — — — — — — — — — — — — — — — — — — — — — — — —
# 5a. FLYE (repeat graph)
# Repeat-aware graph explicitly models repetitive sequences during assembly.
# Designed for the ~178× coverage generated here. The --pacbio-corr mode
# adjusts internal error-tolerance for corrected PacBio reads.
# Selected as the best assembly based on post-run BUSCO/QUAST comparison.
# — — — — — — — — — — — — — — — — — — — — — — — — — —

log_step "Assembly: Flye (repeat graph)"

mamba run -n flye292 flye \
    --pacbio-corr "${PACBIO_CLEAN}" \
    --genome-size 40m \
    --out-dir "${ASM_DIR}/flye" \
    --threads "${THREADS}" \
    > "${LOG_DIR}/SRR10848482_flye.stdout.log" \
    2> "${LOG_DIR}/SRR10848482_flye.stderr.log"

log_info "Flye assembly → ${ASM_DIR}/flye/assembly.fasta"

# — — — — — — — — — — — — — — — — — — — — — — — — — —
# 5b. RAVEN (string graph / OLC)
# Fast alternative with internal Racon polishing. Memory-efficient and
# often produces high contiguity on fungal genomes, providing a strong
# benchmark for comparison.
# — — — — — — — — — — — — — — — — — — — — — — — — — —

log_step "Assembly: Raven (string graph)"

mamba run -n raven raven \
    --threads "${THREADS}" \
    "${PACBIO_CLEAN}" \
    > "${ASM_DIR}/raven/raven_assembly.fasta" \
    2> "${LOG_DIR}/SRR10848482_raven.stderr.log"

log_info "Raven assembly → ${ASM_DIR}/raven/raven_assembly.fasta"

# — — — — — — — — — — — — — — — — — — — — — — — — — —
# 5c. WTDBG2 (fuzzy de Bruijn graph)
# Structurally distinct from Flye and Raven, increasing the probability
# of identifying the best assembly through empirical evaluation.
# Consensus generation (wtpoa-cns) is run separately after graph construction.
# -x rs : PacBio RS2/Sequel error-tolerance preset
# — — — — — — — — — — — — — — — — — — — — — — — — — —

log_step "Assembly: wtdbg2 (fuzzy de Bruijn graph)"

mamba run -n wtdbg2 wtdbg2 \
    -x rs -g 40m -t "${THREADS}" \
    -i "${PACBIO_CLEAN}" \
    -fo "${ASM_DIR}/wtdbg2/wtdbg2" \
    > "${LOG_DIR}/SRR10848482_wtdbg2.stdout.log" \
    2> "${LOG_DIR}/SRR10848482_wtdbg2.stderr.log"

mamba run -n wtdbg2 wtpoa-cns \
    -t "${THREADS}" \
    -i "${ASM_DIR}/wtdbg2/wtdbg2.ctg.lay.gz" \
    -fo "${ASM_DIR}/wtdbg2/wtdbg2_assembly.fasta" \
    >> "${LOG_DIR}/SRR10848482_wtdbg2.stdout.log" \
    2>> "${LOG_DIR}/SRR10848482_wtdbg2.stderr.log"

log_info "wtdbg2 assembly → ${ASM_DIR}/wtdbg2/wtdbg2_assembly.fasta"

#==============================================================================
# STEP 6 — POLISHING
#==============================================================================
# Polishing strategy (Racon ×2 → Polypolish → NextPolish ×2):
#
#   Long-read polishing first (Racon): corrects structural errors (indels,
#   substitutions) that would otherwise confuse short-read alignments.
#
#   Short-read polishing second (Polypolish then NextPolish): corrects
#   base-level SNPs and indels using high-accuracy Illumina reads.
#
# The selected assembly for polishing is FLYE (best BUSCO/QUAST result).
# Run 07_evaluate_raw.sh before this step to confirm the selection.

FLYE_ASM="${ASM_DIR}/flye/assembly.fasta"

# — — — — — — — — — — — — — — — — — — — — — — — — — —
# 6a. RACON ROUND 1 (minimap2 → Racon)
# minimap2 -x map-pb : preset calibrated for PacBio CLR gap penalties.
# First Racon round corrects the largest indel/substitution errors.
# — — — — — — — — — — — — — — — — — — — — — — — — — —

log_step "Polishing Round 1: minimap2 + Racon (long reads)"

mamba run -n minimap2 minimap2 \
    -t "${THREADS}" -x map-pb \
    "${FLYE_ASM}" \
    "${PACBIO_CLEAN}" \
    > "${POLISH_DIR}/pacbio_flye_alignments_r1.paf" \
    2> "${LOG_DIR}/minimap2_racon_r1.stderr.log"

mamba run -n racon1420 racon \
    -t "${THREADS}" \
    "${PACBIO_CLEAN}" \
    "${POLISH_DIR}/pacbio_flye_alignments_r1.paf" \
    "${FLYE_ASM}" \
    > "${POLISH_DIR}/pacbio_flye_assembly.racon1.fasta" \
    2> "${LOG_DIR}/racon_r1.stderr.log"

log_info "Racon round 1 → ${POLISH_DIR}/pacbio_flye_assembly.racon1.fasta"

# — — — — — — — — — — — — — — — — — — — — — — — — — —
# 6b. RACON ROUND 2 (minimap2 → Racon)
# Second round maps to the round-1 output and catches residual errors
# in regions with poor alignment coverage in round 1.
# A third round yields diminishing returns and is not applied.
# — — — — — — — — — — — — — — — — — — — — — — — — — —

log_step "Polishing Round 2: minimap2 + Racon (long reads)"

mamba run -n minimap2 minimap2 \
    -t "${THREADS}" -x map-pb \
    "${POLISH_DIR}/pacbio_flye_assembly.racon1.fasta" \
    "${PACBIO_CLEAN}" \
    > "${POLISH_DIR}/pacbio_flye_alignments_r2.paf" \
    2> "${LOG_DIR}/minimap2_racon_r2.stderr.log"

mamba run -n racon1420 racon \
    -t "${THREADS}" \
    "${PACBIO_CLEAN}" \
    "${POLISH_DIR}/pacbio_flye_alignments_r2.paf" \
    "${POLISH_DIR}/pacbio_flye_assembly.racon1.fasta" \
    > "${POLISH_DIR}/pacbio_flye_assembly.racon2.fasta" \
    2> "${LOG_DIR}/racon_r2.stderr.log"

log_info "Racon round 2 → ${POLISH_DIR}/pacbio_flye_assembly.racon2.fasta"

# — — — — — — — — — — — — — — — — — — — — — — — — — —
# 6c. POLYPOLISH (BWA-MEM -a → Polypolish)
# Polypolish uses ALL alignments (bwa mem -a) rather than only uniquely
# mapping reads, enabling polishing of repetitive regions that Pilon skips.
# Both Illumina runs (SRR10848483 + SRR10848484) are aligned separately
# and jointly applied for maximum coverage (~60×).
# — — — — — — — — — — — — — — — — — — — — — — — — — —

log_step "Polishing: Polypolish (short reads — all alignments)"

RACON2="${POLISH_DIR}/pacbio_flye_assembly.racon2.fasta"

# Index the assembly for BWA
mamba run -n bwa bwa index "${RACON2}"

# Align each Illumina library, each mate separately (required by Polypolish)
for srr in SRR10848483 SRR10848484; do
    for read in 1 2; do
        log_info "Aligning ${srr} R${read}..."
        mamba run -n bwa bwa mem \
            -t "${THREADS}" -a \
            "${RACON2}" \
            "${DECON_DIR}/${srr}_unclassified_${read}.fastq.gz" \
            > "${POLISH_DIR}/alignments_${srr}_${read}.sam" \
            2> "${LOG_DIR}/bwa_${srr}_${read}.stderr.log"
    done
done

# Filter alignments (removes spurious multi-mappers)
for srr in SRR10848483 SRR10848484; do
    mamba run -n polypolish polypolish filter \
        --in1  "${POLISH_DIR}/alignments_${srr}_1.sam" \
        --in2  "${POLISH_DIR}/alignments_${srr}_2.sam" \
        --out1 "${POLISH_DIR}/filtered_${srr}_1.sam" \
        --out2 "${POLISH_DIR}/filtered_${srr}_2.sam"
done

# Polish using all four filtered SAMs jointly
mamba run -n polypolish polypolish polish \
    "${RACON2}" \
    "${POLISH_DIR}/filtered_SRR10848483_1.sam" \
    "${POLISH_DIR}/filtered_SRR10848483_2.sam" \
    "${POLISH_DIR}/filtered_SRR10848484_1.sam" \
    "${POLISH_DIR}/filtered_SRR10848484_2.sam" \
    > "${POLISH_DIR}/pacbio_flye_assembly.racon2.polished.fasta" \
    2> "${LOG_DIR}/polypolish.stderr.log"

# Remove SAMs and BWA index files to save disk space
rm -f "${POLISH_DIR}"/alignments_*.sam \
      "${POLISH_DIR}"/filtered_*.sam \
      "${RACON2}".{amb,ann,bwt,pac,sa}

log_info "Polypolish → ${POLISH_DIR}/pacbio_flye_assembly.racon2.polished.fasta"

# — — — — — — — — — — — — — — — — — — — — — — — — — —
# 6d. NEXTPOLISH ROUNDS 1 AND 2 (final short-read correction)
# NextPolish is more accurate than Pilon for indel correction in homopolymer
# runs — a known PacBio error source. Applied after Polypolish as a final
# refinement. Two rounds applied given initial QV 34.4 and E:1.6% BUSCO
# internal stop codon rate, both indicating residual polishing potential.
# task = best : runs all available polishing sub-tasks.
# — — — — — — — — — — — — — — — — — — — — — — — — — —

log_step "Polishing: NextPolish (×2 rounds — final short-read correction)"

POLYPOLISHED="${POLISH_DIR}/pacbio_flye_assembly.racon2.polished.fasta"

# Build file-of-filenames for Illumina reads
cat > "${POLISH_DIR}/sgs.fofn" << EOF
${DECON_DIR}/SRR10848483_unclassified_1.fastq.gz ${DECON_DIR}/SRR10848483_unclassified_2.fastq.gz
${DECON_DIR}/SRR10848484_unclassified_1.fastq.gz ${DECON_DIR}/SRR10848484_unclassified_2.fastq.gz
EOF

# Round 1
cat > "${POLISH_DIR}/run_r1.cfg" << EOF
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 6
multithread_jobs = 3
genome = ${POLISH_DIR}/pacbio_flye_assembly.racon2.polished.fasta
genome_size = auto
workdir = ${POLISH_DIR}/nextpolish_work_r1
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ${POLISH_DIR}/sgs.fofn
sgs_options = -max_depth 100 -bwa
EOF

mamba run -n nextpolish39 nextPolish \
    "${POLISH_DIR}/run_r1.cfg" \
    > "${LOG_DIR}/nextpolish_r1.stdout.log" \
    2> "${LOG_DIR}/nextpolish_r1.stderr.log"

log_info "NextPolish round 1 → ${POLISH_DIR}/nextpolish_work_r1/genome.nextpolish.fasta"

# Round 2 (input = round 1 output)
cat > "${POLISH_DIR}/run_r2.cfg" << EOF
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 6
multithread_jobs = 3
genome = ${POLISH_DIR}/nextpolish_work_r1/genome.nextpolish.fasta
genome_size = auto
workdir = ${POLISH_DIR}/nextpolish_work_r2
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ${POLISH_DIR}/sgs.fofn
sgs_options = -max_depth 100 -bwa
EOF

mamba run -n nextpolish39 nextPolish \
    "${POLISH_DIR}/run_r2.cfg" \
    > "${LOG_DIR}/nextpolish_r2.stdout.log" \
    2> "${LOG_DIR}/nextpolish_r2.stderr.log"

FINAL_POLISHED="${POLISH_DIR}/nextpolish_work_r2/genome.nextpolish.fasta"

log_info "NextPolish round 2 → ${FINAL_POLISHED}"

#==============================================================================
# SUMMARY
#==============================================================================

echo ""
echo "========================================================================"
echo "STAGES 2–6 COMPLETE"
echo "========================================================================"
echo "Filtered reads         : ${FILTERED_DIR}/"
echo "Filtered QC reports    : ${QC_FILTERED_DIR}/"
echo "Contaminant reports    : ${CONTAMINANTS_DIR}/"
echo "Decontaminated reads   : ${DECON_DIR}/"
echo "Raw assemblies         : ${ASM_DIR}/"
echo "Final polished assembly: ${FINAL_POLISHED}"
echo ""
echo "Next steps:"
echo "  1. Run 03_evaluate.sh  — BUSCO + QUAST + Merqury on all assemblies"
echo "  2. Run 04_scaffold_mask.sh — RagTag scaffolding + RepeatModeler2/Masker"
echo "========================================================================"