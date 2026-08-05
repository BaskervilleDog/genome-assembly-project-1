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
#
# Thin orchestrator — see lib/*.sh for the actual tool invocations.
#==============================================================================

set -euo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"

for f in "$PROJECT_DIR"/lib/*.sh; do
    # shellcheck source=/dev/null
    source "$f"
done

#==============================================================================
# CONFIGURATION
#==============================================================================

THREADS=18
KRAKEN_DB="${HOME}/databases/kraken2_db"

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
# STEP 1 — FILTER PACBIO READS (fastplong)
#==============================================================================
# fastplong is long-read-aware; standard short-read trimmers (fastp,
# Trimmomatic) are not appropriate for PacBio reads.
#
# Quality filtering is explicitly DISABLED (-Q) because quality scores are
# absent in this SRA dataset (Q=0 placeholders). Only length filtering is
# applied: reads <1 kb are too short to contribute to assembly contiguity
# and may introduce false overlaps in the assembler.

io::log_step "Filtering PacBio long reads (fastplong)"

filter::pacbio \
    "${DOWNLOAD_DIR}/SRR10848482_1.fastq.gz" \
    "${FILTERED_DIR}/pacbio/SRR10848482_1_filtered.fastq.gz" \
    1000 \
    "${QC_FILTERED_DIR}/SRR10848482_1.fastplong.html" \
    "${QC_FILTERED_DIR}/SRR10848482_1.fastplong.json" \
    "${LOG_DIR}/SRR10848482_1.fastplong.stdout.log" \
    "${LOG_DIR}/SRR10848482_1.fastplong.stderr.log"

io::log_info "PacBio filtering complete → ${FILTERED_DIR}/pacbio/"

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

io::log_step "Filtering Illumina short reads (fastp)"

for srr in SRR10848483 SRR10848484; do
    io::log_info "Filtering ${srr}..."
    filter::illumina \
        "${DOWNLOAD_DIR}/${srr}_1.fastq.gz" "${DOWNLOAD_DIR}/${srr}_2.fastq.gz" \
        "${FILTERED_DIR}/illumina/${srr}_1_filtered.fastq.gz" \
        "${FILTERED_DIR}/illumina/${srr}_2_filtered.fastq.gz" \
        "${THREADS}" 30 150 \
        "${QC_FILTERED_DIR}/${srr}.fastp.html" "${QC_FILTERED_DIR}/${srr}.fastp.json" \
        "${LOG_DIR}/${srr}.fastp.stdout.log" "${LOG_DIR}/${srr}.fastp.stderr.log"
done

io::log_info "Illumina filtering complete → ${FILTERED_DIR}/illumina/"

#==============================================================================
# STEP 3 — QC: FILTERED READS (SeqKit + FastQC + MultiQC)
#==============================================================================

io::log_step "QC: Filtered reads — SeqKit stats"

# PacBio filtered reads
qc::seqkit_stats \
    "${FILTERED_DIR}/pacbio/SRR10848482_1_filtered.fastq.gz" \
    "${QC_FILTERED_DIR}/SRR10848482_1_filtered_stats.txt" \
    "${LOG_DIR}/seqkit_filtered.stderr.log"

# Illumina filtered reads
for srr in SRR10848483 SRR10848484; do
    for read in 1 2; do
        qc::seqkit_stats \
            "${FILTERED_DIR}/illumina/${srr}_${read}_filtered.fastq.gz" \
            "${QC_FILTERED_DIR}/${srr}_${read}_filtered_stats.txt" \
            "${LOG_DIR}/seqkit_filtered.stderr.log"
    done
done

io::log_step "QC: Filtered Illumina reads — FastQC + MultiQC"

qc::fastqc \
    "${QC_FILTERED_DIR}/" "${THREADS}" \
    "${LOG_DIR}/fastqc_filtered.stdout.log" "${LOG_DIR}/fastqc_filtered.stderr.log" \
    "${FILTERED_DIR}/illumina/SRR10848483_1_filtered.fastq.gz" \
    "${FILTERED_DIR}/illumina/SRR10848483_2_filtered.fastq.gz" \
    "${FILTERED_DIR}/illumina/SRR10848484_1_filtered.fastq.gz" \
    "${FILTERED_DIR}/illumina/SRR10848484_2_filtered.fastq.gz"

qc::multiqc \
    "${QC_FILTERED_DIR}/" "${QC_FILTERED_DIR}/multiqc_report" \
    "${LOG_DIR}/multiqc_filtered.stdout.log" "${LOG_DIR}/multiqc_filtered.stderr.log"

io::log_info "Filtered QC reports written to ${QC_FILTERED_DIR}/"

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

io::log_step "Decontamination: Kraken2 taxonomic classification"

# PacBio long reads (single-end)
io::log_info "Classifying PacBio reads..."
decontaminate::kraken2_single \
    "${KRAKEN_DB}" "${THREADS}" \
    "${FILTERED_DIR}/pacbio/SRR10848482_1_filtered.fastq.gz" \
    "${CONTAMINANTS_DIR}/pacbio_report.k2report" \
    "${DECON_DIR}/SRR10848482_unclassified.fastq" \
    "${CONTAMINANTS_DIR}/SRR10848482_classified.fastq" \
    "${CONTAMINANTS_DIR}/pacbio_output.kraken2" \
    "${LOG_DIR}/kraken2_pacbio.stderr.log"

# Illumina paired-end reads (both runs)
for srr in SRR10848483 SRR10848484; do
    io::log_info "Classifying ${srr}..."
    decontaminate::kraken2_paired \
        "${KRAKEN_DB}" "${THREADS}" \
        "${FILTERED_DIR}/illumina/${srr}_1_filtered.fastq.gz" \
        "${FILTERED_DIR}/illumina/${srr}_2_filtered.fastq.gz" \
        "${CONTAMINANTS_DIR}/${srr}_report.k2report" \
        "${DECON_DIR}/${srr}_unclassified#.fastq" \
        "${CONTAMINANTS_DIR}/${srr}_classified#.fastq" \
        "${CONTAMINANTS_DIR}/${srr}_output.kraken2" \
        "${LOG_DIR}/kraken2_${srr}.stderr.log"
done

# Compress all decontaminated FASTQ files
io::log_info "Compressing decontaminated reads..."
pigz -p "${THREADS}" "${DECON_DIR}"/*.fastq

io::log_info "Decontaminated reads → ${DECON_DIR}/"
io::log_info "Contaminant reports → ${CONTAMINANTS_DIR}/"

#==============================================================================
# STEP 5 — GENOME ASSEMBLY (Flye + Raven + wtdbg2)
#==============================================================================
# Three assemblers with fundamentally different graph approaches are run on
# the same decontaminated PacBio reads. BUSCO and QUAST are applied to all
# three outputs (in 03_quality_assessment.sh) before selecting the best for
# polishing.
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

io::log_step "Assembly: Flye (repeat graph)"

assemble::flye \
    "${PACBIO_CLEAN}" 40m "${ASM_DIR}/flye" "${THREADS}" \
    "${LOG_DIR}/SRR10848482_flye.stdout.log" "${LOG_DIR}/SRR10848482_flye.stderr.log"

io::log_info "Flye assembly → ${ASM_DIR}/flye/assembly.fasta"

# — — — — — — — — — — — — — — — — — — — — — — — — — —
# 5b. RAVEN (string graph / OLC)
# Fast alternative with internal Racon polishing. Memory-efficient and
# often produces high contiguity on fungal genomes, providing a strong
# benchmark for comparison.
# — — — — — — — — — — — — — — — — — — — — — — — — — —

io::log_step "Assembly: Raven (string graph)"

assemble::raven \
    "${PACBIO_CLEAN}" "${THREADS}" \
    "${ASM_DIR}/raven/raven_assembly.fasta" \
    "${LOG_DIR}/SRR10848482_raven.stderr.log"

io::log_info "Raven assembly → ${ASM_DIR}/raven/raven_assembly.fasta"

# — — — — — — — — — — — — — — — — — — — — — — — — — —
# 5c. WTDBG2 (fuzzy de Bruijn graph)
# Structurally distinct from Flye and Raven, increasing the probability
# of identifying the best assembly through empirical evaluation.
# Consensus generation (wtpoa-cns) is run separately after graph construction.
# -x rs : PacBio RS2/Sequel error-tolerance preset
# — — — — — — — — — — — — — — — — — — — — — — — — — —

io::log_step "Assembly: wtdbg2 (fuzzy de Bruijn graph)"

assemble::wtdbg2 \
    "${PACBIO_CLEAN}" 40m "${THREADS}" \
    "${ASM_DIR}/wtdbg2/wtdbg2" "${ASM_DIR}/wtdbg2/wtdbg2_assembly.fasta" \
    "${LOG_DIR}/SRR10848482_wtdbg2.stdout.log" "${LOG_DIR}/SRR10848482_wtdbg2.stderr.log"

io::log_info "wtdbg2 assembly → ${ASM_DIR}/wtdbg2/wtdbg2_assembly.fasta"

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
# Run 03_quality_assessment.sh's raw-assembly evaluation before this step
# to confirm the selection.

FLYE_ASM="${ASM_DIR}/flye/assembly.fasta"

# — — — — — — — — — — — — — — — — — — — — — — — — — —
# 6a. RACON ROUND 1 (minimap2 → Racon)
# minimap2 -x map-pb : preset calibrated for PacBio CLR gap penalties.
# First Racon round corrects the largest indel/substitution errors.
# — — — — — — — — — — — — — — — — — — — — — — — — — —

io::log_step "Polishing Round 1: minimap2 + Racon (long reads)"

RACON1="${POLISH_DIR}/pacbio_flye_assembly.racon1.fasta"
polish::racon_round \
    "${FLYE_ASM}" "${PACBIO_CLEAN}" map-pb "${THREADS}" \
    "${POLISH_DIR}/pacbio_flye_alignments_r1.paf" "${RACON1}" \
    "${LOG_DIR}/minimap2_racon_r1.stderr.log" "${LOG_DIR}/racon_r1.stderr.log"

io::log_info "Racon round 1 → ${RACON1}"

# — — — — — — — — — — — — — — — — — — — — — — — — — —
# 6b. RACON ROUND 2 (minimap2 → Racon)
# Second round maps to the round-1 output and catches residual errors
# in regions with poor alignment coverage in round 1.
# A third round yields diminishing returns and is not applied.
# — — — — — — — — — — — — — — — — — — — — — — — — — —

io::log_step "Polishing Round 2: minimap2 + Racon (long reads)"

RACON2="${POLISH_DIR}/pacbio_flye_assembly.racon2.fasta"
polish::racon_round \
    "${RACON1}" "${PACBIO_CLEAN}" map-pb "${THREADS}" \
    "${POLISH_DIR}/pacbio_flye_alignments_r2.paf" "${RACON2}" \
    "${LOG_DIR}/minimap2_racon_r2.stderr.log" "${LOG_DIR}/racon_r2.stderr.log"

io::log_info "Racon round 2 → ${RACON2}"

# — — — — — — — — — — — — — — — — — — — — — — — — — —
# 6c. POLYPOLISH (BWA-MEM -a → Polypolish)
# Polypolish uses ALL alignments (bwa mem -a) rather than only uniquely
# mapping reads, enabling polishing of repetitive regions that Pilon skips.
# Both Illumina runs (SRR10848483 + SRR10848484) are aligned separately
# and jointly applied for maximum coverage (~60×).
# — — — — — — — — — — — — — — — — — — — — — — — — — —

io::log_step "Polishing: Polypolish (short reads — all alignments)"

POLYPOLISHED="${POLISH_DIR}/pacbio_flye_assembly.racon2.polished.fasta"
polish::polypolish \
    "${RACON2}" "${THREADS}" "${POLISH_DIR}" "${POLYPOLISHED}" "${LOG_DIR}" \
    "${DECON_DIR}/SRR10848483_unclassified_1.fastq.gz" \
    "${DECON_DIR}/SRR10848483_unclassified_2.fastq.gz" SRR10848483 \
    "${DECON_DIR}/SRR10848484_unclassified_1.fastq.gz" \
    "${DECON_DIR}/SRR10848484_unclassified_2.fastq.gz" SRR10848484

io::log_info "Polypolish → ${POLYPOLISHED}"

# — — — — — — — — — — — — — — — — — — — — — — — — — —
# 6d. NEXTPOLISH ROUNDS 1 AND 2 (final short-read correction)
# NextPolish is more accurate than Pilon for indel correction in homopolymer
# runs — a known PacBio error source. Applied after Polypolish as a final
# refinement. Two rounds applied given initial QV 34.4 and E:1.6% BUSCO
# internal stop codon rate, both indicating residual polishing potential.
# task = best : runs all available polishing sub-tasks.
# — — — — — — — — — — — — — — — — — — — — — — — — — —

io::log_step "Polishing: NextPolish (×2 rounds — final short-read correction)"

# Build file-of-filenames for Illumina reads
cat > "${POLISH_DIR}/sgs.fofn" << EOF
${DECON_DIR}/SRR10848483_unclassified_1.fastq.gz ${DECON_DIR}/SRR10848483_unclassified_2.fastq.gz
${DECON_DIR}/SRR10848484_unclassified_1.fastq.gz ${DECON_DIR}/SRR10848484_unclassified_2.fastq.gz
EOF

polish::nextpolish_round \
    "${POLYPOLISHED}" "${POLISH_DIR}/sgs.fofn" \
    "${POLISH_DIR}/nextpolish_work_r1" "${POLISH_DIR}/run_r1.cfg" \
    "${LOG_DIR}/nextpolish_r1.stdout.log" "${LOG_DIR}/nextpolish_r1.stderr.log"

io::log_info "NextPolish round 1 → ${POLISH_DIR}/nextpolish_work_r1/genome.nextpolish.fasta"

# Round 2 (input = round 1 output)
polish::nextpolish_round \
    "${POLISH_DIR}/nextpolish_work_r1/genome.nextpolish.fasta" "${POLISH_DIR}/sgs.fofn" \
    "${POLISH_DIR}/nextpolish_work_r2" "${POLISH_DIR}/run_r2.cfg" \
    "${LOG_DIR}/nextpolish_r2.stdout.log" "${LOG_DIR}/nextpolish_r2.stderr.log"

FINAL_POLISHED="${POLISH_DIR}/nextpolish_work_r2/genome.nextpolish.fasta"

io::log_info "NextPolish round 2 → ${FINAL_POLISHED}"

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
echo "Next step: run 03_quality_assessment.sh"
echo "  1. Evaluate raw assemblies (BUSCO/QUAST) — confirm Flye is still the"
echo "     best choice before this polishing run's output is scaffolded"
echo "  2. Evaluate + scaffold + repeat-mask the polished assembly"
echo "========================================================================"
