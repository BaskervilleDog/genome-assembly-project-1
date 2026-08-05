#!/usr/bin/env bash
#==============================================================================
# STAGE 1 — DATA DOWNLOAD & RAW READ QC
# Genome Assembly Pipeline: Trichoderma harzianum TW11
#==============================================================================
# Author : Gianlucca de Urzêda Alves
# Date   : 11-02-2026
#
# Purpose:
#   - Download raw FASTQ files from NCBI SRA
#   - Quality control raw PacBio reads (SeqKit + LongReadSum)
#   - Quality control raw Illumina reads (FastQC + MultiQC)
#
# Tools (mamba environments):
#   ncbi_tools  | seqfetcher
#   seqkit      | seqkit
#   longreadsum | longreadsum
#   fastqc      | fastqc
#   multiqc     | multiqc
#
# Outputs: 00_downloads/, 01_qc_raw/
#
# This script is intentionally thin: the actual tool invocations live in
# lib/*.sh (one function per tool call, module-prefixed - see the naming
# convention note in README.md) so they're unit-testable via tests/*.bats
# without needing 20 real mamba environments installed. This script only
# wires config + paths together and calls them in order. See
# scripts/check_deps.sh to verify the required environments exist before
# running this.
#==============================================================================

set -euo pipefail

BASE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"

for f in "$BASE_DIR"/lib/*.sh; do
    # shellcheck source=/dev/null
    source "$f"
done

#==============================================================================
# CONFIGURATION
#==============================================================================

THREADS=18

ACCESSIONS=(
  SRR10848482   # PacBio RS2/Sequel — long reads (7.1 Gb, ~178x)
  SRR10848483   # Illumina HiSeq 4000 — short reads run 1 (1.5 Gb, ~38x)
  SRR10848484   # Illumina HiSeq 4000 — short reads run 2 (0.95 Gb, ~23x)
)

SEQFETCHER_DIR="${BASE_DIR}/seqfetcher"
DOWNLOAD_DIR="${BASE_DIR}/00_downloads"
QC_RAW_DIR="${BASE_DIR}/01_qc_raw"
LOG_DIR="${BASE_DIR}/logs"

#==============================================================================
# DIRECTORY SETUP
#==============================================================================

mkdir -p \
    "${DOWNLOAD_DIR}" \
    "${QC_RAW_DIR}" \
    "${LOG_DIR}"

#==============================================================================
# STEP 1 — DOWNLOAD SEQUENCING DATA (SeqFetcher)
#==============================================================================
# SeqFetcher wraps fasterq-dump and related NCBI utilities into a single
# command that handles multiple accessions simultaneously, ensuring consistent
# file naming across all three runs.

io::log_step "Downloading SRA data (SeqFetcher)"

download::ensure_seqfetcher "${SEQFETCHER_DIR}"

io::log_info "Downloading: ${ACCESSIONS[*]}"
download::sra_reads "${SEQFETCHER_DIR}" "${DOWNLOAD_DIR}" "${ACCESSIONS[@]}"

io::log_info "All files downloaded to ${DOWNLOAD_DIR}/"

#==============================================================================
# STEP 2 — QC: RAW READS (SeqKit stats)
#==============================================================================
# SeqKit is used as a lightweight sanity check across all libraries before
# running heavier QC tools. The -a flag enables all extended metrics
# including N50 and Q-score statistics.
#
# NOTE: SeqKit will report Q=0 for PacBio reads downloaded from SRA.
# Quality scores are stripped during SRA archiving. This is expected.
# Quality-based filtering is NOT applied to long reads (see Script 2).

io::log_step "QC: Raw reads — SeqKit stats (all libraries)"

# PacBio long reads
io::log_info "SeqKit stats — PacBio long reads"
qc::seqkit_stats \
    "${DOWNLOAD_DIR}/SRR10848482_1.fastq.gz" \
    "${QC_RAW_DIR}/SRR10848482_1_raw_stats.txt" \
    "${LOG_DIR}/seqkit_raw_pacbio.stderr.log"

# Illumina short reads (both runs, both mates)
io::log_info "SeqKit stats — Illumina short reads"
for srr in SRR10848483 SRR10848484; do
    for read in 1 2; do
        qc::seqkit_stats \
            "${DOWNLOAD_DIR}/${srr}_${read}.fastq.gz" \
            "${QC_RAW_DIR}/${srr}_${read}_raw_stats.txt" \
            "${LOG_DIR}/seqkit_raw_illumina.stderr.log"
    done
done

io::log_info "SeqKit stats written to ${QC_RAW_DIR}/"

#==============================================================================
# STEP 3 — QC: RAW PACBIO READS (LongReadSum)
#==============================================================================
# For PacBio data where quality scores are absent, read-length distribution
# is the primary QC metric. LongReadSum provides read-length histograms and
# N50 distribution plots that FastQC cannot produce meaningfully for long reads.
# Used to confirm adequate long reads (>5 kb) for contiguous assembly.

io::log_step "QC: Raw PacBio reads — LongReadSum (read-length profiling)"

qc::longreadsum \
    "${DOWNLOAD_DIR}/SRR10848482_1.fastq.gz" \
    "${QC_RAW_DIR}/" \
    "${LOG_DIR}/longreadsum_raw.stdout.log" \
    "${LOG_DIR}/longreadsum_raw.stderr.log"

io::log_info "LongReadSum HTML report written to ${QC_RAW_DIR}/"

#==============================================================================
# STEP 4 — QC: RAW ILLUMINA READS (FastQC + MultiQC)
#==============================================================================
# FastQC detects adapter contamination, base quality drop-offs, and per-tile
# quality issues. Its HTML reports are consumed directly by MultiQC.
#
# MultiQC aggregates all four Illumina FASTQ reports into a single interactive
# report, making cross-run batch effects immediately apparent.

io::log_step "QC: Raw Illumina reads — FastQC"

qc::fastqc \
    "${QC_RAW_DIR}/" "${THREADS}" \
    "${LOG_DIR}/fastqc_raw.stdout.log" "${LOG_DIR}/fastqc_raw.stderr.log" \
    "${DOWNLOAD_DIR}/SRR10848483_1.fastq.gz" \
    "${DOWNLOAD_DIR}/SRR10848483_2.fastq.gz" \
    "${DOWNLOAD_DIR}/SRR10848484_1.fastq.gz" \
    "${DOWNLOAD_DIR}/SRR10848484_2.fastq.gz"

io::log_step "QC: Raw Illumina reads — MultiQC (aggregate report)"

qc::multiqc \
    "${QC_RAW_DIR}/" "${QC_RAW_DIR}/multiqc_report" \
    "${LOG_DIR}/multiqc_raw.stdout.log" "${LOG_DIR}/multiqc_raw.stderr.log"

io::log_info "MultiQC report written to ${QC_RAW_DIR}/multiqc_report/"

#==============================================================================
# SUMMARY
#==============================================================================

echo ""
echo "========================================================================"
echo "STAGE 1 COMPLETE"
echo "========================================================================"
echo "Downloaded reads : ${DOWNLOAD_DIR}/"
echo "Raw QC reports   : ${QC_RAW_DIR}/"
echo "Logs             : ${LOG_DIR}/"
echo ""
echo "Next step: run 02_assembly.sh"
echo "========================================================================"
