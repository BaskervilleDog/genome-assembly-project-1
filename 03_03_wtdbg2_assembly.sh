#!/usr/bin/env bash
#==============================================================================
# WTDBG2 (REDBEAN) ASSEMBLY PIPELINE
#==============================================================================
# Author: Gianlucca de Urzêda Alves
# Date: 30/01/2026
# Description: Long-read genome assembly using wtdbg2
#==============================================================================

set -Eeuo pipefail

#==============================================================================
# CONFIGURATION
#==============================================================================

PACBIO_READS="./filtered_data/SRR10848482_subreads.fastq.gz"

GENOME_SIZE="40m"
THREADS=12
PRESET="rs"          # rs=PacBio RSII | sq=Sequel | ont=Nanopore
MIN_READ_LENGTH=5000
SUBSAMPLING=1
SPECIES="Trichoderma_harzianum"

OUTPUT_DIR="./wtdbg2_assembly"
LOG_DIR="${OUTPUT_DIR}/logs"
QC_DIR="${OUTPUT_DIR}/quality_assessment"

PREFIX="${OUTPUT_DIR}/assembly"
ASSEMBLY="${OUTPUT_DIR}/wtdbg2_assembly.fasta"
FILTERED_READS="${OUTPUT_DIR}/filtered_reads.fastq.gz"

#==============================================================================
# LOGGING FUNCTIONS
#==============================================================================

log_info()   { echo "[INFO]  $(date '+%F %T')  $*"; }
log_warn()   { echo "[WARN]  $(date '+%F %T')  $*" >&2; }
log_error()  { echo "[ERROR] $(date '+%F %T')  $*" >&2; }
log_section() {
    echo
    echo "========================================================================"
    echo "$*"
    echo "========================================================================"
}

#==============================================================================
# INITIAL CHECKS
#==============================================================================

log_section "WTDBG2 ASSEMBLY PIPELINE — START"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${QC_DIR}"

[[ -f "${PACBIO_READS}" ]] || {
    log_error "Input reads not found: ${PACBIO_READS}"
    exit 1
}

command -v mamba >/dev/null || {
    log_error "mamba not available"
    exit 1
}

mamba run -n wtdbg2 wtdbg2 --version &>/dev/null || {
    log_error "wtdbg2 not found in environment"
    exit 1
}

#==============================================================================
# STEP 1 — OPTIONAL READ FILTERING
#==============================================================================

log_section "STEP 1/5 — READ FILTERING (OPTIONAL)"

if [[ -f "${FILTERED_READS}" ]]; then
    log_info "Filtered reads already exist → skipping"
elif mamba run -n filtlong filtlong --version &>/dev/null; then
    log_info "Filtering reads with filtlong"
    mamba run -n filtlong filtlong \
        --min_length "${MIN_READ_LENGTH}" \
        --keep_percent 90 \
        "${PACBIO_READS}" | gzip > "${FILTERED_READS}"
else
    log_warn "Filtlong not available → using original reads"
    FILTERED_READS="${PACBIO_READS}"
fi

#==============================================================================
# STEP 2 — WTDBG2 GRAPH CONSTRUCTION
#==============================================================================

log_section "STEP 2/5 — GRAPH CONSTRUCTION"

if [[ -f "${PREFIX}.ctg.lay.gz" ]]; then
    log_info "Graph already exists → skipping"
else
    mamba run -n wtdbg2 wtdbg2 \
        -x "${PRESET}" \
        -g "${GENOME_SIZE}" \
        -i "${FILTERED_READS}" \
        -t "${THREADS}" \
        -L "${MIN_READ_LENGTH}" \
        -S "${SUBSAMPLING}" \
        -e 3 \
        -fo "${PREFIX}" \
        2>&1 | tee "${LOG_DIR}/wtdbg2_assembly.log"
fi

#==============================================================================
# STEP 3 — CONSENSUS GENERATION
#==============================================================================

log_section "STEP 3/5 — CONSENSUS GENERATION"

if [[ -f "${ASSEMBLY}" ]]; then
    log_info "Consensus already exists → skipping"
else
    mamba run -n wtdbg2 wtpoa-cns \
        -t "${THREADS}" \
        -i "${PREFIX}.ctg.lay.gz" \
        -fo "${ASSEMBLY}" \
        2>&1 | tee "${LOG_DIR}/wtpoa_consensus.log"
fi

#==============================================================================
# BASIC ASSEMBLY STATS
#==============================================================================

log_section "ASSEMBLY STATISTICS"

CONTIGS=$(grep -c "^>" "${ASSEMBLY}")
SIZE=$(awk '!/^>/ {sum+=length} END {print sum}' "${ASSEMBLY}")
LARGEST=$(awk '!/^>/ {print length}' "${ASSEMBLY}" | sort -nr | head -1)

log_info "Contigs:        ${CONTIGS}"
log_info "Total size:    ${SIZE} bp"
log_info "Largest contig:${LARGEST} bp"

#==============================================================================
# STEP 4 — QUALITY ASSESSMENT
#==============================================================================

log_section "STEP 4/5 — QUALITY ASSESSMENT"

if mamba run -n assembly_stats assembly-stats --version &>/dev/null; then
    mamba run -n assembly_stats assembly-stats "${ASSEMBLY}" \
        > "${QC_DIR}/assembly_stats.txt"
fi

if [[ ! -d "${QC_DIR}/quast" ]]; then
    mamba run -n quast quast \
        "${ASSEMBLY}" \
        -o "${QC_DIR}/quast" \
        --threads "${THREADS}" \
        --fungus \
        2>&1 | tee "${LOG_DIR}/quast.log"
fi

if [[ ! -d "${QC_DIR}/busco" ]]; then
    mamba run -n busco busco \
        -i "${ASSEMBLY}" \
        -o "${QC_DIR}/busco" \
        -m genome \
        -l hypocreaceae_odb12 \
        -c "${THREADS}" \
        --offline -f \
        2>&1 | tee "${LOG_DIR}/busco.log"
fi

#==============================================================================
# STEP 5 — SUMMARY REPORT
#==============================================================================

log_section "STEP 5/5 — REPORT"

REPORT="${OUTPUT_DIR}/WTDBG2_ASSEMBLY_REPORT.txt"

cat > "${REPORT}" <<EOF
WTDBG2 (REDBEAN) ASSEMBLY REPORT
========================================
Species: ${SPECIES}
Date: $(date)

INPUT:
PacBio reads: ${PACBIO_READS}

PARAMETERS:
Genome size: ${GENOME_SIZE}
Preset: ${PRESET}
Threads: ${THREADS}

ASSEMBLY:
Contigs: ${CONTIGS}
Total size: ${SIZE} bp
Largest contig: ${LARGEST} bp

NOTES:
This is a DRAFT assembly.
Polishing with Pilon is REQUIRED.
EOF

log_info "Report saved to ${REPORT}"

#==============================================================================
# FINAL
#==============================================================================

log_section "PIPELINE COMPLETE"

log_info "Assembly: ${ASSEMBLY}"
log_info "QC:       ${QC_DIR}"
log_info "Report:   ${REPORT}"

exit 0
