#!/usr/bin/env bash
#==============================================================================
# GENOME SIZE ESTIMATION — ILLUMINA READS
#==============================================================================
# Author: Gianlucca de Urzêda Alves
# Date: 02/02/2026
#
# Tools: jellyfish | GenomeScope2
# Description: Estimate genome size using k-mer spectra from Illumina PE libraries
#==============================================================================

set -euo pipefail

#==============================================================================
# CONFIGURATION
#==============================================================================

THREADS=12
KMER=21
HASH_SIZE=5G
MIN_COV=10
READ_LEN=150

DATA_DIR="./filtered_data"
OUT_DIR="./genome_size_estimation"
LOG_DIR="${OUT_DIR}/logs"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

# Libraries: sample_name R1 R2
LIBRARIES=(
  "SRR10848483 SRR10848483_1.fastq.gz SRR10848483_2.fastq.gz"
  "SRR10848484 SRR10848484_1.fastq.gz SRR10848484_2.fastq.gz"
)

#==============================================================================
# FUNCTIONS
#==============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
}

log_section() {
    echo ""
    echo "========================================================================"
    echo "$*"
    echo "========================================================================"
}

run_jellyfish() {
    local sample=$1
    local r1=$2
    local r2=$3
    
    log_info "Running Jellyfish for ${sample}"
    
    # Decompress to temporary directory
    local tmp_dir="${OUT_DIR}/tmp_${sample}"
    mkdir -p "${tmp_dir}"
    
    gunzip -c "${DATA_DIR}/${r1}" > "${tmp_dir}/R1.fastq"
    gunzip -c "${DATA_DIR}/${r2}" > "${tmp_dir}/R2.fastq"
    
    mamba run -n jellyfish jellyfish count \
        -C \
        -m "${KMER}" \
        -s "${HASH_SIZE}" \
        -t "${THREADS}" \
        -o "${OUT_DIR}/${sample}.jf" \
        "${tmp_dir}/R1.fastq" \
        "${tmp_dir}/R2.fastq" \
        2>&1 | tee "${LOG_DIR}/${sample}_count.log"
    
    mamba run -n jellyfish jellyfish histo \
        -t "${THREADS}" \
        "${OUT_DIR}/${sample}.jf" \
        > "${OUT_DIR}/${sample}.histo"
    
    # Clean up
    rm -rf "${tmp_dir}"
    rm -f "${OUT_DIR}/${sample}.jf"
}

estimate_genome_size() {
    local sample=$1
    log_info "Estimating genome size (k-mer sum ≥ ${MIN_COV}) for ${sample}"
    
    awk -v mincov="${MIN_COV}" '
        $1 >= mincov { sum += $1 * $2 }
        END { print sum }
    ' "${OUT_DIR}/${sample}.histo"
}

run_genomescope() {
    local sample=$1
    log_info "Running GenomeScope2 for ${sample}"
    
    mamba run -n genomescope genomescope2 \
        -i "${OUT_DIR}/${sample}.histo" \
        -o "${OUT_DIR}/${sample}_genomescope" \
        -k "${KMER}" \
        -p 1 \
        -l "${READ_LEN}" \
        --fitted_hist \
        2>&1 | tee "${LOG_DIR}/${sample}_genomescope.log"
}

#==============================================================================
# PIPELINE
#==============================================================================

log_section "GENOME SIZE ESTIMATION PIPELINE - START"

for lib in "${LIBRARIES[@]}"; do
    read -r SAMPLE R1 R2 <<< "${lib}"
    
    log_info "Processing library: ${SAMPLE}"
    
    # Skip if already processed
    if [[ -f "${OUT_DIR}/${SAMPLE}.histo" ]]; then
        log_info "Histogram already exists, skipping Jellyfish"
    else
        run_jellyfish "${SAMPLE}" "${R1}" "${R2}"
    fi
    
    TOTAL_KMERS=$(estimate_genome_size "${SAMPLE}")
    log_info "Total k-mers (≥${MIN_COV}x): ${TOTAL_KMERS}"
    
    if [[ -d "${OUT_DIR}/${SAMPLE}_genomescope" ]]; then
        log_info "GenomeScope already run, skipping"
    else
        run_genomescope "${SAMPLE}"
    fi
    
    echo "------------------------------------------------------------"
done

#==============================================================================
# SUMMARY REPORT
#==============================================================================

log_section "GENOME SIZE ESTIMATION - SUMMARY"

REPORT="${OUT_DIR}/GENOME_SIZE_REPORT.txt"

cat > "${REPORT}" << EOF
================================================================================
GENOME SIZE ESTIMATION REPORT
================================================================================
Date: $(date)
K-mer size: ${KMER}
Read length: ${READ_LEN}
Min coverage: ${MIN_COV}

RESULTS:
--------
EOF

for lib in "${LIBRARIES[@]}"; do
    read -r SAMPLE R1 R2 <<< "${lib}"
    
    cat >> "${REPORT}" << EOF

Library: ${SAMPLE}
-----------------
Histogram: ${OUT_DIR}/${SAMPLE}.histo
GenomeScope: ${OUT_DIR}/${SAMPLE}_genomescope/

EOF
    
    if [[ -f "${OUT_DIR}/${SAMPLE}_genomescope/summary.txt" ]]; then
        cat "${OUT_DIR}/${SAMPLE}_genomescope/summary.txt" >> "${REPORT}"
    fi
done

cat >> "${REPORT}" << EOF

================================================================================
INTERPRETATION NOTES
================================================================================
Genome size estimation formula:
  Genome size = Total k-mers / Coverage peak

Coverage peak identified from GenomeScope2 output:
  - model.txt
  - summary.txt
  - k-mer spectra plots

Compare estimates across libraries for consistency.

Next steps:
  1. Review GenomeScope2 plots
  2. Check for bimodal distributions (heterozygosity)
  3. Use estimated genome size for assembly parameters
================================================================================
EOF

cat "${REPORT}"

log_info "Report saved to: ${REPORT}"
log_info "✓ Genome size estimation complete"