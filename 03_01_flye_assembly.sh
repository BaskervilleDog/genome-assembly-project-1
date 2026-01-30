#!/bin/bash
#==============================================================================
# FLYE ASSEMBLY + OPTIONAL PILON POLISHING PIPELINE
#==============================================================================
# Author: Gianlucca de Urzêda Alves
# Date: 30/01/2026
# Description:
#   - Flye assembly (PacBio mandatory)
#   - Pilon polishing (Illumina optional, 2 rounds)
#   - QUAST + BUSCO quality assessment
#==============================================================================

set -e
set -u
set -o pipefail

#==============================================================================
# CONFIGURATION
#==============================================================================

# ------------------
# Mandatory input
# ------------------
PACBIO_READS="./filtered_data/SRR10848482_subreads.fastq.gz"

# ------------------
# Optional Illumina inputs (leave empty to skip Pilon)
# ------------------
ILLUMINA_R1="${ILLUMINA_R1:-}"
ILLUMINA_R2="${ILLUMINA_R2:-}"
ILLUMINA2_R1="${ILLUMINA2_R1:-}"
ILLUMINA2_R2="${ILLUMINA2_R2:-}"

# ------------------
# Parameters
# ------------------
GENOME_SIZE="40m"
THREADS=12
PILON_MEMORY="32G"
SPECIES="Trichoderma_harzianum"

# ------------------
# Output
# ------------------
OUTPUT_DIR="./flye_assembly"
LOG_DIR="${OUTPUT_DIR}/logs"

#==============================================================================
# FUNCTIONS
#==============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

log_section() {
    echo ""
    echo "========================================================================"
    echo "$*"
    echo "========================================================================"
}

check_file() {
    if [[ ! -f "$1" ]]; then
        log_error "File not found: $1"
        exit 1
    fi
}

#==============================================================================
# INITIAL CHECKS
#==============================================================================

log_section "FLYE + OPTIONAL PILON PIPELINE - START"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

log_info "Checking mandatory input..."
check_file "${PACBIO_READS}"
log_info "✓ PacBio reads found"

#==============================================================================
# DETECT OPTIONAL PILON STEP
#==============================================================================

RUN_PILON=false

if [[ -n "${ILLUMINA_R1}" && -n "${ILLUMINA_R2}" && \
      -n "${ILLUMINA2_R1}" && -n "${ILLUMINA2_R2}" ]]; then

    log_info "Illumina reads detected → Pilon polishing ENABLED"

    check_file "${ILLUMINA_R1}"
    check_file "${ILLUMINA_R2}"
    check_file "${ILLUMINA2_R1}"
    check_file "${ILLUMINA2_R2}"

    RUN_PILON=true
else
    log_info "No (or incomplete) Illumina reads → Pilon polishing DISABLED"
fi

#==============================================================================
# STEP 1: FLYE ASSEMBLY (ALWAYS RUNS)
#==============================================================================

log_section "STEP 1/5: FLYE ASSEMBLY"

FLYE_OUT="${OUTPUT_DIR}/flye_output"
FLYE_ASM="${OUTPUT_DIR}/assembly.fasta"

if [[ -f "${FLYE_ASM}" ]]; then
    log_info "✓ Flye assembly already exists, skipping..."
else
    log_info "Running Flye assembly..."
    mamba run -n flye flye \
        --pacbio-raw "${PACBIO_READS}" \
        --out-dir "${FLYE_OUT}" \
        --genome-size "${GENOME_SIZE}" \
        --threads "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/flye.log"

    cp "${FLYE_OUT}/assembly.fasta" "${FLYE_ASM}"
fi

CONTIGS=$(grep -c "^>" "${FLYE_ASM}")
SIZE=$(grep -v "^>" "${FLYE_ASM}" | tr -d '\n' | wc -c)
log_info "Flye assembly: ${CONTIGS} contigs, ${SIZE} bp"

#==============================================================================
# OPTIONAL PILON POLISHING
#==============================================================================

FINAL_ASSEMBLY="${FLYE_ASM}"

if [[ "${RUN_PILON}" == true ]]; then

    #--------------------------------------------------------------------------
    # PILON ROUND 1
    #--------------------------------------------------------------------------
    log_section "STEP 2/5: PILON POLISHING - ROUND 1"

    ROUND1_DIR="${OUTPUT_DIR}/pilon_round1"
    mkdir -p "${ROUND1_DIR}"

    mamba run -n bwa bwa index "${FINAL_ASSEMBLY}" \
        2>&1 | tee "${LOG_DIR}/bwa_index_r1.log"

    for LIB in 1 2; do
        R1_VAR="ILLUMINA${LIB}_R1"
        R2_VAR="ILLUMINA${LIB}_R2"

        mamba run -n bwa bwa mem -t "${THREADS}" "${FINAL_ASSEMBLY}" \
            "${!R1_VAR}" "${!R2_VAR}" | \
            mamba run -n samtools samtools view -bS - | \
            mamba run -n samtools samtools sort -@ "${THREADS}" \
            -o "${ROUND1_DIR}/lib${LIB}.bam" -
    done

    mamba run -n samtools samtools merge -@ "${THREADS}" \
        "${ROUND1_DIR}/merged.bam" \
        "${ROUND1_DIR}"/lib*.bam

    mamba run -n samtools samtools index "${ROUND1_DIR}/merged.bam"

    export JAVA_TOOL_OPTIONS="-Xmx${PILON_MEMORY}"

    mamba run -n pilon pilon \
        --genome "${FINAL_ASSEMBLY}" \
        --frags "${ROUND1_DIR}/merged.bam" \
        --output pilon_round1 \
        --outdir "${ROUND1_DIR}" \
        --threads "${THREADS}" \
        --changes \
        2>&1 | tee "${LOG_DIR}/pilon_round1.log"

    FINAL_ASSEMBLY="${ROUND1_DIR}/pilon_round1.fasta"

    #--------------------------------------------------------------------------
    # PILON ROUND 2
    #--------------------------------------------------------------------------
    log_section "STEP 3/5: PILON POLISHING - ROUND 2"

    ROUND2_DIR="${OUTPUT_DIR}/pilon_round2"
    mkdir -p "${ROUND2_DIR}"

    mamba run -n bwa bwa index "${FINAL_ASSEMBLY}" \
        2>&1 | tee "${LOG_DIR}/bwa_index_r2.log"

    for LIB in 1 2; do
        R1_VAR="ILLUMINA${LIB}_R1"
        R2_VAR="ILLUMINA${LIB}_R2"

        mamba run -n bwa bwa mem -t "${THREADS}" "${FINAL_ASSEMBLY}" \
            "${!R1_VAR}" "${!R2_VAR}" | \
            mamba run -n samtools samtools view -bS - | \
            mamba run -n samtools samtools sort -@ "${THREADS}" \
            -o "${ROUND2_DIR}/lib${LIB}.bam" -
    done

    mamba run -n samtools samtools merge -@ "${THREADS}" \
        "${ROUND2_DIR}/merged.bam" \
        "${ROUND2_DIR}"/lib*.bam

    mamba run -n samtools samtools index "${ROUND2_DIR}/merged.bam"

    mamba run -n pilon pilon \
        --genome "${FINAL_ASSEMBLY}" \
        --frags "${ROUND2_DIR}/merged.bam" \
        --output pilon_round2 \
        --outdir "${ROUND2_DIR}" \
        --threads "${THREADS}" \
        --changes \
        2>&1 | tee "${LOG_DIR}/pilon_round2.log"

    FINAL_ASSEMBLY="${ROUND2_DIR}/pilon_round2.fasta"
fi

#==============================================================================
# STEP 4: QUALITY ASSESSMENT
#==============================================================================

log_section "STEP 4/5: QUALITY ASSESSMENT"

QC_DIR="${OUTPUT_DIR}/quality_assessment"
mkdir -p "${QC_DIR}"

mamba run -n quast quast \
    "${FINAL_ASSEMBLY}" \
    -o "${QC_DIR}/quast" \
    --threads "${THREADS}" \
    --fungus \
    2>&1 | tee "${LOG_DIR}/quast.log"

mamba run -n busco busco \
    -i "${FINAL_ASSEMBLY}" \
    -o "${QC_DIR}/busco" \
    -m genome \
    -l hypocreaceae_odb12 \
    -c "${THREADS}" \
    --offline \
    -f \
    2>&1 | tee "${LOG_DIR}/busco.log"

#==============================================================================
# STEP 5: SUMMARY REPORT
#==============================================================================

log_section "STEP 5/5: SUMMARY REPORT"

REPORT="${OUTPUT_DIR}/ASSEMBLY_REPORT.txt"

cat > "${REPORT}" << EOF
================================================================================
FLYE + OPTIONAL PILON ASSEMBLY REPORT
================================================================================
Species: ${SPECIES}
Date: $(date)

PacBio reads: ${PACBIO_READS}
Pilon polishing: $( [[ "${RUN_PILON}" == true ]] && echo "ENABLED" || echo "DISABLED" )

Final assembly:
${FINAL_ASSEMBLY}

QUAST:
${QC_DIR}/quast/report.html

BUSCO:
${QC_DIR}/busco/

Logs:
${LOG_DIR}
================================================================================
EOF

cat "${REPORT}"

log_section "PIPELINE COMPLETE"
log_info "Final assembly: ${FINAL_ASSEMBLY}"
log_info "Report: ${REPORT}"
