#!/bin/bash
#==============================================================================
# RAVEN ASSEMBLY PIPELINE (OPTIONAL PILON POLISHING)
#==============================================================================
# Author: Gianlucca de Urzêda Alves
# Date: 30/01/2026
# Description:
#   PacBio assembly with Raven.
#   Optional Illumina polishing with Pilon if short reads are provided.
#==============================================================================

set -euo pipefail

#==============================================================================
# CONFIGURATION
#==============================================================================

# --- Input data ---
PACBIO_READS="./filtered_data/SRR10848482_subreads.fastq.gz"

# Optional Illumina (leave empty or unset to skip Pilon)
ILLUMINA_R1="${ILLUMINA_R1:-./filtered_data/SRR10848483_1.fastq.gz}"
ILLUMINA_R2="${ILLUMINA_R2:-./filtered_data/SRR10848483_2.fastq.gz}"
ILLUMINA2_R1="${ILLUMINA2_R1:-./filtered_data/SRR10848484_1.fastq.gz}"
ILLUMINA2_R2="${ILLUMINA2_R2:-./filtered_data/SRR10848484_2.fastq.gz}"

# --- Parameters ---
THREADS=12
RAVEN_POLISHING_ROUNDS=3
KMER_SIZE=15
WINDOW_SIZE=5

PILON_ROUNDS=2
PILON_MEMORY="32G"

SPECIES="Trichoderma_harzianum"

# --- Output ---
OUTPUT_DIR="./raven_assembly"
LOG_DIR="${OUTPUT_DIR}/logs"

RAVEN_ASSEMBLY="${OUTPUT_DIR}/raven_assembly.fasta"
FINAL_ASSEMBLY=""

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
    [[ -f "$1" ]] || return 1
}

#==============================================================================
# STEP 0: INITIAL CHECKS
#==============================================================================

log_section "RAVEN ASSEMBLY PIPELINE - START"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

check_file "${PACBIO_READS}" || {
    log_error "PacBio reads not found: ${PACBIO_READS}"
    exit 1
}

mamba run -n raven raven --version &> /dev/null || {
    log_error "Raven not installed"
    exit 1
}

log_info "✓ Initial checks passed"

#==============================================================================
# STEP 1: RAVEN ASSEMBLY
#==============================================================================

log_section "STEP 1/5: RAVEN ASSEMBLY"

FILTERED_READS="${PACBIO_READS}"

if [[ -f "${RAVEN_ASSEMBLY}" ]]; then
    log_info "Raven assembly already exists — skipping"
else
    log_info "Running Raven assembly"

    mamba run -n raven raven \
        --threads "${THREADS}" \
        --polishing-rounds "${RAVEN_POLISHING_ROUNDS}" \
        -k "${KMER_SIZE}" \
        -w "${WINDOW_SIZE}" \
        --disable-checkpoints \
        "${FILTERED_READS}" \
        > "${RAVEN_ASSEMBLY}" \
        2> "${LOG_DIR}/raven.log"

    log_info "✓ Raven assembly complete"
fi

#==============================================================================
# STEP 2: CHECK IF PILON SHOULD RUN
#==============================================================================

RUN_PILON=false

if check_file "${ILLUMINA_R1}" && check_file "${ILLUMINA_R2}" && \
   check_file "${ILLUMINA2_R1}" && check_file "${ILLUMINA2_R2}"; then
    RUN_PILON=true
    log_info "Illumina reads detected → Pilon ENABLED"
else
    log_info "Illumina reads missing → Pilon DISABLED"
fi

#==============================================================================
# STEP 3: OPTIONAL PILON POLISHING
#==============================================================================

CURRENT_ASSEMBLY="${RAVEN_ASSEMBLY}"

if [[ "${RUN_PILON}" == true ]]; then

    log_section "STEP 3/5: PILON POLISHING"

    PILON_DIR="${OUTPUT_DIR}/pilon"
    mkdir -p "${PILON_DIR}"

    for ROUND in $(seq 1 "${PILON_ROUNDS}"); do
        ROUND_DIR="${PILON_DIR}/round_${ROUND}"
        mkdir -p "${ROUND_DIR}"

        OUT_FASTA="${ROUND_DIR}/pilon_round${ROUND}.fasta"

        if [[ -f "${OUT_FASTA}" ]]; then
            log_info "Pilon round ${ROUND} already complete"
            CURRENT_ASSEMBLY="${OUT_FASTA}"
            continue
        fi

        log_info "Starting Pilon round ${ROUND}"

        [[ -f "${CURRENT_ASSEMBLY}.bwt" ]] || \
            mamba run -n bwa bwa index "${CURRENT_ASSEMBLY}"

        # Map Illumina libraries
        mamba run -n bwa bwa mem -t "${THREADS}" "${CURRENT_ASSEMBLY}" \
            "${ILLUMINA_R1}" "${ILLUMINA_R2}" | \
            mamba run -n samtools samtools sort -@ "${THREADS}" \
            -o "${ROUND_DIR}/lib1.bam"

        mamba run -n bwa bwa mem -t "${THREADS}" "${CURRENT_ASSEMBLY}" \
            "${ILLUMINA2_R1}" "${ILLUMINA2_R2}" | \
            mamba run -n samtools samtools sort -@ "${THREADS}" \
            -o "${ROUND_DIR}/lib2.bam"

        mamba run -n samtools samtools merge -@ "${THREADS}" \
            "${ROUND_DIR}/merged.bam" \
            "${ROUND_DIR}/lib1.bam" "${ROUND_DIR}/lib2.bam"

        mamba run -n samtools samtools index "${ROUND_DIR}/merged.bam"

        JAVA_TOOL_OPTIONS="-Xmx${PILON_MEMORY}" \
        mamba run -n pilon pilon \
            --genome "${CURRENT_ASSEMBLY}" \
            --frags "${ROUND_DIR}/merged.bam" \
            --output "pilon_round${ROUND}" \
            --outdir "${ROUND_DIR}" \
            --threads "${THREADS}" \
            --changes \
            2>&1 | tee "${LOG_DIR}/pilon_round${ROUND}.log"

        CURRENT_ASSEMBLY="${OUT_FASTA}"
    done
fi

FINAL_ASSEMBLY="${CURRENT_ASSEMBLY}"

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
    --offline -f \
    2>&1 | tee "${LOG_DIR}/busco.log"

#==============================================================================
# STEP 5: SUMMARY REPORT
#==============================================================================

log_section "STEP 5/5: SUMMARY"

CONTIGS=$(grep -c "^>" "${FINAL_ASSEMBLY}")
SIZE=$(grep -v "^>" "${FINAL_ASSEMBLY}" | tr -d '\n' | wc -c)

REPORT="${OUTPUT_DIR}/RAVEN_ASSEMBLY_REPORT.txt"

cat > "${REPORT}" << EOF
RAVEN ASSEMBLY REPORT
====================
Species: ${SPECIES}

Final assembly:
${FINAL_ASSEMBLY}

Assembly statistics:
- Contigs: ${CONTIGS}
- Genome size: ${SIZE} bp

Raven polishing rounds: ${RAVEN_POLISHING_ROUNDS}
Pilon polishing rounds: $( [[ "${RUN_PILON}" == true ]] && echo "${PILON_ROUNDS}" || echo "0 (skipped)" )

QC outputs:
- QUAST: ${QC_DIR}/quast/report.html
- BUSCO: ${QC_DIR}/busco/
EOF

log_info "Pipeline completed successfully"
log_info "Final assembly: ${FINAL_ASSEMBLY}"
log_info "Report: ${REPORT}"
