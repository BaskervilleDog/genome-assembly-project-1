#!/bin/bash
#==============================================================================
# FLYE ASSEMBLY + PILON POLISHING PIPELINE
#==============================================================================
# Author: Gianlucca de Urzêda Alves
# Date: 02/02/2026
#
# Description:
#   - Flye assembly (PacBio mandatory)
#   - Memory-efficient Pilon polishing (Illumina optional, contig-by-contig)
#   - QUAST + BUSCO quality assessment
#==============================================================================

set -euo pipefail

#==============================================================================
# CONFIGURATION
#==============================================================================

# Input data
PACBIO_READS="./filtered_data/SRR10848482_subreads.fastq.gz"

# Optional Illumina (leave empty or unset to skip Pilon)
ILLUMINA_R1="${ILLUMINA_R1:-./filtered_data/SRR10848483_1.fastq.gz}"
ILLUMINA_R2="${ILLUMINA_R2:-./filtered_data/SRR10848483_2.fastq.gz}"
ILLUMINA2_R1="${ILLUMINA2_R1:-./filtered_data/SRR10848484_1.fastq.gz}"
ILLUMINA2_R2="${ILLUMINA2_R2:-./filtered_data/SRR10848484_2.fastq.gz}"

# Parameters
GENOME_SIZE="40m"
THREADS=12
PILON_ROUNDS=2
PILON_MEMORY="8G"  # Per contig for memory efficiency

SPECIES="Trichoderma_harzianum"

# Output
OUTPUT_DIR="./flye_assembly"
LOG_DIR="${OUTPUT_DIR}/logs"

FLYE_OUT="${OUTPUT_DIR}/flye_output"
FLYE_ASM="${OUTPUT_DIR}/assembly.fasta"
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

cleanup_bwa_index() {
    local fasta="$1"
    if [[ -f "${fasta}.bwt" ]]; then
        log_info "Removing BWA index files for: $(basename ${fasta})"
        rm -f "${fasta}".{amb,ann,bwt,pac,sa}
    fi
}

force_cleanup() {
    sync
    sleep 2
}

#==============================================================================
# INITIAL CHECKS
#==============================================================================

log_section "FLYE ASSEMBLY PIPELINE - START"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

check_file "${PACBIO_READS}" || {
    log_error "PacBio reads not found: ${PACBIO_READS}"
    exit 1
}

log_info "✓ Initial checks passed"

#==============================================================================
# DETECT PILON
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
# STEP 1: FLYE ASSEMBLY
#==============================================================================

log_section "STEP 1/4: FLYE ASSEMBLY"

# Check if reads are pre-filtered
if [[ "${PACBIO_READS}" == *"filtered_data"* ]] || [[ "${PACBIO_READS}" == *"filtered"* ]]; then
    log_info "Input reads are pre-filtered (from filtered_data/ directory)"
fi

log_info "Reads for assembly: ${PACBIO_READS}"

if [[ -f "${FLYE_ASM}" ]]; then
    log_info "Flye assembly already exists, skipping..."
else
    log_info "Running Flye assembly..."
    
    mamba run -n flye flye \
        --pacbio-raw "${PACBIO_READS}" \
        --out-dir "${FLYE_OUT}" \
        --genome-size "${GENOME_SIZE}" \
        --threads "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/flye.log"
    
    cp "${FLYE_OUT}/assembly.fasta" "${FLYE_ASM}"
    log_info "✓ Flye assembly complete"
fi

CONTIGS=$(grep -c "^>" "${FLYE_ASM}")
SIZE=$(grep -v "^>" "${FLYE_ASM}" | tr -d '\n' | wc -c)

log_info "Flye assembly: ${CONTIGS} contigs, ${SIZE} bp"

#==============================================================================
# STEP 2: PILON POLISHING (STANDARD MODE FOR FLYE)
#==============================================================================

CURRENT_ASSEMBLY="${FLYE_ASM}"

if [[ "${RUN_PILON}" == true ]]; then

    log_section "STEP 2/4: PILON POLISHING"

    for ROUND in 1 2; do
        ROUND_DIR="${OUTPUT_DIR}/pilon_round${ROUND}"
        mkdir -p "${ROUND_DIR}"
        
        OUT_FASTA="${ROUND_DIR}/pilon_round${ROUND}.fasta"
        
        if [[ -f "${OUT_FASTA}" ]]; then
            log_info "Pilon round ${ROUND} already complete"
            CURRENT_ASSEMBLY="${OUT_FASTA}"
            continue
        fi
        
        log_info "Starting Pilon round ${ROUND}/${PILON_ROUNDS}"
        
        # Index assembly
        cleanup_bwa_index "${CURRENT_ASSEMBLY}"
        log_info "Indexing assembly with BWA"
        mamba run -n bwa bwa index "${CURRENT_ASSEMBLY}" \
            2>&1 | tee -a "${LOG_DIR}/bwa_index_r${ROUND}.log"
        
        # Map both libraries
        log_info "Mapping Illumina library 1"
        mamba run -n bwa bwa mem -t "${THREADS}" "${CURRENT_ASSEMBLY}" \
            "${ILLUMINA_R1}" "${ILLUMINA_R2}" 2>> "${LOG_DIR}/pilon_round${ROUND}.log" | \
            mamba run -n samtools samtools sort -@ "${THREADS}" \
            -o "${ROUND_DIR}/lib1.bam" 2>> "${LOG_DIR}/pilon_round${ROUND}.log"
        mamba run -n samtools samtools index "${ROUND_DIR}/lib1.bam"
        
        log_info "Mapping Illumina library 2"
        mamba run -n bwa bwa mem -t "${THREADS}" "${CURRENT_ASSEMBLY}" \
            "${ILLUMINA2_R1}" "${ILLUMINA2_R2}" 2>> "${LOG_DIR}/pilon_round${ROUND}.log" | \
            mamba run -n samtools samtools sort -@ "${THREADS}" \
            -o "${ROUND_DIR}/lib2.bam" 2>> "${LOG_DIR}/pilon_round${ROUND}.log"
        mamba run -n samtools samtools index "${ROUND_DIR}/lib2.bam"
        
        # Merge
        log_info "Merging BAM files"
        mamba run -n samtools samtools merge -f -@ "${THREADS}" \
            "${ROUND_DIR}/merged.bam" \
            "${ROUND_DIR}/lib1.bam" "${ROUND_DIR}/lib2.bam" \
            2>> "${LOG_DIR}/pilon_round${ROUND}.log"
        
        # Clean up individual BAMs
        rm -f "${ROUND_DIR}/lib1.bam" "${ROUND_DIR}/lib1.bam.bai"
        rm -f "${ROUND_DIR}/lib2.bam" "${ROUND_DIR}/lib2.bam.bai"
        
        mamba run -n samtools samtools index "${ROUND_DIR}/merged.bam"
        
        # Run Pilon
        log_info "Running Pilon (this may take a while)"
        unset JAVA_TOOL_OPTIONS
        
        JAVA_TOOL_OPTIONS="-Xmx${PILON_MEMORY}" \
            mamba run -n pilon pilon \
            --genome "${CURRENT_ASSEMBLY}" \
            --frags "${ROUND_DIR}/merged.bam" \
            --output "pilon_round${ROUND}" \
            --outdir "${ROUND_DIR}" \
            --fix bases \
            --changes \
            2>&1 | tee -a "${LOG_DIR}/pilon_round${ROUND}.log"
        
        # Cleanup
        log_info "Cleaning up round ${ROUND}"
        rm -f "${ROUND_DIR}/merged.bam" "${ROUND_DIR}/merged.bam.bai"
        cleanup_bwa_index "${CURRENT_ASSEMBLY}"
        force_cleanup
        
        CURRENT_ASSEMBLY="${OUT_FASTA}"
        log_info "✓ Pilon round ${ROUND} complete"
    done
    
    log_info "✓ All Pilon polishing rounds complete"
else
    log_section "STEP 2/4: PILON POLISHING (SKIPPED)"
fi

FINAL_ASSEMBLY="${CURRENT_ASSEMBLY}"

#==============================================================================
# STEP 3: QUALITY ASSESSMENT
#==============================================================================

log_section "STEP 3/4: QUALITY ASSESSMENT"

QC_DIR="${OUTPUT_DIR}/quality_assessment"
mkdir -p "${QC_DIR}"

log_info "Running QUAST"
mamba run -n quast quast \
    "${FINAL_ASSEMBLY}" \
    -o "${QC_DIR}/quast" \
    --threads "${THREADS}" \
    --fungus \
    2>&1 | tee "${LOG_DIR}/quast.log"

log_info "Running BUSCO"
mamba run -n busco busco \
    -i "${FINAL_ASSEMBLY}" \
    -o "${QC_DIR}/busco" \
    -m genome \
    -l hypocreaceae_odb12 \
    -c "${THREADS}" \
    --offline -f \
    2>&1 | tee "${LOG_DIR}/busco.log"

log_info "✓ Quality assessment complete"

#==============================================================================
# STEP 4: SUMMARY REPORT
#==============================================================================

log_section "STEP 4/4: SUMMARY"

FINAL_CONTIGS=$(grep -c "^>" "${FINAL_ASSEMBLY}")
FINAL_SIZE=$(grep -v "^>" "${FINAL_ASSEMBLY}" | tr -d '\n' | wc -c)

REPORT="${OUTPUT_DIR}/FLYE_ASSEMBLY_REPORT.txt"

cat > "${REPORT}" << EOF
================================================================================
FLYE ASSEMBLY REPORT
================================================================================
Species: ${SPECIES}
Date: $(date '+%Y-%m-%d %H:%M:%S')

INPUT:
- PacBio reads: ${PACBIO_READS}
$( [[ "${RUN_PILON}" == true ]] && cat << INPUTS
- Illumina R1: ${ILLUMINA_R1}
- Illumina R2: ${ILLUMINA_R2}
- Illumina2 R1: ${ILLUMINA2_R1}
- Illumina2 R2: ${ILLUMINA2_R2}
INPUTS
)

PARAMETERS:
- Genome size: ${GENOME_SIZE}
- Threads: ${THREADS}
- Pilon rounds: $( [[ "${RUN_PILON}" == true ]] && echo "${PILON_ROUNDS}" || echo "0 (skipped)" )
- Pilon memory: ${PILON_MEMORY}

DRAFT ASSEMBLY (Flye):
- Contigs: ${CONTIGS}
- Total size: ${SIZE} bp

FINAL ASSEMBLY (after Pilon):
- File: ${FINAL_ASSEMBLY}
- Contigs: ${FINAL_CONTIGS}
- Total size: ${FINAL_SIZE} bp

QC OUTPUTS:
- QUAST: ${QC_DIR}/quast/report.html
- BUSCO: ${QC_DIR}/busco/

NOTES:
$( [[ "${RUN_PILON}" == true ]] && echo "Assembly polished with Pilon" || echo "This is a DRAFT assembly - Pilon polishing recommended" )
================================================================================
EOF

cat "${REPORT}"

log_info "Report saved to: ${REPORT}"
log_info "✓ Pipeline completed successfully"
log_info "Final assembly: ${FINAL_ASSEMBLY}"

log_section "PIPELINE COMPLETE"