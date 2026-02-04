#!/bin/bash
#==============================================================================
# RAVEN ASSEMBLY PIPELINE (CONTIG-BY-CONTIG PILON POLISHING)
#==============================================================================
# Author: Gianlucca de Urzêda Alves
# Date: 30/01/2026
# Modified: 02/02/2026 - Ultra memory-efficient contig-by-contig processing
# Description:
#   PacBio assembly with Raven.
#   Ultra-low memory Illumina polishing by processing each contig separately.
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

# Process contigs separately to minimize memory
PILON_ROUNDS=2
PILON_MEMORY="8G"  # Much lower since we process one contig at a time

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

split_fasta_by_contig() {
    local input_fasta="$1"
    local output_dir="$2"
    
    mkdir -p "${output_dir}"
    
    log_info "Splitting assembly into individual contigs"
    
    mamba run -n samtools samtools faidx "${input_fasta}"
    
    # Get list of contig names
    cut -f1 "${input_fasta}.fai" > "${output_dir}/contig_list.txt"
    
    local contig_count=$(wc -l < "${output_dir}/contig_list.txt")
    log_info "Found ${contig_count} contigs to process"
    
    echo "${contig_count}"
}

extract_contig() {
    local input_fasta="$1"
    local contig_name="$2"
    local output_fasta="$3"
    
    mamba run -n samtools samtools faidx "${input_fasta}" "${contig_name}" > "${output_fasta}"
}

merge_polished_contigs() {
    local contig_list="$1"
    local contig_dir="$2"
    local output_fasta="$3"
    
    log_info "Merging polished contigs into final assembly"
    
    > "${output_fasta}"
    
    while IFS= read -r contig; do
        local polished_contig="${contig_dir}/${contig}_polished.fasta"
        if [[ -f "${polished_contig}" ]]; then
            cat "${polished_contig}" >> "${output_fasta}"
        else
            log_error "Missing polished contig: ${contig}"
            exit 1
        fi
    done < "${contig_list}"
    
    log_info "✓ Merged $(wc -l < ${contig_list}) contigs"
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

# Check if reads are already filtered (multiple detection methods)
if [[ "${PACBIO_READS}" == *"filtered_data"* ]] || [[ "${PACBIO_READS}" == *"filtered"* ]]; then
    log_info "Input reads are pre-filtered (from filtered_data/ directory)"
    log_info "Raven will use these reads directly"
    FILTERED_READS="${PACBIO_READS}"
else
    log_info "Using PacBio reads: ${PACBIO_READS}"
    log_info "Note: Raven performs internal quality filtering during assembly"
    FILTERED_READS="${PACBIO_READS}"
fi

log_info "Reads for assembly: ${FILTERED_READS}"

if [[ -f "${RAVEN_ASSEMBLY}" ]]; then
    log_info "Raven assembly already exists — skipping"
else
    log_info "Running Raven assembly with ${RAVEN_POLISHING_ROUNDS} polishing rounds"

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
# STEP 3: CONTIG-BY-CONTIG PILON POLISHING
#==============================================================================

CURRENT_ASSEMBLY="${RAVEN_ASSEMBLY}"

if [[ "${RUN_PILON}" == true ]]; then

    log_section "STEP 3/5: CONTIG-BY-CONTIG PILON POLISHING"

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

        log_info "Starting Pilon round ${ROUND}/${PILON_ROUNDS}"

        # Split assembly into contigs
        CONTIG_DIR="${ROUND_DIR}/contigs"
        CONTIG_COUNT=$(split_fasta_by_contig "${CURRENT_ASSEMBLY}" "${CONTIG_DIR}")
        
        # Create BAM file ONCE for all contigs
        cleanup_bwa_index "${CURRENT_ASSEMBLY}"
        
        log_info "Indexing assembly with BWA"
        mamba run -n bwa bwa index "${CURRENT_ASSEMBLY}" \
            2>&1 | tee -a "${LOG_DIR}/pilon_round${ROUND}.log"

        log_info "Mapping Illumina reads (this will be used for all contigs)"
        
        # Map both libraries
        mamba run -n bwa bwa mem -t "${THREADS}" "${CURRENT_ASSEMBLY}" \
            "${ILLUMINA_R1}" "${ILLUMINA_R2}" 2>> "${LOG_DIR}/pilon_round${ROUND}.log" | \
            mamba run -n samtools samtools sort -@ "${THREADS}" \
            -o "${ROUND_DIR}/lib1.bam" 2>> "${LOG_DIR}/pilon_round${ROUND}.log"
        mamba run -n samtools samtools index "${ROUND_DIR}/lib1.bam"

        mamba run -n bwa bwa mem -t "${THREADS}" "${CURRENT_ASSEMBLY}" \
            "${ILLUMINA2_R1}" "${ILLUMINA2_R2}" 2>> "${LOG_DIR}/pilon_round${ROUND}.log" | \
            mamba run -n samtools samtools sort -@ "${THREADS}" \
            -o "${ROUND_DIR}/lib2.bam" 2>> "${LOG_DIR}/pilon_round${ROUND}.log"
        mamba run -n samtools samtools index "${ROUND_DIR}/lib2.bam"

        # Merge BAM files
        mamba run -n samtools samtools merge -f -@ "${THREADS}" \
            "${ROUND_DIR}/merged.bam" \
            "${ROUND_DIR}/lib1.bam" "${ROUND_DIR}/lib2.bam" \
            2>> "${LOG_DIR}/pilon_round${ROUND}.log"
        
        # Clean up individual BAMs immediately
        rm -f "${ROUND_DIR}/lib1.bam" "${ROUND_DIR}/lib1.bam.bai"
        rm -f "${ROUND_DIR}/lib2.bam" "${ROUND_DIR}/lib2.bam.bai"
        
        mamba run -n samtools samtools index "${ROUND_DIR}/merged.bam"
        
        log_info "✓ Mapping complete, now processing ${CONTIG_COUNT} contigs individually"

        # Process each contig separately
        CONTIG_NUM=0
        while IFS= read -r CONTIG_NAME; do
            CONTIG_NUM=$((CONTIG_NUM + 1))
            
            log_info "Processing contig ${CONTIG_NUM}/${CONTIG_COUNT}: ${CONTIG_NAME}"
            
            CONTIG_FASTA="${CONTIG_DIR}/${CONTIG_NAME}.fasta"
            CONTIG_BAM="${CONTIG_DIR}/${CONTIG_NAME}.bam"
            POLISHED_FASTA="${CONTIG_DIR}/${CONTIG_NAME}_polished.fasta"
            
            # Extract contig sequence
            extract_contig "${CURRENT_ASSEMBLY}" "${CONTIG_NAME}" "${CONTIG_FASTA}"
            
            # Extract reads mapping to this contig
            mamba run -n samtools samtools view -b -h "${ROUND_DIR}/merged.bam" \
                "${CONTIG_NAME}" > "${CONTIG_BAM}"
            mamba run -n samtools samtools index "${CONTIG_BAM}"
            
            # Run Pilon on this contig only
            unset JAVA_TOOL_OPTIONS
            
            JAVA_TOOL_OPTIONS="-Xmx${PILON_MEMORY}" \
                mamba run -n pilon pilon \
                --genome "${CONTIG_FASTA}" \
                --frags "${CONTIG_BAM}" \
                --output "${CONTIG_NAME}_polished" \
                --outdir "${CONTIG_DIR}" \
                --fix bases \
                --changes \
                2>&1 | tee -a "${LOG_DIR}/pilon_round${ROUND}_${CONTIG_NAME}.log" || {
                    log_error "Pilon failed for contig ${CONTIG_NAME}"
                    # If Pilon fails, keep original contig
                    cp "${CONTIG_FASTA}" "${POLISHED_FASTA}"
                }
            
            # Clean up intermediate files for this contig
            rm -f "${CONTIG_FASTA}" "${CONTIG_BAM}" "${CONTIG_BAM}.bai"
            
        done < "${CONTIG_DIR}/contig_list.txt"
        
        log_info "✓ All contigs polished, merging results"
        
        # Merge all polished contigs
        merge_polished_contigs "${CONTIG_DIR}/contig_list.txt" "${CONTIG_DIR}" "${OUT_FASTA}"
        
        # Clean up
        log_info "Cleaning up round ${ROUND}"
        rm -f "${ROUND_DIR}/merged.bam" "${ROUND_DIR}/merged.bam.bai"
        rm -rf "${CONTIG_DIR}"
        
        if [[ ${ROUND} -gt 1 ]]; then
            PREV_ASSEMBLY="${PILON_DIR}/round_$((ROUND-1))/pilon_round$((ROUND-1)).fasta"
            cleanup_bwa_index "${PREV_ASSEMBLY}"
        fi
        
        cleanup_bwa_index "${CURRENT_ASSEMBLY}"
        force_cleanup
        
        CURRENT_ASSEMBLY="${OUT_FASTA}"
        log_info "✓ Pilon round ${ROUND} complete"
    done

    log_info "✓ All Pilon polishing rounds complete"
fi

FINAL_ASSEMBLY="${CURRENT_ASSEMBLY}"

#==============================================================================
# STEP 4: QUALITY ASSESSMENT
#==============================================================================

log_section "STEP 4/5: QUALITY ASSESSMENT"

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
Date: $(date '+%Y-%m-%d %H:%M:%S')

Final assembly:
${FINAL_ASSEMBLY}

Assembly statistics:
- Contigs: ${CONTIGS}
- Genome size: ${SIZE} bp

Pipeline parameters:
- Raven polishing rounds: ${RAVEN_POLISHING_ROUNDS}
- Pilon polishing rounds: $( [[ "${RUN_PILON}" == true ]] && echo "${PILON_ROUNDS} (contig-by-contig mode)" || echo "0 (skipped)" )
- Pilon memory limit: ${PILON_MEMORY} per contig
- Threads: ${THREADS}

QC outputs:
- QUAST: ${QC_DIR}/quast/report.html
- BUSCO: ${QC_DIR}/busco/

Input files:
- PacBio: ${PACBIO_READS}
$( [[ "${RUN_PILON}" == true ]] && cat << INPUTS
- Illumina R1: ${ILLUMINA_R1}
- Illumina R2: ${ILLUMINA_R2}
- Illumina2 R1: ${ILLUMINA2_R1}
- Illumina2 R2: ${ILLUMINA2_R2}
INPUTS
)
EOF

log_info "✓ Pipeline completed successfully"
log_info "Final assembly: ${FINAL_ASSEMBLY}"
log_info "Report: ${REPORT}"

log_section "PIPELINE COMPLETE"