#!/usr/bin/env bash
#==============================================================================
# WTDBG2 (REDBEAN) ASSEMBLY PIPELINE WITH PILON POLISHING
#==============================================================================
# Author: Gianlucca de Urzêda Alves
# Date: 30/01/2026
# Modified: 02/02/2026 - Added contig-by-contig Pilon polishing
# Description: Long-read genome assembly using wtdbg2 with optional Illumina polishing
#==============================================================================

set -Eeuo pipefail

#==============================================================================
# CONFIGURATION
#==============================================================================

PACBIO_READS="./filtered_data/SRR10848482_subreads.fastq.gz"

# Optional Illumina (leave empty or unset to skip Pilon)
ILLUMINA_R1="${ILLUMINA_R1:-./filtered_data/SRR10848483_1.fastq.gz}"
ILLUMINA_R2="${ILLUMINA_R2:-./filtered_data/SRR10848483_2.fastq.gz}"
ILLUMINA2_R1="${ILLUMINA2_R1:-./filtered_data/SRR10848484_1.fastq.gz}"
ILLUMINA2_R2="${ILLUMINA2_R2:-./filtered_data/SRR10848484_2.fastq.gz}"

GENOME_SIZE="40m"
THREADS=12
PRESET="rs"          # rs=PacBio RSII | sq=Sequel | ont=Nanopore
MIN_READ_LENGTH=5000
SUBSAMPLING=1
SPECIES="Trichoderma_harzianum"

# Pilon parameters (contig-by-contig mode for low memory)
PILON_ROUNDS=2
PILON_MEMORY="8G"  # Per contig, not total

OUTPUT_DIR="./wtdbg2_assembly"
LOG_DIR="${OUTPUT_DIR}/logs"
QC_DIR="${OUTPUT_DIR}/quality_assessment"

PREFIX="${OUTPUT_DIR}/assembly"
ASSEMBLY="${OUTPUT_DIR}/wtdbg2_assembly.fasta"
FILTERED_READS="${OUTPUT_DIR}/filtered_reads.fastq.gz"
FINAL_ASSEMBLY=""

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

log_section "STEP 1/6 — READ FILTERING (OPTIONAL)"

# Check if reads are already filtered (multiple detection methods)
if [[ "${PACBIO_READS}" == *"filtered_data"* ]] || [[ "${PACBIO_READS}" == *"filtered"* ]]; then
    log_info "Input reads are pre-filtered (from filtered_data/ directory)"
    log_info "Using reads directly: ${PACBIO_READS}"
    FILTERED_READS="${PACBIO_READS}"
elif [[ -f "${FILTERED_READS}" ]]; then
    log_info "Filtered reads already exist at: ${FILTERED_READS}"
    log_info "Skipping filtlong"
elif mamba run -n filtlong filtlong --version &>/dev/null; then
    log_info "Filtering reads with filtlong (min_length=${MIN_READ_LENGTH}, keep_percent=90)"
    mamba run -n filtlong filtlong \
        --min_length "${MIN_READ_LENGTH}" \
        --keep_percent 90 \
        "${PACBIO_READS}" | gzip > "${FILTERED_READS}"
    log_info "✓ Filtering complete: ${FILTERED_READS}"
else
    log_warn "Filtlong not available → using original reads"
    log_warn "Consider running the QC/filtering pipeline first for better results"
    FILTERED_READS="${PACBIO_READS}"
fi

log_info "Reads for assembly: ${FILTERED_READS}"

#==============================================================================
# STEP 2 — WTDBG2 GRAPH CONSTRUCTION
#==============================================================================

log_section "STEP 2/6 — GRAPH CONSTRUCTION"

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

log_section "STEP 3/6 — CONSENSUS GENERATION"

# Check if consensus exists and is valid
if [[ -f "${ASSEMBLY}" ]] && [[ -s "${ASSEMBLY}" ]]; then
    log_info "Consensus already exists and is valid → skipping"
else
    if [[ -f "${ASSEMBLY}" ]] && [[ ! -s "${ASSEMBLY}" ]]; then
        log_warn "Consensus file exists but is empty, regenerating"
        rm -f "${ASSEMBLY}"
    fi
    
    log_info "Generating consensus sequences"
    mamba run -n wtdbg2 wtpoa-cns \
        -t "${THREADS}" \
        -i "${PREFIX}.ctg.lay.gz" \
        -fo "${ASSEMBLY}" \
        2>&1 | tee "${LOG_DIR}/wtpoa_consensus.log"
    
    # Verify consensus was created successfully
    if [[ ! -s "${ASSEMBLY}" ]]; then
        log_error "Consensus generation failed - output file is empty"
        log_error "Check the log: ${LOG_DIR}/wtpoa_consensus.log"
        exit 1
    fi
    
    log_info "✓ Consensus generation complete"
fi

#==============================================================================
# BASIC ASSEMBLY STATS
#==============================================================================

log_section "ASSEMBLY STATISTICS"

# Verify assembly file exists and is not empty
if [[ ! -f "${ASSEMBLY}" ]]; then
    log_error "Assembly file not found: ${ASSEMBLY}"
    exit 1
fi

if [[ ! -s "${ASSEMBLY}" ]]; then
    log_error "Assembly file is empty: ${ASSEMBLY}"
    log_error "Consensus generation may have failed. Check ${LOG_DIR}/wtpoa_consensus.log"
    exit 1
fi

# Calculate statistics with error handling
CONTIGS=$(grep -c "^>" "${ASSEMBLY}" || echo "0")
SIZE=$(awk '!/^>/ {sum+=length} END {print sum}' "${ASSEMBLY}" || echo "0")
LARGEST=$(awk '!/^>/ {print length}' "${ASSEMBLY}" | sort -nr | head -1 || echo "0")

if [[ "${CONTIGS}" == "0" ]]; then
    log_error "No contigs found in assembly file"
    log_error "The consensus generation step may have failed"
    log_error "Check the log: ${LOG_DIR}/wtpoa_consensus.log"
    exit 1
fi

log_info "Contigs:        ${CONTIGS}"
log_info "Total size:     ${SIZE} bp"
log_info "Largest contig: ${LARGEST} bp"

#==============================================================================
# STEP 4 — CHECK IF PILON SHOULD RUN
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
# STEP 5 — CONTIG-BY-CONTIG PILON POLISHING
#==============================================================================

CURRENT_ASSEMBLY="${ASSEMBLY}"

if [[ "${RUN_PILON}" == true ]]; then

    log_section "STEP 5/6 — CONTIG-BY-CONTIG PILON POLISHING"

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
                    log_warn "Pilon failed for contig ${CONTIG_NAME}, keeping original"
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
else
    log_section "STEP 5/6 — PILON POLISHING (SKIPPED)"
fi

FINAL_ASSEMBLY="${CURRENT_ASSEMBLY}"

#==============================================================================
# STEP 6 — QUALITY ASSESSMENT
#==============================================================================

log_section "STEP 6/6 — QUALITY ASSESSMENT"

if mamba run -n assembly_stats assembly-stats --version &>/dev/null; then
    mamba run -n assembly_stats assembly-stats "${FINAL_ASSEMBLY}" \
        > "${QC_DIR}/assembly_stats.txt"
fi

if [[ ! -d "${QC_DIR}/quast" ]]; then
    log_info "Running QUAST"
    mamba run -n quast quast \
        "${FINAL_ASSEMBLY}" \
        -o "${QC_DIR}/quast" \
        --threads "${THREADS}" \
        --fungus \
        2>&1 | tee "${LOG_DIR}/quast.log"
fi

if [[ ! -d "${QC_DIR}/busco" ]]; then
    log_info "Running BUSCO"
    mamba run -n busco busco \
        -i "${FINAL_ASSEMBLY}" \
        -o "${QC_DIR}/busco" \
        -m genome \
        -l hypocreaceae_odb12 \
        -c "${THREADS}" \
        --offline -f \
        2>&1 | tee "${LOG_DIR}/busco.log"
fi

#==============================================================================
# STEP 7 — SUMMARY REPORT
#==============================================================================

log_section "SUMMARY REPORT"

FINAL_CONTIGS=$(grep -c "^>" "${FINAL_ASSEMBLY}")
FINAL_SIZE=$(awk '!/^>/ {sum+=length} END {print sum}' "${FINAL_ASSEMBLY}")
FINAL_LARGEST=$(awk '!/^>/ {print length}' "${FINAL_ASSEMBLY}" | sort -nr | head -1)

REPORT="${OUTPUT_DIR}/WTDBG2_ASSEMBLY_REPORT.txt"

cat > "${REPORT}" <<EOF
WTDBG2 (REDBEAN) ASSEMBLY REPORT
========================================
Species: ${SPECIES}
Date: $(date)

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
- Preset: ${PRESET}
- Threads: ${THREADS}
- Pilon rounds: $( [[ "${RUN_PILON}" == true ]] && echo "${PILON_ROUNDS} (contig-by-contig mode)" || echo "0 (skipped)" )
- Pilon memory: ${PILON_MEMORY} per contig

DRAFT ASSEMBLY (wtdbg2):
- Contigs: ${CONTIGS}
- Total size: ${SIZE} bp
- Largest contig: ${LARGEST} bp

FINAL ASSEMBLY (after Pilon):
- File: ${FINAL_ASSEMBLY}
- Contigs: ${FINAL_CONTIGS}
- Total size: ${FINAL_SIZE} bp
- Largest contig: ${FINAL_LARGEST} bp

QC OUTPUTS:
- QUAST: ${QC_DIR}/quast/report.html
- BUSCO: ${QC_DIR}/busco/

NOTES:
$( [[ "${RUN_PILON}" == true ]] && echo "Assembly polished with Pilon (contig-by-contig for memory efficiency)" || echo "This is a DRAFT assembly - Pilon polishing recommended" )
EOF

log_info "Report saved to ${REPORT}"

#==============================================================================
# FINAL
#==============================================================================

log_section "PIPELINE COMPLETE"

log_info "Final assembly: ${FINAL_ASSEMBLY}"
log_info "QC directory:   ${QC_DIR}"
log_info "Report:         ${REPORT}"

exit 0