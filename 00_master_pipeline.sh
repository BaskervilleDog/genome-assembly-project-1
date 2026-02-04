#!/bin/bash
#==============================================================================
# MASTER GENOME ASSEMBLY PIPELINE - Trichoderma harzianum
#==============================================================================
# Author: Gianlucca de Urzêda Alves
# Date: 02/02/2026
#
# Description:
#   Complete genome assembly pipeline from raw data to final assembly selection
#
# Pipeline Steps:
#   1. Download & QC
#   2. Genome size estimation
#   3. Assembly (Flye, Raven, wtdbg2)
#   4. Scaffolding (RagTag)
#   5. Assembly comparison
#
# Usage:
#   ./00_master_pipeline.sh [--skip-download] [--skip-genome-size] [--assemblers flye,raven,wtdbg2]
#==============================================================================

set -euo pipefail

#==============================================================================
# CONFIGURATION
#==============================================================================

SPECIES="Trichoderma_harzianum"
THREADS=12

# Default: run all steps
SKIP_DOWNLOAD=false
SKIP_GENOME_SIZE=false
ASSEMBLERS="flye,raven,wtdbg2"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-download)
            SKIP_DOWNLOAD=true
            shift
            ;;
        --skip-genome-size)
            SKIP_GENOME_SIZE=true
            shift
            ;;
        --assemblers)
            ASSEMBLERS="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--skip-download] [--skip-genome-size] [--assemblers flye,raven,wtdbg2]"
            exit 1
            ;;
    esac
done

#==============================================================================
# FUNCTIONS
#==============================================================================

log_section() {
    echo ""
    echo "========================================================================"
    echo "MASTER PIPELINE: $*"
    echo "========================================================================"
    echo ""
}

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

check_script() {
    if [[ ! -f "$1" ]]; then
        log_error "Required script not found: $1"
        exit 1
    fi
    chmod +x "$1"
}

#==============================================================================
# PIPELINE START
#==============================================================================

log_section "PIPELINE START"
log_info "Species: ${SPECIES}"
log_info "Threads: ${THREADS}"
log_info "Assemblers: ${ASSEMBLERS}"
echo ""

START_TIME=$(date +%s)

#==============================================================================
# STEP 1: DOWNLOAD & QC
#==============================================================================

if [[ "${SKIP_DOWNLOAD}" == false ]]; then
    log_section "STEP 1/5: DOWNLOAD & QC"
    
    check_script "01_download_and_qc.sh"
    
    log_info "Running download and quality control pipeline..."
    ./01_download_and_qc.sh
    
    log_info "✓ Download and QC complete"
else
    log_info "Skipping download step (--skip-download)"
fi

#==============================================================================
# STEP 2: GENOME SIZE ESTIMATION
#==============================================================================

if [[ "${SKIP_GENOME_SIZE}" == false ]]; then
    log_section "STEP 2/5: GENOME SIZE ESTIMATION"
    
    check_script "02_genome_size_estimation.sh"
    
    log_info "Running genome size estimation..."
    ./02_genome_size_estimation.sh
    
    log_info "✓ Genome size estimation complete"
    log_info "Review results in: ./genome_size_estimation/"
else
    log_info "Skipping genome size estimation (--skip-genome-size)"
fi

#==============================================================================
# STEP 3: ASSEMBLIES
#==============================================================================

log_section "STEP 3/5: GENOME ASSEMBLY"

IFS=',' read -ra ASSEMBLER_ARRAY <<< "${ASSEMBLERS}"

for assembler in "${ASSEMBLER_ARRAY[@]}"; do
    case "${assembler}" in
        flye)
            log_info "Running Flye assembly..."
            check_script "03_01_flye_assembly.sh"
            ./03_01_flye_assembly.sh
            log_info "✓ Flye assembly complete"
            ;;
        raven)
            log_info "Running Raven assembly..."
            check_script "03_02_raven_assembly.sh"
            ./03_02_raven_assembly.sh
            log_info "✓ Raven assembly complete"
            ;;
        wtdbg2)
            log_info "Running wtdbg2 assembly..."
            check_script "03_03_wtdbg2_assembly.sh"
            ./03_03_wtdbg2_assembly.sh
            log_info "✓ wtdbg2 assembly complete"
            ;;
        *)
            log_error "Unknown assembler: ${assembler}"
            log_error "Valid options: flye, raven, wtdbg2"
            exit 1
            ;;
    esac
done

#==============================================================================
# STEP 4: SCAFFOLDING
#==============================================================================

log_section "STEP 4/5: REFERENCE-GUIDED SCAFFOLDING"

check_script "04_ragtag_scaffolding.sh"

for assembler in "${ASSEMBLER_ARRAY[@]}"; do
    case "${assembler}" in
        flye)
            if [[ -f "./flye_assembly/pilon_round2/pilon_round2.fasta" ]]; then
                log_info "Scaffolding Flye assembly..."
                ./04_ragtag_scaffolding.sh flye_assembly pilon_round2/pilon_round2.fasta
            elif [[ -f "./flye_assembly/assembly.fasta" ]]; then
                log_info "Scaffolding Flye assembly (no Pilon)..."
                ./04_ragtag_scaffolding.sh flye_assembly assembly.fasta
            fi
            ;;
        raven)
            if [[ -f "./raven_assembly/pilon/round_2/pilon_round2.fasta" ]]; then
                log_info "Scaffolding Raven assembly..."
                ./04_ragtag_scaffolding.sh raven_assembly pilon/round_2/pilon_round2.fasta
            elif [[ -f "./raven_assembly/raven_assembly.fasta" ]]; then
                log_info "Scaffolding Raven assembly (no Pilon)..."
                ./04_ragtag_scaffolding.sh raven_assembly raven_assembly.fasta
            fi
            ;;
        wtdbg2)
            if [[ -f "./wtdbg2_assembly/pilon/round_2/pilon_round2.fasta" ]]; then
                log_info "Scaffolding wtdbg2 assembly..."
                ./04_ragtag_scaffolding.sh wtdbg2_assembly pilon/round_2/pilon_round2.fasta
            elif [[ -f "./wtdbg2_assembly/wtdbg2_assembly.fasta" ]]; then
                log_info "Scaffolding wtdbg2 assembly (no Pilon)..."
                ./04_ragtag_scaffolding.sh wtdbg2_assembly wtdbg2_assembly.fasta
            fi
            ;;
    esac
done

log_info "✓ Scaffolding complete"

#==============================================================================
# STEP 5: ASSEMBLY COMPARISON
#==============================================================================

log_section "STEP 5/5: ASSEMBLY COMPARISON"

check_script "05_assembly_comparison.sh"

log_info "Comparing all assemblies..."
./05_assembly_comparison.sh

log_info "✓ Assembly comparison complete"

#==============================================================================
# PIPELINE COMPLETE
#==============================================================================

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(( (ELAPSED % 3600) / 60 ))

log_section "PIPELINE COMPLETE"

cat << EOF
================================================================================
GENOME ASSEMBLY PIPELINE SUMMARY
================================================================================
Species: ${SPECIES}
Total time: ${HOURS}h ${MINUTES}m

Outputs:
--------
1. Filtered reads:      ./filtered_data/
2. Genome size:         ./genome_size_estimation/
3. Assemblies:          ./flye_assembly/, ./raven_assembly/, ./wtdbg2_assembly/
4. Scaffolded:          ./<assembler>/ragtag/
5. Final comparison:    ./assembly_comparison/

Next Steps:
-----------
1. Review assembly comparison: ./assembly_comparison/FINAL_COMPARISON_REPORT.txt
2. Open QUAST report:          ./assembly_comparison/quast_comparison/report.html
3. Check BUSCO scores:         ./assembly_comparison/busco_comparison/
4. Select best assembly for annotation

================================================================================
EOF

log_info "Pipeline log available in: ./pipeline_$(date +%Y%m%d).log"