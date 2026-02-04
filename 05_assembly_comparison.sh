#!/bin/bash
#==============================================================================
# ASSEMBLY COMPARISON PIPELINE (AUTO-DISCOVERY)
#==============================================================================
# Author: Gianlucca de Urzêda Alves
# Date: 30/01/2026
# Description:
#   - Compares all available assemblies and scaffolded versions
#   - Runs QUAST comparative analysis
#   - Collects BUSCO results
#   - Generates a final decision-ready report
#==============================================================================

set -euo pipefail

#==============================================================================
# CONFIGURATION
#==============================================================================

THREADS=12
SPECIES="Trichoderma_harzianum"

OUTPUT_DIR="./assembly_comparison"
LOG_DIR="${OUTPUT_DIR}/logs"
BUSCO_DIR="${OUTPUT_DIR}/busco_comparison"

# Declare assemblies here (easy to extend)
# label | assembly_fasta | busco_dir
ASSEMBLY_TABLE=(
  "Flye|./flye_assembly/pilon_round2/pilon_round2.fasta|./flye_assembly/quality_assessment/busco"
  "Flye_Scaffolded|./scaffolded_assembly/scaffold/ragtag.scaffold.fasta|./scaffolded_assembly/quality_assessment/busco"
  "Raven|./raven_assembly/raven_assembly.fasta|./raven_assembly/quality_assessment/busco"
  "Raven_Scaffolded|./raven_assembly/ragtag/scaffold/ragtag.scaffold.fasta|./raven_assembly/ragtag/quality_assessment/busco"
  "wtdbg2|./wtdbg2_assembly/wtdbg2_assembly.fasta|./wtdbg2_assembly/quality_assessment/busco"
  "wtdbg2_Scaffolded|./wtdbg2_assembly/ragtag/scaffold/ragtag.scaffold.fasta|./wtdbg2_assembly/ragtag/quality_assessment/busco"
)

#==============================================================================
# FUNCTIONS
#==============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
}

log_warn() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARN: $*"
}

log_section() {
    echo ""
    echo "========================================================================"
    echo "$*"
    echo "========================================================================"
}

#==============================================================================
# INITIAL SETUP
#==============================================================================

log_section "ASSEMBLY COMPARISON PIPELINE - START"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${BUSCO_DIR}"

ASSEMBLIES=()
LABELS=()
BUSCO_PATHS=()

log_info "Scanning available assemblies..."

for entry in "${ASSEMBLY_TABLE[@]}"; do
    IFS="|" read -r label fasta busco <<< "${entry}"

    if [[ -f "${fasta}" ]]; then
        ASSEMBLIES+=("${fasta}")
        LABELS+=("${label}")
        BUSCO_PATHS+=("${busco}")
        log_info "✓ Found assembly: ${label}"
    else
        log_warn "Missing assembly FASTA: ${label} (${fasta})"
    fi
done

if [[ ${#ASSEMBLIES[@]} -lt 2 ]]; then
    log_warn "Only ${#ASSEMBLIES[@]} assembly(ies) found"
    log_warn "At least 2 assemblies are required for comparison"
    exit 1
fi

log_info "Total assemblies to compare: ${#ASSEMBLIES[@]}"

#==============================================================================
# STEP 1: QUAST COMPARATIVE ANALYSIS
#==============================================================================

log_section "STEP 1/3: QUAST COMPARATIVE ANALYSIS"

log_info "Running QUAST..."
log_info "Assemblies: ${LABELS[*]}"

mamba run -n quast quast \
    "${ASSEMBLIES[@]}" \
    -o "${OUTPUT_DIR}/quast_comparison" \
    --labels "$(IFS=,; echo "${LABELS[*]}")" \
    --threads "${THREADS}" \
    --fungus \
    2>&1 | tee "${LOG_DIR}/quast_comparison.log"

log_info "✓ QUAST finished"
log_info "Report: ${OUTPUT_DIR}/quast_comparison/report.html"

#==============================================================================
# STEP 2: COLLECT BUSCO RESULTS
#==============================================================================

log_section "STEP 2/3: COLLECTING BUSCO RESULTS"

for i in "${!LABELS[@]}"; do
    label="${LABELS[$i]}"
    busco_path="${BUSCO_PATHS[$i]}"

    if [[ -d "${busco_path}" ]]; then
        summary=$(find "${busco_path}" -name "short_summary*.txt" | head -1 || true)

        if [[ -n "${summary}" && -f "${summary}" ]]; then
            cp "${summary}" "${BUSCO_DIR}/${label}_busco_summary.txt"
            log_info "✓ BUSCO collected for ${label}"
        else
            log_warn "BUSCO summary not found for ${label}"
        fi
    else
        log_warn "BUSCO directory missing for ${label}"
    fi
done

#==============================================================================
# STEP 3: FINAL REPORT
#==============================================================================

log_section "STEP 3/3: GENERATING FINAL REPORT"

REPORT="${OUTPUT_DIR}/FINAL_COMPARISON_REPORT.txt"

cat > "${REPORT}" << EOF
================================================================================
GENOME ASSEMBLY COMPARISON REPORT
================================================================================
Species: ${SPECIES}
Date: $(date)

Assemblies compared:
--------------------------------------------------------------------------------
EOF

for i in "${!LABELS[@]}"; do
    echo "$((i+1)). ${LABELS[$i]} : ${ASSEMBLIES[$i]}" >> "${REPORT}"
done

cat >> "${REPORT}" << EOF

================================================================================
QUAST SUMMARY
================================================================================

EOF

if [[ -f "${OUTPUT_DIR}/quast_comparison/report.tsv" ]]; then
    head -1 "${OUTPUT_DIR}/quast_comparison/report.tsv" >> "${REPORT}"
    grep -E "^# contigs|^Total length|^N50|^L50|^GC|# N's per 100 kbp" \
        "${OUTPUT_DIR}/quast_comparison/report.tsv" >> "${REPORT}"
else
    echo "QUAST TSV report not found." >> "${REPORT}"
fi

cat >> "${REPORT}" << EOF

================================================================================
BUSCO SUMMARY
================================================================================

EOF

for label in "${LABELS[@]}"; do
    busco_file="${BUSCO_DIR}/${label}_busco_summary.txt"

    if [[ -f "${busco_file}" ]]; then
        echo "${label}:" >> "${REPORT}"
        grep "C:" "${busco_file}" >> "${REPORT}"
        echo "" >> "${REPORT}"
    else
        echo "${label}: BUSCO not available" >> "${REPORT}"
        echo "" >> "${REPORT}"
    fi
done

#==============================================================================
# AUTOMATIC RECOMMENDATION
#==============================================================================

cat >> "${REPORT}" << EOF

================================================================================
AUTOMATIC RECOMMENDATION (BUSCO-BASED)
================================================================================

EOF

best_label=""
best_busco=0

for label in "${LABELS[@]}"; do
    busco_file="${BUSCO_DIR}/${label}_busco_summary.txt"

    if [[ -f "${busco_file}" ]]; then
        value=$(grep "C:" "${busco_file}" | awk -F'[:%]' '{print $2}' | tr -d ' ' | sed 's/,/./')

        if [[ -n "${value}" ]]; then
            if (( $(echo "${value} > ${best_busco}" | bc -l) )); then
                best_busco="${value}"
                best_label="${label}"
            fi
        fi
    fi
done

if [[ -n "${best_label}" ]]; then
    cat >> "${REPORT}" << EOF
Recommended assembly:
  → ${best_label} (BUSCO Complete: ${best_busco}%)

Please confirm using:
  - QUAST N50
  - Number of contigs/scaffolds
  - Duplication levels

EOF
else
    echo "Unable to automatically recommend an assembly." >> "${REPORT}"
fi

#==============================================================================
# FINALIZATION
#==============================================================================

log_info "Final report written to: ${REPORT}"
log_info "QUAST HTML report: ${OUTPUT_DIR}/quast_comparison/report.html}"

echo ""
log_section "ASSEMBLY COMPARISON COMPLETE"
echo "Next steps:"
echo "  1) Open QUAST report"
echo "  2) Review BUSCO completeness"
echo "  3) Choose best assembly for annotation"
echo ""

exit 0
