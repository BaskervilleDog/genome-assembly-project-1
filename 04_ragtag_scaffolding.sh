#!/usr/bin/env bash
#==============================================================================
# RAGTAG REFERENCE-GUIDED SCAFFOLDING PIPELINE
#==============================================================================
# Author: Gianlucca de Urzêda Alves
# Date: 30/01/2026
# Description:
#   Reference-guided scaffolding and optional gap filling using RagTag
#==============================================================================

set -Eeuo pipefail

#==============================================================================
# CONFIGURATION
#==============================================================================

# REQUIRED: assembler directory (flye_assembly, wtdbg2_assembly, raven_assembly)
ASSEMBLY_DIR="${1:?Usage: $0 <assembly_dir>}"

# REQUIRED: input assembly fasta (relative to ASSEMBLY_DIR)
ASSEMBLY_FASTA="${2:?Usage: $0 <assembly_dir> <assembly_fasta>}"

REFERENCE_GENOME="./downloads/GCF_003025095.1/ncbi_dataset/data/GCF_003025095.1/GCF_003025095.1_Triha_v1.0_genomic.fna"

THREADS=12
SPECIES="Trichoderma_harzianum"

# Output inside assembler folder
OUTPUT_DIR="${ASSEMBLY_DIR}/ragtag"
LOG_DIR="${OUTPUT_DIR}/logs"
SCAFFOLD_DIR="${OUTPUT_DIR}/scaffold"
PATCH_DIR="${OUTPUT_DIR}/patched"
QC_DIR="${OUTPUT_DIR}/quality_assessment"

INPUT_ASSEMBLY="${ASSEMBLY_DIR}/${ASSEMBLY_FASTA}"
FINAL_ASSEMBLY="${SCAFFOLD_DIR}/ragtag.scaffold.fasta"
PATCHED_ASSEMBLY="${PATCH_DIR}/ragtag.patch.fasta"

#==============================================================================
# FUNCTIONS
#==============================================================================

log() {
    echo "[$(date '+%F %T')] $*"
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

section() {
    echo ""
    echo "========================================================================"
    echo "$*"
    echo "========================================================================"
}

#==============================================================================
# STEP 0: INITIAL CHECKS
#==============================================================================

section "RAGTAG SCAFFOLDING PIPELINE — START"

mkdir -p "${LOG_DIR}" "${SCAFFOLD_DIR}" "${PATCH_DIR}" "${QC_DIR}"

[[ -f "${INPUT_ASSEMBLY}" ]]  \
    || die "Input assembly not found: ${INPUT_ASSEMBLY}"
[[ -f "${REFERENCE_GENOME}" ]] || die "Reference genome not found: ${REFERENCE_GENOME}"

mamba run -n ragtag ragtag.py --help &> /dev/null \
    || die "RagTag not available in environment"

log "✓ Inputs and dependencies verified"
log "Input assembly   : ${INPUT_ASSEMBLY}"
log "Reference genome : ${REFERENCE_GENOME}"

#==============================================================================
# STEP 1: RAGTAG SCAFFOLDING
#==============================================================================

section "STEP 1/4 — RAGTAG SCAFFOLDING"

if [[ -f "${FINAL_ASSEMBLY}" ]]; then
    log "✓ Scaffolded assembly already exists — skipping"
else
    START=$(date +%s)

    mamba run -n ragtag ragtag.py scaffold \
        "${REFERENCE_GENOME}" \
        "${INPUT_ASSEMBLY}" \
        -o "${SCAFFOLD_DIR}" \
        -t "${THREADS}" \
        -u \
        2>&1 | tee "${LOG_DIR}/ragtag_scaffold.log"

    END=$(date +%s)
    log "✓ Scaffolding complete in $((END-START)) seconds"
fi

# Basic contig statistics
INPUT_CONTIGS=$(grep -c "^>" "${INPUT_ASSEMBLY}")
SCAFFOLD_CONTIGS=$(grep -c "^>" "${FINAL_ASSEMBLY}")

log "Input contigs     : ${INPUT_CONTIGS}"
log "Output scaffolds  : ${SCAFFOLD_CONTIGS}"
log "Contig reduction  : $((INPUT_CONTIGS - SCAFFOLD_CONTIGS))"

#==============================================================================
# STEP 2: OPTIONAL GAP FILLING (PATCH)
#==============================================================================

section "STEP 2/4 — GAP FILLING (RAGTAG PATCH)"

PATCHED_ASSEMBLY="${PATCH_DIR}/ragtag.patch.fasta"

if [[ -f "${PATCHED_ASSEMBLY}" ]]; then
    log "✓ Patched assembly already exists — skipping"
else
    START=$(date +%s)

    mamba run -n ragtag ragtag.py patch \
        "${FINAL_ASSEMBLY}" \
        "${INPUT_ASSEMBLY}" \
        -o "${PATCH_DIR}" \
        -t "${THREADS}" \
        2>&1 | tee "${LOG_DIR}/ragtag_patch.log"

    END=$(date +%s)
    log "✓ Gap filling complete in $((END-START)) seconds"
fi

#==============================================================================
# STEP 3: QUALITY ASSESSMENT
#==============================================================================

section "STEP 3/4 — QUALITY ASSESSMENT"

# Assembly-stats comparison
if command -v assembly-stats &> /dev/null; then
    {
        echo "INPUT ASSEMBLY"
        assembly-stats "${INPUT_ASSEMBLY}"
        echo ""
        echo "SCAFFOLDED ASSEMBLY"
        assembly-stats "${FINAL_ASSEMBLY}"
    } > "${QC_DIR}/assembly_stats_comparison.txt"
fi

# QUAST (comparative)
if [[ ! -d "${QC_DIR}/quast" ]]; then
    mamba run -n quast quast \
        "${INPUT_ASSEMBLY}" \
        "${FINAL_ASSEMBLY}" \
        -o "${QC_DIR}/quast" \
        --labels "Input,Scaffolded" \
        -r "${REFERENCE_GENOME}" \
        --threads "${THREADS}" \
        --fungus \
        2>&1 | tee "${LOG_DIR}/quast.log"
fi

# BUSCO
if [[ ! -d "${QC_DIR}/busco" ]]; then
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
# STEP 4: SUMMARY REPORT
#==============================================================================

section "STEP 4/4 — SUMMARY REPORT"

REPORT="${OUTPUT_DIR}/RAGTAG_SCAFFOLDING_REPORT.txt"

cat > "${REPORT}" <<EOF
================================================================================
RAGTAG SCAFFOLDING REPORT
================================================================================
Species: ${SPECIES}
Date: $(date)

INPUT:
------
Assembly: ${INPUT_ASSEMBLY}
Reference: ${REFERENCE_GENOME}

RESULTS:
--------
Input contigs: ${INPUT_CONTIGS}
Scaffolds:     ${SCAFFOLD_CONTIGS}
Reduction:     $((INPUT_CONTIGS - SCAFFOLD_CONTIGS))

OUTPUT FILES:
-------------
Scaffolded assembly : ${FINAL_ASSEMBLY}
Patched assembly    : ${PATCHED_ASSEMBLY}
QUAST report        : ${QC_DIR}/quast/report.html
BUSCO results       : ${QC_DIR}/busco/

NOTES:
------
- ragtag.scaffold.fasta is the PRIMARY OUTPUT
- ragtag.patch.fasta is optional and may introduce local misjoins
- BUSCO completeness should remain stable
- Small increases in Ns are expected

================================================================================
EOF

log "✓ Report written to ${REPORT}"

#==============================================================================
# FINAL
#==============================================================================

section "RAGTAG SCAFFOLDING PIPELINE COMPLETE"

log "Final assembly : ${FINAL_ASSEMBLY}"
log "Reports        : ${QC_DIR}"
log "Done at        : $(date)"
