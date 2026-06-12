#!/usr/bin/env bash
#==============================================================================
# STAGES 7–10 — EVALUATE, SCAFFOLD & REPEAT MASK
# Genome Assembly Pipeline: Trichoderma harzianum TW11
#==============================================================================
# Author : Gianlucca de Urzêda Alves
# Date   : 11-02-2026
#
# Purpose:
#   - Evaluate raw assemblies and final polished assembly
#     (BUSCO, QUAST, Merqury)
#   - Scaffold polished assembly against chromosome-level reference (RagTag)
#   - Build de novo repeat library (RepeatModeler2)
#   - Soft-mask assembly for annotation (RepeatMasker)
#
# Tools (mamba environments):
#   busco        | busco
#   quast        | quast
#   merqury      | merqury, meryl
#   ragtag       | ragtag
#   repeatmodeler| RepeatModeler, RepeatMasker
#
# Outputs: 07_assembly_evaluation/, 09_ragtag/, 10_repeat_masking/
#==============================================================================

set -euo pipefail

#==============================================================================
# CONFIGURATION
#==============================================================================

THREADS=18
PROJECT_DIR="$(pwd)"

# Key paths
ASM_DIR="${PROJECT_DIR}/06_assemblies"
POLISH_DIR="${PROJECT_DIR}/08_polishing/longreads"
EVAL_DIR="${PROJECT_DIR}/07_assembly_evaluation"
RAGTAG_DIR="${PROJECT_DIR}/09_ragtag"
MASK_DIR="${PROJECT_DIR}/10_repeat_masking"
DECON_DIR="${PROJECT_DIR}/05_decontaminated"
LOG_DIR="${PROJECT_DIR}/logs"

# Final polished assembly (output of Script 2 NextPolish round 2)
POLISHED="${POLISH_DIR}/nextpolish_work_r2/genome.nextpolish.fasta"

# Reference genome for scaffolding (chromosome-level T. harzianum)
REF="${PROJECT_DIR}/00_downloads/reference/GCA_019097725.1/ncbi_dataset/data/\
GCA_019097725.1/GCA_019097725.1_SYAU_Tha_1.0_genomic.fna"

# BUSCO lineage (Hypocreales — most specific available for Trichoderma)
BUSCO_LINEAGE="hypocreaceae_odb12"

#==============================================================================
# DIRECTORY SETUP
#==============================================================================

mkdir -p \
    "${EVAL_DIR}/raw_assembly/busco" \
    "${EVAL_DIR}/raw_assembly/quast" \
    "${EVAL_DIR}/polished_assembly/busco" \
    "${EVAL_DIR}/polished_assembly/quast_reference" \
    "${EVAL_DIR}/polished_assembly/merqury" \
    "${RAGTAG_DIR}/scaffold" \
    "${MASK_DIR}/repeatmodeler" \
    "${MASK_DIR}/repeatmasker" \
    "${LOG_DIR}"

#==============================================================================
# UTILITY FUNCTIONS
#==============================================================================

log_step() {
    echo ""
    echo "========================================================================"
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP: $*"
    echo "========================================================================"
}

log_info() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] [INFO] $*"
}

#==============================================================================
# STEP 1 — EVALUATE RAW ASSEMBLIES (BUSCO + QUAST)
#==============================================================================
# BUSCO and QUAST are applied to all three raw assembler outputs to select
# the best assembly for polishing. Two complementary metrics are used because
# each measures a different quality dimension — gene completeness vs contiguity.
#
# BUSCO lineage hypocreaceae_odb12 (4,323 BUSCOs) is specific to Hypocreales,
# the order containing Trichoderma, making it the most precise reference set.
#
# The E: metric (BUSCOs with internal stop codons) is particularly informative
# post-polishing as it reflects residual frameshift errors from indels.

log_step "Evaluate raw assemblies: BUSCO (gene completeness)"

for assembler in flye raven wtdbg2; do
    case "${assembler}" in
        flye)   asm="${ASM_DIR}/flye/assembly.fasta" ;;
        raven)  asm="${ASM_DIR}/raven/raven_assembly.fasta" ;;
        wtdbg2) asm="${ASM_DIR}/wtdbg2/wtdbg2_assembly.fasta" ;;
    esac

    log_info "BUSCO — ${assembler}"
    mamba run -n busco busco \
        -i "${asm}" -f -m genome \
        -l "${BUSCO_LINEAGE}" -c "${THREADS}" \
        -o "${EVAL_DIR}/raw_assembly/busco/${assembler}" \
        > "${LOG_DIR}/busco_raw_${assembler}.stdout.log" \
        2> "${LOG_DIR}/busco_raw_${assembler}.stderr.log"
done

log_step "Evaluate raw assemblies: QUAST (contiguity)"

# Compare all three raw assemblies in a single QUAST run
mamba run -n quast quast \
    "${ASM_DIR}/flye/assembly.fasta" \
    "${ASM_DIR}/raven/raven_assembly.fasta" \
    "${ASM_DIR}/wtdbg2/wtdbg2_assembly.fasta" \
    --labels "Flye,Raven,wtdbg2" \
    -t "${THREADS}" \
    -o "${EVAL_DIR}/raw_assembly/quast" \
    > "${LOG_DIR}/quast_raw.stdout.log" \
    2> "${LOG_DIR}/quast_raw.stderr.log"

log_info "Raw assembly evaluation complete."
log_info "Review ${EVAL_DIR}/raw_assembly/ to confirm Flye is the best assembly before proceeding."

#==============================================================================
# STEP 2 — EVALUATE POLISHED ASSEMBLY (BUSCO + QUAST + Merqury)
#==============================================================================
# Three complementary tools measure different quality dimensions:
#
#   BUSCO   — gene space completeness (E: metric = residual frameshifts)
#   QUAST   -- structural contiguity; with -r computes misassembly rates
#              and NGA50. --fragmented prevents false misassembly calls at
#              reference scaffold boundaries.
#   Merqury — reference-free base accuracy (QV). Compares k-mer spectrum
#              of Illumina reads against assembly k-mers; absent k-mers
#              indicate assembly errors. QV ≥ 40 is annotation-grade.
#              Both Illumina runs are combined (k=21) for ~60× coverage.

log_step "Evaluate polished assembly: BUSCO"

mamba run -n busco busco \
    -i "${POLISHED}" -f -m genome \
    -l "${BUSCO_LINEAGE}" -c "${THREADS}" \
    -o "${EVAL_DIR}/polished_assembly/busco/polished" \
    > "${LOG_DIR}/busco_polished.stdout.log" \
    2> "${LOG_DIR}/busco_polished.stderr.log"

log_step "Evaluate polished assembly: QUAST (with reference)"

mamba run -n quast quast \
    "${POLISHED}" \
    -r "${REF}" \
    -t "${THREADS}" \
    --fragmented \
    -o "${EVAL_DIR}/polished_assembly/quast_reference" \
    > "${LOG_DIR}/quast_reference.stdout.log" \
    2> "${LOG_DIR}/quast_reference.stderr.log"

log_step "Evaluate polished assembly: Merqury (reference-free QV)"

# Build meryl k-mer database from both Illumina runs (combined ~60× coverage)
mamba run -n merqury meryl \
    count k=21 \
    output "${EVAL_DIR}/polished_assembly/merqury/reads.meryl" \
    "${DECON_DIR}/SRR10848483_unclassified_1.fastq.gz" \
    "${DECON_DIR}/SRR10848483_unclassified_2.fastq.gz" \
    "${DECON_DIR}/SRR10848484_unclassified_1.fastq.gz" \
    "${DECON_DIR}/SRR10848484_unclassified_2.fastq.gz" \
    > "${LOG_DIR}/meryl_polished.stdout.log" \
    2> "${LOG_DIR}/meryl_polished.stderr.log"

# Run Merqury (must be run from within the Merqury output directory)
MERQURY_DIR="${EVAL_DIR}/polished_assembly/merqury"
cd "${MERQURY_DIR}"

mamba run -n merqury merqury.sh \
    reads.meryl \
    "${POLISHED}" \
    merqury_output \
    >> "${LOG_DIR}/merqury_polished.stdout.log" \
    2>> "${LOG_DIR}/merqury_polished.stderr.log"

cd "${PROJECT_DIR}"

log_info "Polished assembly evaluation complete → ${EVAL_DIR}/polished_assembly/"

#==============================================================================
# STEP 3 — REFERENCE-GUIDED SCAFFOLDING (RagTag)
#==============================================================================
# RagTag uses minimap2 internally to align contigs to a chromosome-level
# reference, ordering and orienting them into pseudochromosomes.
#
# Reference: GCA_019097725.1 (T. harzianum SYAU_Tha_1.0) — chromosome-level.
#
# -u : retain all unplaced contigs — do NOT discard them. TW11 may carry
#      strain-specific sequences absent from the reference strain.
#
# IMPORTANT: Scaffolding is performed AFTER polishing and BEFORE repeat
# masking. RepeatMasker's soft-masking (lowercase) inhibits minimap2
# alignment, so masking must always be the final step before annotation.
#
# Gaps between joined contigs are represented as 100-N runs in the FASTA
# and as W/U records in the AGP file.

log_step "Reference-guided scaffolding: RagTag"

mamba run -n ragtag ragtag.py scaffold \
    "${REF}" \
    "${POLISHED}" \
    -o "${RAGTAG_DIR}/scaffold" \
    -t "${THREADS}" \
    -u \
    > "${LOG_DIR}/ragtag_scaffold.stdout.log" \
    2> "${LOG_DIR}/ragtag_scaffold.stderr.log"

SCAFFOLDED="${RAGTAG_DIR}/scaffold/ragtag.scaffold.fasta"

log_info "Scaffolded assembly → ${SCAFFOLDED}"
log_info "Key outputs:"
log_info "  ragtag.scaffold.fasta  — pseudochromosome assembly"
log_info "  ragtag.scaffold.agp    — contig placement and orientation"
log_info "  ragtag.scaffold.stats  — placed vs unplaced contig summary"

# Evaluate scaffolded assembly with BUSCO (confirm no completeness lost)
log_step "Evaluate scaffolded assembly: BUSCO"

mamba run -n busco busco \
    -i "${SCAFFOLDED}" -f -m genome \
    -l "${BUSCO_LINEAGE}" -c "${THREADS}" \
    -o "${EVAL_DIR}/polished_assembly/busco/scaffolded" \
    > "${LOG_DIR}/busco_scaffolded.stdout.log" \
    2> "${LOG_DIR}/busco_scaffolded.stderr.log"

#==============================================================================
# STEP 4 — DE NOVO REPEAT LIBRARY (RepeatModeler2)
#==============================================================================
# Relying solely on generic Dfam databases would miss T. harzianum-specific
# TEs, LTRs, and tandem repeats. RepeatModeler2 builds a custom library
# directly from the assembly using RECON and RepeatScout, then classifies
# families by homology to known repeats.
#
# -LTRStruct : activates structural LTR retrotransposon discovery via
#              LTR_Harvest and LTR_retriever — more comprehensive than
#              sequence-based detection alone.
#
# -pa 4 : parallel jobs. Each RMBlast job uses 4 cores;
#         18 cores ÷ 4 = 4 jobs (16 cores used). Passing more jobs than
#         core capacity over-subscribes the machine.
#
# NOTE: RepeatModeler2 on ~39 Mb typically runs 12–24 hours.
# Run inside tmux or screen. Recovery: -recoverDir RM_*/ directory.

log_step "De novo repeat library: RepeatModeler2"

# Build sequence database
mamba run -n repeatmodeler BuildDatabase \
    -name "${MASK_DIR}/repeatmodeler/trichoderma_harzianum" \
    "${SCAFFOLDED}" \
    > "${LOG_DIR}/repeatmodeler_builddatabase.stdout.log" \
    2> "${LOG_DIR}/repeatmodeler_builddatabase.stderr.log"

# Run RepeatModeler2 from within the repeatmodeler directory
cd "${MASK_DIR}/repeatmodeler"

mamba run -n repeatmodeler RepeatModeler \
    -database trichoderma_harzianum \
    -pa 4 \
    -LTRStruct \
    > "${LOG_DIR}/repeatmodeler.stdout.log" \
    2> "${LOG_DIR}/repeatmodeler.stderr.log"

cd "${PROJECT_DIR}"

REPEAT_LIB="${MASK_DIR}/repeatmodeler/trichoderma_harzianum-families.fa"

log_info "Custom repeat library → ${REPEAT_LIB}"

#==============================================================================
# STEP 5 — SOFT-MASK ASSEMBLY (RepeatMasker)
#==============================================================================
# RepeatMasker applies the custom library to identify all repeat instances
# and soft-mask them (convert to lowercase).
#
# -xsmall : soft-masking (lowercase) — REQUIRED by BRAKER3 and Funannotate.
#           Hard-masking (N replacement) is NOT appropriate here; it fragments
#           gene models at repeat boundaries.
#
# -lib : custom RepeatModeler2 library; ensures strain-specific repeats
#        are masked.
#
# -gff : produce a GFF repeat annotation track for downstream visualisation.
#
# -pa 4 : 4 jobs × 4 cores each = 16 cores.
#
# CRITICAL: RepeatMasker requires an ABSOLUTE path for -dir.
# Relative paths cause it to silently write output to a temporary RM_*/
# directory instead of the specified location.

log_step "Soft-masking assembly: RepeatMasker"

mamba run -n repeatmodeler RepeatMasker \
    -lib "${REPEAT_LIB}" \
    -pa 4 \
    -xsmall \
    -gff \
    -dir "${PROJECT_DIR}/${MASK_DIR}/repeatmasker" \
    "${SCAFFOLDED}" \
    > "${LOG_DIR}/repeatmasker.stdout.log" \
    2> "${LOG_DIR}/repeatmasker.stderr.log"

MASKED_ASM="${MASK_DIR}/repeatmasker/ragtag.scaffold.fasta.masked"

log_info "Soft-masked assembly → ${MASKED_ASM}"
log_info "Key outputs:"
log_info "  ragtag.scaffold.fasta.masked  — soft-masked assembly (input for annotation)"
log_info "  ragtag.scaffold.fasta.tbl     — repeat content summary (expected 5–12%)"
log_info "  ragtag.scaffold.fasta.gff     — GFF3 repeat annotation track"
log_info "  ragtag.scaffold.fasta.out     — full tab-delimited repeat hit table"

#==============================================================================
# SUMMARY
#==============================================================================

echo ""
echo "========================================================================"
echo "PIPELINE COMPLETE"
echo "========================================================================"
echo "Assembly evaluation     : ${EVAL_DIR}/"
echo "Scaffolded assembly     : ${SCAFFOLDED}"
echo "Custom repeat library   : ${REPEAT_LIB}"
echo "Soft-masked assembly    : ${MASKED_ASM}"
echo ""
echo "The soft-masked assembly is ready for structural annotation."
echo ""
echo "Recommended next steps:"
echo "  1. BRAKER3    — RNA-seq + protein homology (Augustus/GeneMark)"
echo "  2. Funannotate — end-to-end fungal annotation pipeline"
echo ""
echo "Functional annotation:"
echo "  - InterProScan  (Pfam, PANTHER, TIGRFAM domains)"
echo "  - eggNOG-mapper (COG/KEGG/GO via orthology)"
echo "  - antiSMASH     (secondary metabolite gene clusters — fungal mode)"
echo "  - dbCAN / MEROPS (CAZymes and proteases)"
echo "========================================================================"