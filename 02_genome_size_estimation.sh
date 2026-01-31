#!/usr/bin/env bash
#==============================================================================
# GENOME SIZE ESTIMATION — ILLUMINA READS
#==============================================================================
# Date: 28/01/2026
# Author: Gianlucca de Urzêda Alves
#
# Tools:
#   - jellyfish
#   - GenomeScope2
#
# Description:
#   Estimate genome size using k-mer spectra from Illumina PE libraries
#==============================================================================

set -Eeuo pipefail

#==============================================================================
# CONFIGURATION
#==============================================================================

THREADS=10
KMER=21
HASH_SIZE=5G
MIN_COV=10
READ_LEN=150

DATA_DIR="./filtered_data"
OUT_DIR="./genome_size_estimation"

mkdir -p "${OUT_DIR}"

# Libraries: sample_name R1 R2
LIBRARIES=(
  "SRR10848483 SRR10848483_1.fastq.gz SRR10848483_2.fastq.gz"
  "SRR10848484 SRR10848484_1.fastq.gz SRR10848484_2.fastq.gz"
)

#==============================================================================
# FUNCTIONS
#==============================================================================

log() {
  echo "[INFO] $(date '+%F %T') $*"
}

run_jellyfish() {
  local sample=$1
  local r1=$2
  local r2=$3
  
  log "Running Jellyfish for ${sample}"
  
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
    "${tmp_dir}/R2.fastq"
  
  mamba run -n jellyfish jellyfish histo \
    -t "${THREADS}" \
    "${OUT_DIR}/${sample}.jf" \
    > "${OUT_DIR}/${sample}.histo"
  
  # Clean up
  rm -rf "${tmp_dir}"
}

estimate_genome_size() {
  local sample=$1

  log "Estimating genome size (k-mer sum ≥ ${MIN_COV}) for ${sample}"

  awk -v mincov="${MIN_COV}" '
    $1 >= mincov { sum += $1 * $2 }
    END { print sum }
  ' "${OUT_DIR}/${sample}.histo"
}

run_genomescope() {
  local sample=$1

  log "Running GenomeScope2 for ${sample}"

  mamba run -n genomescope genomescope2 \
    -i "${OUT_DIR}/${sample}.histo" \
    -o "${OUT_DIR}/${sample}_genomescope" \
    -k "${KMER}" \
    -p 1 \
    -l "${READ_LEN}" \
    --fitted_hist
}

#==============================================================================
# PIPELINE
#==============================================================================

log "Starting genome size estimation pipeline"

for lib in "${LIBRARIES[@]}"; do
  read SAMPLE R1 R2 <<< "${lib}"

  log "Processing library: ${SAMPLE}"

  run_jellyfish "${SAMPLE}" "${R1}" "${R2}"

  TOTAL_KMERS=$(estimate_genome_size "${SAMPLE}")
  log "Total k-mers (≥${MIN_COV}x): ${TOTAL_KMERS}"

  run_genomescope "${SAMPLE}"

  echo "------------------------------------------------------------"
done

log "Pipeline finished"

#==============================================================================
# NOTES (MANUAL INTERPRETATION)
#==============================================================================

cat <<EOF

INTERPRETATION NOTES
===================

Genome size estimation:
  Genome size = Total k-mers / Coverage peak

Coverage peak must be identified manually from:
  - ${OUT_DIR}/*_genomescope/model.txt
  - ${OUT_DIR}/*_genomescope/summary.txt
  - k-mer spectra plots

Compare estimates across libraries for consistency.

EOF