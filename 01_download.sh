#!/bin/bash
#==============================================================================
# GENOME DATA DOWNLOAD
#==============================================================================
# Author : Gianlucca de Urzêda Alves
# Date   : 28/01/2026
#
# Purpose:
# - Download genome sequencing data from SRA
# - Perform quality control on Illumina and PacBio reads
#
# Tools:
# mamba | seqfetcher
#==============================================================================

set -euo pipefail

#==============================================================================
# CONFIGURATION
#==============================================================================

THREADS=12

ACCESSIONS=(
  SRR10848482
  SRR10848483
  SRR10848484
)

BASE_DIR="$(pwd)"
SEQFETCHER_DIR="${BASE_DIR}/seqfetcher"
DOWNLOAD_DIR="${SEQFETCHER_DIR}/downloads/fastq"

QC_DIR="${BASE_DIR}/qc_report"
FILTERED_DIR="${BASE_DIR}/filtered_data"
FILTERED_QC_DIR="${FILTERED_DIR}/qc_report"

#==============================================================================
# STEP 1: DOWNLOAD TOOLS AND DATA
#==============================================================================

echo "========== STEP 1: INSTALLING SEQFETCHER =========="

if [[ ! -d "${SEQFETCHER_DIR}" ]]; then
    git clone https://github.com/BaskervilleDog/seqfetcher.git
fi

cd "${SEQFETCHER_DIR}"
chmod +x install.sh
mamba run -n ncbi_tools ./install.sh

echo "========== STEP 1.1: DOWNLOADING SRA DATA =========="

printf "%s\n" "${ACCESSIONS[@]}" > accession_list.txt

mamba run -n ncbi_tools seqfetcher download \
    --sra-method ena \
    --sra-accession-file accession_list.txt

cd "${BASE_DIR}"
