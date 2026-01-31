#!/bin/bash
#==============================================================================
# GENOME DATA DOWNLOAD + QC PIPELINE
#==============================================================================
# Author : Gianlucca de Urzêda Alves
# Date   : 28/01/2026
#
# Purpose:
# - Download genome sequencing data from SRA
# - Perform quality control on Illumina and PacBio reads
#
# Tools:
# mamba | seqfetcher | sra-tools | seqkit | fastqc | multiqc
# fastp | nanoplot | longQC | filtlong
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

#==============================================================================
# STEP 2: INITIAL QC (RAW DATA)
#==============================================================================

echo "========== STEP 2: INITIAL QC =========="

mkdir -p "${QC_DIR}"

#----------------------------------------------------
# 2.1 Basic statistics (seqkit)
#----------------------------------------------------

mamba run -n seqkit seqkit stats "${DOWNLOAD_DIR}"/*.fastq.gz

#----------------------------------------------------
# 2.2 FASTQC + MultiQC
#----------------------------------------------------

mamba run -n fastqc fastqc \
    "${DOWNLOAD_DIR}"/*.fastq.gz \
    -o "${QC_DIR}" \
    --threads ${THREADS}

mamba run -n multiqc multiqc "${QC_DIR}"

#----------------------------------------------------
# 2.3 Long-read QC
#----------------------------------------------------

## NanoPlot (PacBio)
mamba run -n nanoplot NanoPlot \
    --fastq "${DOWNLOAD_DIR}"/*subreads.fastq.gz \
    -t ${THREADS} \
    -o "${QC_DIR}/nanoplot_output"

## LongQC (PacBio RS II)
#mkdir -p "${QC_DIR}/longqc_output"

#for f in "${DOWNLOAD_DIR}"/*subreads.fastq.gz; do
#    sample=$(basename "$f" .subreads.fastq.gz)

#    sudo docker run --rm \
#        -v "${DOWNLOAD_DIR}:/input" \
#        -v "${QC_DIR}/longqc_output:/output" \
#        longqc sampleqc \
#            -x pb-rs2 \
#            -p 4 \
#            -o "/output/${sample}" \
#            "/input/$(basename "$f")"
#done

#==============================================================================
# STEP 3: READ FILTERING
#==============================================================================

echo "========== STEP 3: READ FILTERING =========="

mkdir -p "${FILTERED_DIR}"

#----------------------------------------------------
# 3.1 Illumina filtering (fastp)
#----------------------------------------------------

for r1 in "${DOWNLOAD_DIR}"/*_1.fastq.gz; do
    r2="${r1/_1.fastq.gz/_2.fastq.gz}"
    sample=$(basename "$r1" _1.fastq.gz)

    mamba run -n fastp fastp \
        -w ${THREADS} \
        -i "$r1" \
        -I "$r2" \
        -o "${FILTERED_DIR}/${sample}_1.fastq.gz" \
        -O "${FILTERED_DIR}/${sample}_2.fastq.gz" \
        -q 20 \
        -l 100 \
        -D
done

#----------------------------------------------------
# 3.2 PacBio filtering (Filtlong)
#----------------------------------------------------

for f in "${DOWNLOAD_DIR}"/*subreads.fastq.gz; do
    sample=$(basename "$f" .fastq.gz)

    mamba run -n filtlong filtlong \
        --min_length 1000 \
        --keep_percent 90 \
        "$f" | gzip > "${FILTERED_DIR}/${sample}.fastq.gz"
done

#==============================================================================
# STEP 4: QC ON FILTERED DATA
#==============================================================================

echo "========== STEP 4: QC ON FILTERED DATA =========="

mkdir -p "${FILTERED_QC_DIR}"

#----------------------------------------------------
# 4.1 FASTQC + MultiQC
#----------------------------------------------------

mamba run -n fastqc fastqc \
    "${FILTERED_DIR}"/*.fastq.gz \
    -o "${FILTERED_QC_DIR}" \
    --threads ${THREADS}

mamba run -n multiqc multiqc "${FILTERED_QC_DIR}"

#----------------------------------------------------
# 4.2 Long-read QC (filtered PacBio)
#----------------------------------------------------

mamba run -n nanoplot NanoPlot \
    --fastq "${FILTERED_DIR}"/*subreads.fastq.gz \
    -t ${THREADS} \
    -o "${FILTERED_QC_DIR}/nanoplot_output"

#mkdir -p "${FILTERED_QC_DIR}/longqc_output"

#for f in "${FILTERED_DIR}"/*subreads.fastq.gz; do
#    sample=$(basename "$f" .subreads.fastq.gz)

#    sudo docker run --rm \
#        -v "${FILTERED_DIR}:/input" \
#        -v "${FILTERED_QC_DIR}/longqc_output:/output" \
#        longqc sampleqc \
#            -x pb-rs2 \
#            -p 4 \
#            -o "/output/${sample}" \
#            "/input/$(basename "$f")"
#done

#==============================================================================
# END
#==============================================================================

echo "========== PIPELINE COMPLETE =========="
echo "Raw QC reports      : ${QC_DIR}"
echo "Filtered data       : ${FILTERED_DIR}"
echo "Filtered QC reports : ${FILTERED_QC_DIR}"
