#!/usr/bin/env bash

# Comparing mapping of sequencing reads to assembled genome

OUTDIR="assembly_output"
SAMPLE_ID="T_harzianum_TW11"

PACBIO_READS="seqfetcher/downloads/fastq/SRR10848482_subreads.fastq.gz"
ILLUMINA_R1_FILES="seqfetcher/downloads/fastq/SRR10848483_1.fastq.gz,seqfetcher/downloads/fastq/SRR10848484_1.fastq.gz"  # Comma-separated list
ILLUMINA_R2_FILES="seqfetcher/downloads/fastq/SRR10848483_2.fastq.gz,seqfetcher/downloads/fastq/SRR10848484_2.fastq.gz"  # Comma-separated list

#==============================================================================
# ENVIRONMENT DEFINITION
#==============================================================================

ENV_BWA="bwa"
ENV_SAMTOOLS="samtools"

#==============================================================================
# ASSEMBLY DEFINITION
#==============================================================================

SCAFFOLD_DIR="${OUTDIR}/${SAMPLE_ID}/scaffolded_assemblies"

ASN_FLYE="${SCAFFOLD_DIR}/flye/ragtag_output/ragtag.scaffold.fasta"
ASN_RAVEN="${SCAFFOLD_DIR}/raven/ragtag_output/ragtag.scaffold.fasta"
ASN_WTDBG2="${SCAFFOLD_DIR}/wtdbg2/ragtag_output/ragtag.scaffold.fasta"

#==============================================================================
# FOLDER STRUCTURE
#==============================================================================

LOG_DIR="${OUTDIR}/${SAMPLE_ID}/logs"
METRICS_DIR="${OUTDIR}/${SAMPLE_ID}/assembly_metrics"

mkdir -p \
    "${LOG_DIR}/bwa" \
    "${METRICS_DIR}/alignment"

#==============================================================================
# UTILITY FUNCTIONS
#==============================================================================

log_step() {
  echo ""
  echo "========================================================================"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP: $*"
  echo "========================================================================"
}

#==============================================================================
# INDEXING ASSEMBLIES
#==============================================================================

# Select the best assembly using BUSCO and QUAST

log_step "Indexing Assemblies"

# Flye assembly
mamba run -n ${ENV_BWA} bwa \
    index \
    ${ASN_FLYE} \
    > "${LOG_DIR}/bwa/${SAMPLE_ID}_flye_index.stdout.log" \
    2> "${LOG_DIR}/bwa/${SAMPLE_ID}_flye_index.stderr.log"

# Raven assembly
#mamba run -n ${ENV_BWA} bwa \
#    index \
#    ${ASN_RAVEN} \
#    > "${LOG_DIR}/bwa/${SAMPLE_ID}_raven_index.stdout.log" \
#    2> "${LOG_DIR}/bwa/${SAMPLE_ID}_raven_index.stderr.log"

# wtdbg2 assembly
#mamba run -n ${ENV_BWA} bwa \
#    index \
#    ${ASN_WTDBG2} \
#    > "${LOG_DIR}/bwa/${SAMPLE_ID}_wtdbg2_index.stdout.log" \
#    2> "${LOG_DIR}/bwa/${SAMPLE_ID}_wtdbg2_index.stderr.log"

#==============================================================================
# MAPPING SEQUENCING READS TO ASSEMBLIES
#==============================================================================

# Bwa stdout generates a SAM file

mamba run -n ${ENV_BWA} bwa \
    mem \
    -t 12 \
    -x pacbio \
    ${ASN_FLYE} \
    ${PACBIO_READS} \
    > "${METRICS_DIR}/alignment/${SAMPLE_ID}_flye_alignment_pacbio.sam" \
    2> "${LOG_DIR}/bwa/${SAMPLE_ID}_flye_alignment_pacbio.stderr.log"

# Convert to BAM file

mamba run -n ${ENV_SAMTOOLS} samtools \
    view \
    --threads 12 \
    -b \
    "${METRICS_DIR}/alignment/${SAMPLE_ID}_flye_alignment_pacbio.sam" \
    -o "${METRICS_DIR}/alignment/${SAMPLE_ID}_flye_alignment_pacbio.bam"

# Sort the BAM file generated

mamba run -n ${ENV_SAMTOOLS} samtools \
    sort \
    --threads 12 \
    "${METRICS_DIR}/alignment/${SAMPLE_ID}_flye_alignment_pacbio.bam" \
    -o "${METRICS_DIR}/alignment/${SAMPLE_ID}_flye_alignment_pacbio_sorted.bam"

mamba run -n bamtools bamtools stats -in "${METRICS_DIR}/alignment/${SAMPLE_ID}_flye_alignment_pacbio_sorted.bam"

#**********************************************
#Stats for BAM file(s): 
#**********************************************

#Total reads:       908895
#Mapped reads:      884450   (97.3105%)
#Forward strand:    469966   (51.7074%)
#Reverse strand:    438929   (48.2926%)
#Failed QC:         0    (0%)
#Duplicates:        0    (0%)
#Paired-end reads:  0    (0%)

mamba run -n qualimap qualimap \
    bamqc \
    -bam "${METRICS_DIR}/alignment/${SAMPLE_ID}_flye_alignment_pacbio_sorted.bam" \
    -outdir assembly_output/T_harzianum_TW11/assembly_metrics/alignment \
    --java-mem-size=16G