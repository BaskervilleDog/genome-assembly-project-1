#!/usr/bin/env bash
#==============================================================================
# HYBRID GENOME ASSEMBLY PIPELINE (Long Reads + Illumina Polishing)
#
# Purpose:
#   De novo genome assembly using PacBio long reads with three independent
#   assemblers (Flye, Raven, wtdbg2), followed by Illumina short-read
#   polishing with Pilon. Produces high-quality genome assemblies suitable
#   for publication.
#
# Biological Context:
#   - Long reads (PacBio RS2/Sequel) provide contiguity across repeats
#   - Illumina reads provide base-level accuracy for polishing
#   - Multiple assemblers mitigate algorithm-specific biases
#
# Input Requirements:
#   - PacBio subreads (FASTA/FASTQ, gzipped OK)
#   - Illumina paired-end libraries (FASTQ.gz, 1 or more library pairs)
#   - Estimated genome size (for wtdbg2 and coverage calculations)
#
# Outputs:
#   - QC reports (FastQC, MultiQC, seqkit)
#   - Filtered reads
#   - Three draft assemblies (Flye, Raven, wtdbg2)
#   - Three polished assemblies
#   - Assembly quality metrics (QUAST, BUSCO)
#   - Complete run metadata and tool versions
#
# Dependencies (conda/mamba environments):
#   seqkit, filtlong, flye, raven, wtdbg2, bwa, samtools, pilon,
#   fastqc, multiqc, fastp, quast, busco
#
# Usage:
#   ./assembly_pipeline.sh \
#     --pacbio reads.fastq.gz \
#     --illumina-r1 lib1_R1.fq.gz,lib2_R1.fq.gz \
#     --illumina-r2 lib1_R2.fq.gz,lib2_R2.fq.gz \
#     --genome-size 4.5m \
#     --sample-id SAMPLE001 \
#     --outdir results \
#     --threads 16
#
# Author: [Your Name]
# Version: 2.0
# Date: 2024-XX-XX
#==============================================================================

#set -euo pipefail

#==============================================================================
# DEFAULT CONFIGURATION
#==============================================================================

# Sample/project metadata
SAMPLE_ID="T_harzianum_TW11"
PROJECT_ID="GenomeAssembly"

# Input files
PACBIO_READS="seqfetcher/downloads/fastq/SRR10848482_subreads.fastq.gz"
ILLUMINA_R1_FILES="seqfetcher/downloads/fastq/SRR10848483_1.fastq.gz,seqfetcher/downloads/fastq/SRR10848484_1.fastq.gz"  # Comma-separated list
ILLUMINA_R2_FILES="seqfetcher/downloads/fastq/SRR10848483_2.fastq.gz,seqfetcher/downloads/fastq/SRR10848484_2.fastq.gz"  # Comma-separated list

# Genome parameters
GENOME_SIZE="40m"        # Required: e.g., "4.5m" or "40m"
PLOIDY="haploid"      # haploid or diploid

# Output directory
OUTDIR="assembly_output"

# Computational resources
THREADS=$(nproc)
JAVA_MEM="32G"

# Assembly parameters
WTDBG2_PRESET="rs"    # "rs" for RS2/Sequel, "ccs" for HiFi
PILON_ROUNDS=2

# Read filtering (FILTlong)
MIN_LENGTH=1000
TARGET_BASES=0        # 0 = no downsampling; e.g., 60000000000 for 60 Gb

# Trimming (fastp)
FASTP_QUALIFIED_QUALITY_PHRED=20
FASTP_LENGTH_REQUIRED=50

# Assembly selection (run all by default)
RUN_FLYE=true
RUN_RAVEN=true
RUN_WTDBG2=true

# QC options
RUN_BUSCO=true
BUSCO_FORCE=true
BUSCO_LINEAGE="hypocreales_odb10"  # bacteria_odb10, fungi_odb10, etc.
RUN_QUAST=true

# Scaffolding options
RUN_SCAFFOLDING=true
SCAFFOLD_REF="downloads/GCF_003025095.1/ncbi_dataset/data/GCF_003025095.1/GCF_003025095.1_Triha_v1.0_genomic.fna"

# Checkpointing (skip completed steps)
RESUME=true

# Cleanup (remove intermediate files to save space)
CLEANUP_SAMS=true
CLEANUP_BAMS=false    # Keep BAMs for validation

# Conda/mamba environment names
ENV_SEQKIT="seqkit"
ENV_FILTLONG="filtlong"
ENV_FLYE="flye"
ENV_RAVEN="raven"
ENV_WTDBG2="wtdbg2"
ENV_BWA="bwa"
ENV_SAMTOOLS="samtools"
ENV_PILON="pilon"
ENV_FASTQC="fastqc"
ENV_MULTIQC="multiqc"
ENV_FASTP="fastp"
ENV_RAGTAG="ragtag"
ENV_QUAST="quast"
ENV_BUSCO="busco"

#==============================================================================
# USAGE AND HELP
#==============================================================================

usage() {
  cat << EOF
Usage: $0 [OPTIONS]

Required Arguments:
  --pacbio FILE           PacBio subreads (FASTA/FASTQ, can be gzipped)
  --illumina-r1 FILES     Comma-separated R1 files (e.g., lib1_R1.fq.gz,lib2_R1.fq.gz)
  --illumina-r2 FILES     Comma-separated R2 files (e.g., lib1_R2.fq.gz,lib2_R2.fq.gz)
  --genome-size SIZE      Estimated genome size (e.g., 4.5m, 40m, 3.2g)
  --sample-id ID          Sample identifier (used in output filenames)

Optional Arguments:
  --outdir DIR            Output directory (default: assembly_output)
  --threads N             CPU threads (default: all available)
  --java-mem SIZE         Java heap size for Pilon (default: 32G)
  --ploidy TYPE           haploid or diploid (default: haploid)
  --pilon-rounds N        Number of polishing rounds (default: 2)
  --min-length N          Minimum read length for FILTlong (default: 1000)
  --target-bases N        Target bases for downsampling (default: 0 = disabled)
  --wtdbg2-preset PRESET  wtdbg2 preset: rs, ccs, sq, ont (default: rs)
  --skip-flye             Skip Flye assembly
  --skip-raven            Skip Raven assembly
  --skip-wtdbg2           Skip wtdbg2 assembly
  --skip-busco            Skip BUSCO gene completeness check
  --skip-quast            Skip QUAST assembly statistics
  --busco-lineage NAME    BUSCO lineage dataset (default: bacteria_odb10)
  --no-resume             Disable checkpointing (re-run all steps)
  --cleanup-all           Remove all intermediate files (SAMs and BAMs)
  --help                  Show this help message
  --scaffold-ref FILE     Reference genome FASTA for RagTag scaffolding
  --enable-scaffolding    Enable scaffolding step (default: off)
  --skip-scaffolding      Disable scaffolding step (overrides enable)

Example:
  $0 \\
    --pacbio data/pacbio_subreads.fastq.gz \\
    --illumina-r1 data/lib1_R1.fq.gz,data/lib2_R1.fq.gz \\
    --illumina-r2 data/lib1_R2.fq.gz,data/lib2_R2.fq.gz \\
    --genome-size 4.5m \\
    --sample-id EColi_K12 \\
    --threads 32 \\
    --pilon-rounds 2

EOF
  exit 1
}

#==============================================================================
# COMMAND-LINE ARGUMENT PARSING
#==============================================================================

if [[ $# -eq 0 ]]; then
  usage
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pacbio)
      PACBIO_READS="$2"
      shift 2
      ;;
    --illumina-r1)
      ILLUMINA_R1_FILES="$2"
      shift 2
      ;;
    --illumina-r2)
      ILLUMINA_R2_FILES="$2"
      shift 2
      ;;
    --genome-size)
      GENOME_SIZE="$2"
      shift 2
      ;;
    --sample-id)
      SAMPLE_ID="$2"
      shift 2
      ;;
    --outdir)
      OUTDIR="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --java-mem)
      JAVA_MEM="$2"
      shift 2
      ;;
    --ploidy)
      PLOIDY="$2"
      shift 2
      ;;
    --pilon-rounds)
      PILON_ROUNDS="$2"
      shift 2
      ;;
    --min-length)
      MIN_LENGTH="$2"
      shift 2
      ;;
    --target-bases)
      TARGET_BASES="$2"
      shift 2
      ;;
    --wtdbg2-preset)
      WTDBG2_PRESET="$2"
      shift 2
      ;;
    --skip-flye)
      RUN_FLYE=false
      shift
      ;;
    --skip-raven)
      RUN_RAVEN=false
      shift
      ;;
    --skip-wtdbg2)
      RUN_WTDBG2=false
      shift
      ;;
    --skip-busco)
      RUN_BUSCO=false
      shift
      ;;
    --skip-quast)
      RUN_QUAST=false
      shift
      ;;
    --busco-lineage)
      BUSCO_LINEAGE="$2"
      shift 2
      ;;
    --no-resume)
      RESUME=false
      shift
      ;;
    --cleanup-all)
      CLEANUP_SAMS=true
      CLEANUP_BAMS=true
      shift
      ;;
    --enable-scaffolding)
      RUN_SCAFFOLDING=true
      shift
      ;;
    --skip-scaffolding)
      RUN_SCAFFOLDING=false
      shift
      ;;
    --scaffold-ref)
      SCAFFOLD_REF="$2"
      shift 2
      ;;
    --help)
      usage
      ;;
    *)
      echo "[ERROR] Unknown option: $1" >&2
      usage
      ;;
  esac
done

#==============================================================================
# UTILITY FUNCTIONS
#==============================================================================

# Print formatted log messages with timestamps
log_info() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] [INFO] $*"
}

log_error() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2
}

log_step() {
  echo ""
  echo "========================================================================"
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP: $*"
  echo "========================================================================"
}

# Exit with error message
die() {
  log_error "$*"
  exit 1
}

# Check if a conda/mamba environment exists
check_conda_env() {
  local env_name="$1"
  if ! mamba env list | awk '{print $1}' | grep -qx "${env_name}"; then
    die "Required conda environment not found: ${env_name}"
  fi
}

# Check if file exists and is non-empty
check_file() {
  local filepath="$1"
  local description="${2:-File}"
  
  if [[ ! -f "${filepath}" ]]; then
    die "${description} not found: ${filepath}"
  fi
  
  if [[ ! -s "${filepath}" ]]; then
    die "${description} is empty: ${filepath}"
  fi
}

# Run command and capture tool version
run_with_version_log() {
  local env_name="$1"
  local tool_name="$2"
  shift 2
  
  # Log tool version (if available)
  if mamba run -n "${env_name}" "${tool_name}" --version &>/dev/null; then
    local version
    version=$(mamba run -n "${env_name}" "${tool_name}" --version 2>&1 | head -n1)
    echo "${tool_name}: ${version}" >> "${METADATA_FILE}"
  fi
  
  # Execute command
  mamba run -n "${env_name}" "$@"
}

# Check if step should be skipped (resume mode)
should_skip() {
  local output_file="$1"
  local step_name="${2:-this step}"
  
  if [[ "${RESUME}" == "true" ]] && [[ -s "${output_file}" ]]; then
    log_info "Skipping ${step_name} (output exists: ${output_file})"
    return 0
  else
    return 1
  fi
}

# Validate genome size format
validate_genome_size() {
  local size="$1"
  if [[ ! "${size}" =~ ^[0-9]+(\.[0-9]+)?[kmg]$ ]]; then
    die "Invalid genome size format: ${size}. Expected format: 4.5m, 40m, 3.2g, etc."
  fi
}

# Convert genome size to base pairs (approximate)
genome_size_to_bp() {
  local size="$1"
  local number="${size//[kmg]/}"
  local unit="${size: -1}"
  
  case "${unit}" in
    k)
      echo $(echo "${number} * 1000" | bc | cut -d. -f1)
      ;;
    m)
      echo $(echo "${number} * 1000000" | bc | cut -d. -f1)
      ;;
    g)
      echo $(echo "${number} * 1000000000" | bc | cut -d. -f1)
      ;;
    *)
      echo "${number}"
      ;;
  esac
}

#==============================================================================
# INPUT VALIDATION
#==============================================================================

log_step "Validating Inputs and Configuration"

# Check required arguments
[[ -z "${PACBIO_READS}" ]] && die "Missing required argument: --pacbio"
[[ -z "${ILLUMINA_R1_FILES}" ]] && die "Missing required argument: --illumina-r1"
[[ -z "${ILLUMINA_R2_FILES}" ]] && die "Missing required argument: --illumina-r2"
[[ -z "${GENOME_SIZE}" ]] && die "Missing required argument: --genome-size"
[[ -z "${SAMPLE_ID}" ]] && die "Missing required argument: --sample-id"

# Validate genome size format
validate_genome_size "${GENOME_SIZE}"
GENOME_SIZE_BP=$(genome_size_to_bp "${GENOME_SIZE}")
log_info "Genome size: ${GENOME_SIZE} (${GENOME_SIZE_BP} bp)"

# Validate ploidy
if [[ "${PLOIDY}" != "haploid" ]] && [[ "${PLOIDY}" != "diploid" ]]; then
  die "Invalid ploidy: ${PLOIDY}. Must be 'haploid' or 'diploid'."
fi

# Check PacBio reads exist
check_file "${PACBIO_READS}" "PacBio reads file"

# Parse and validate Illumina library files
IFS=',' read -ra R1_ARRAY <<< "${ILLUMINA_R1_FILES}"
IFS=',' read -ra R2_ARRAY <<< "${ILLUMINA_R2_FILES}"

if [[ ${#R1_ARRAY[@]} -ne ${#R2_ARRAY[@]} ]]; then
  die "Number of R1 files (${#R1_ARRAY[@]}) does not match R2 files (${#R2_ARRAY[@]})"
fi

NUM_ILLUMINA_LIBS=${#R1_ARRAY[@]}
log_info "Found ${NUM_ILLUMINA_LIBS} Illumina library pair(s)"

for i in "${!R1_ARRAY[@]}"; do
  check_file "${R1_ARRAY[$i]}" "Illumina R1 library $((i+1))"
  check_file "${R2_ARRAY[$i]}" "Illumina R2 library $((i+1))"
done

# Validate thread count
if [[ ! "${THREADS}" =~ ^[0-9]+$ ]] || [[ "${THREADS}" -lt 1 ]]; then
  die "Invalid thread count: ${THREADS}"
fi

# Check at least one assembler is enabled
if [[ "${RUN_FLYE}" == "false" ]] && [[ "${RUN_RAVEN}" == "false" ]] && [[ "${RUN_WTDBG2}" == "false" ]]; then
  die "At least one assembler must be enabled (Flye, Raven, or wtdbg2)"
fi

# Validate scaffolding settings
if [[ "${RUN_SCAFFOLDING}" == "true" ]]; then
  [[ -z "${SCAFFOLD_REF}" ]] && die "Scaffolding enabled but --scaffold-ref not provided"
  check_file "${SCAFFOLD_REF}" "Scaffolding reference genome"
fi

# Validate conda environments
log_info "Checking conda/mamba environments..."
check_conda_env "${ENV_SEQKIT}"
check_conda_env "${ENV_FILTLONG}"
[[ "${RUN_FLYE}" == "true" ]] && check_conda_env "${ENV_FLYE}"
[[ "${RUN_RAVEN}" == "true" ]] && check_conda_env "${ENV_RAVEN}"
[[ "${RUN_WTDBG2}" == "true" ]] && check_conda_env "${ENV_WTDBG2}"
check_conda_env "${ENV_BWA}"
check_conda_env "${ENV_SAMTOOLS}"
check_conda_env "${ENV_PILON}"
check_conda_env "${ENV_FASTQC}"
check_conda_env "${ENV_MULTIQC}"
check_conda_env "${ENV_FASTP}"
[[ "${RUN_QUAST}" == "true" ]] && check_conda_env "${ENV_QUAST}"
[[ "${RUN_BUSCO}" == "true" ]] && check_conda_env "${ENV_BUSCO}"
[[ "${RUN_SCAFFOLDING}" == "true" ]] && check_conda_env "${ENV_RAGTAG}"

log_info "All environment checks passed"

#==============================================================================
# DIRECTORY SETUP
#==============================================================================

log_step "Creating Output Directory Structure"

# Main directories
QC_DIR="${OUTDIR}/${SAMPLE_ID}/qc"
FILTER_DIR="${OUTDIR}/${SAMPLE_ID}/filtered_reads"
ASM_DIR="${OUTDIR}/${SAMPLE_ID}/assemblies"
POLISH_DIR="${OUTDIR}/${SAMPLE_ID}/polished_assemblies"
METRICS_DIR="${OUTDIR}/${SAMPLE_ID}/assembly_metrics"
LOG_DIR="${OUTDIR}/${SAMPLE_ID}/logs"
TRIM_DIR="${OUTDIR}/${SAMPLE_ID}/trimmed_reads"
SCAFFOLD_DIR="${OUTDIR}/${SAMPLE_ID}/scaffolded_assemblies"

# Create all subdirectories
mkdir -p \
  "${QC_DIR}/longreads/stats" \
  "${QC_DIR}/longreads/lengths" \
  "${QC_DIR}/illumina/fastqc_raw" \
  "${QC_DIR}/illumina/multiqc_raw" \
  "${QC_DIR}/illumina/fastqc_trimmed" \
  "${QC_DIR}/illumina/multiqc_trimmed" \
  "${FILTER_DIR}" \
  "${ASM_DIR}/flye" \
  "${ASM_DIR}/raven" \
  "${ASM_DIR}/wtdbg2" \
  "${POLISH_DIR}/flye" \
  "${POLISH_DIR}/raven" \
  "${POLISH_DIR}/wtdbg2" \
  "${METRICS_DIR}/quast" \
  "${METRICS_DIR}/busco" \
  "${SCAFFOLD_DIR}/flye" \
  "${SCAFFOLD_DIR}/raven" \
  "${SCAFFOLD_DIR}/wtdbg2" \
  "${LOG_DIR}/scaffolding" \
  "${LOG_DIR}/filtlong" \
  "${LOG_DIR}/flye" \
  "${LOG_DIR}/raven" \
  "${LOG_DIR}/wtdbg2" \
  "${LOG_DIR}/qc" \
  "${LOG_DIR}/trimming" \
  "${LOG_DIR}/quast" \
  "${LOG_DIR}/busco" \
  "${LOG_DIR}/polishing" \
  "${TRIM_DIR}"

QUAST_DIR="${METRICS_DIR}/quast"
BUSCO_DIR="${METRICS_DIR}/busco"

# Create metadata file
METADATA_FILE="${OUTDIR}/${SAMPLE_ID}/pipeline_metadata.txt"
cat > "${METADATA_FILE}" << EOF
Pipeline: Hybrid Genome Assembly (PacBio + Illumina)
Sample ID: ${SAMPLE_ID}
Project ID: ${PROJECT_ID}
Run Date: $(date +'%Y-%m-%d %H:%M:%S')
User: $(whoami)
Hostname: $(hostname)
Working Directory: $(pwd)

=== Input Files ===
PacBio Reads: ${PACBIO_READS}
Illumina R1 Files: ${ILLUMINA_R1_FILES}
Illumina R2 Files: ${ILLUMINA_R2_FILES}

=== Genome Parameters ===
Genome Size: ${GENOME_SIZE} (${GENOME_SIZE_BP} bp)
Ploidy: ${PLOIDY}

=== Computational Resources ===
Threads: ${THREADS}
Java Memory: ${JAVA_MEM}

=== Assembly Parameters ===
Run Flye: ${RUN_FLYE}
Run Raven: ${RUN_RAVEN}
Run wtdbg2: ${RUN_WTDBG2}
wtdbg2 Preset: ${WTDBG2_PRESET}
Pilon Polishing Rounds: ${PILON_ROUNDS}

=== Read Filtering (FILTlong) ===
Minimum Length: ${MIN_LENGTH} bp
Target Bases: ${TARGET_BASES}

=== Scaffolding (RagTag) ===
Run Scaffolding: ${RUN_SCAFFOLDING}
Scaffolding Reference: ${SCAFFOLD_REF}

=== Tool Versions ===
EOF

log_info "Directory structure created successfully"

#==============================================================================
# QUALITY CONTROL: RAW PACBIO READS
#==============================================================================

log_step "QC: Raw PacBio Reads (seqkit)"

RAW_STATS="${QC_DIR}/longreads/stats/${SAMPLE_ID}_pacbio_raw_seqkit_stats.txt"
RAW_LENGTHS="${QC_DIR}/longreads/lengths/${SAMPLE_ID}_pacbio_raw_read_lengths.tsv"

if ! should_skip "${RAW_STATS}" "raw PacBio QC"; then
  log_info "Generating summary statistics..."
  run_with_version_log "${ENV_SEQKIT}" seqkit \
    seqkit stats -a "${PACBIO_READS}" \
    > "${RAW_STATS}"
  
  log_info "Extracting read length distribution..."
  mamba run -n "${ENV_SEQKIT}" \
    seqkit fx2tab -n -l "${PACBIO_READS}" \
    > "${RAW_LENGTHS}"
  
  # Calculate and log basic statistics
  local total_reads
  total_reads=$(grep -v '^file' "${RAW_STATS}" | awk '{print $4}' | sed 's/,//g')
  local total_bases
  total_bases=$(grep -v '^file' "${RAW_STATS}" | awk '{print $5}' | sed 's/,//g')
  local mean_length
  mean_length=$(grep -v '^file' "${RAW_STATS}" | awk '{print $7}' | sed 's/,//g')
  
  log_info "Raw PacBio reads: ${total_reads} reads, ${total_bases} bp total, mean length ${mean_length} bp"
  
  # Estimate coverage
  if [[ -n "${total_bases}" ]] && [[ -n "${GENOME_SIZE_BP}" ]]; then
    local coverage
    coverage=$(echo "scale=1; ${total_bases} / ${GENOME_SIZE_BP}" | bc)
    log_info "Estimated coverage: ${coverage}x"
  fi
fi

#==============================================================================
# READ FILTERING WITH FILTLONG
#==============================================================================

log_step "Filtering PacBio Reads (FILTlong)"

FILTERED_READS="${FILTER_DIR}/${SAMPLE_ID}_pacbio_filtlong.fastq"
FILTLONG_LOG="${LOG_DIR}/filtlong/${SAMPLE_ID}_filtlong.stderr.log"

if ! should_skip "${FILTERED_READS}" "FILTlong filtering"; then
  log_info "Filtering parameters: min_length=${MIN_LENGTH} bp, target_bases=${TARGET_BASES}"
  
  # Biological rationale: Remove short reads that cannot span repetitive regions,
  # and optionally downsample to reduce computational requirements while 
  # maintaining sufficient coverage for assembly
  
  if [[ "${TARGET_BASES}" -gt 0 ]]; then
    log_info "Downsampling to ${TARGET_BASES} bases (highest quality reads)"
    mamba run -n "${ENV_FILTLONG}" \
      filtlong \
        --min_length "${MIN_LENGTH}" \
        --target_bases "${TARGET_BASES}" \
        "${PACBIO_READS}" \
      > "${FILTERED_READS}" \
      2> "${FILTLONG_LOG}"
  else
    log_info "No downsampling (keeping all reads ≥${MIN_LENGTH} bp)"
    mamba run -n "${ENV_FILTLONG}" \
      filtlong \
        --min_length "${MIN_LENGTH}" \
        "${PACBIO_READS}" \
      > "${FILTERED_READS}" \
      2> "${FILTLONG_LOG}"
  fi
  
  check_file "${FILTERED_READS}" "Filtered reads"
  log_info "Filtering complete"
fi

#==============================================================================
# QUALITY CONTROL: FILTERED PACBIO READS
#==============================================================================

log_step "QC: Filtered PacBio Reads (seqkit)"

FILTERED_STATS="${QC_DIR}/longreads/stats/${SAMPLE_ID}_pacbio_filtlong_seqkit_stats.txt"
FILTERED_LENGTHS="${QC_DIR}/longreads/lengths/${SAMPLE_ID}_pacbio_filtlong_read_lengths.tsv"

if ! should_skip "${FILTERED_STATS}" "filtered PacBio QC"; then
  log_info "Generating summary statistics for filtered reads..."
  mamba run -n "${ENV_SEQKIT}" \
    seqkit stats -a "${FILTERED_READS}" \
    > "${FILTERED_STATS}"
  
  log_info "Extracting read length distribution..."
  mamba run -n "${ENV_SEQKIT}" \
    seqkit fx2tab -n -l "${FILTERED_READS}" \
    > "${FILTERED_LENGTHS}"
  
  # Log filtering results
  local filt_reads
  filt_reads=$(grep -v '^file' "${FILTERED_STATS}" | awk '{print $4}' | sed 's/,//g')
  local filt_bases
  filt_bases=$(grep -v '^file' "${FILTERED_STATS}" | awk '{print $5}' | sed 's/,//g')
  local filt_n50
  filt_n50=$(grep -v '^file' "${FILTERED_STATS}" | awk '{print $14}' | sed 's/,//g')
  
  log_info "Filtered reads: ${filt_reads} reads, ${filt_bases} bp total, N50 ${filt_n50} bp"
  
  # Estimate post-filter coverage
  if [[ -n "${filt_bases}" ]] && [[ -n "${GENOME_SIZE_BP}" ]]; then
    local coverage
    coverage=$(echo "scale=1; ${filt_bases} / ${GENOME_SIZE_BP}" | bc)
    log_info "Post-filtering coverage: ${coverage}x"
    
    # Warn if coverage is too low for good assembly
    if (( $(echo "${coverage} < 30" | bc -l) )); then
      log_error "WARNING: Coverage below 30x may result in fragmented assembly"
    fi
  fi
  
  # Warn if N50 is suspiciously low for PacBio RS2
  if [[ -n "${filt_n50}" ]] && (( filt_n50 < 3000 )); then
    log_error "WARNING: N50 read length (${filt_n50} bp) is low for PacBio RS2. Expected >5kb."
  fi
fi

#==============================================================================
# QUALITY CONTROL: RAW ILLUMINA READS
#==============================================================================

log_step "QC: Raw Illumina Reads (FastQC + MultiQC)"

FASTQC_RAW_DIR="${QC_DIR}/illumina/fastqc_raw"
MULTIQC_RAW_DIR="${QC_DIR}/illumina/multiqc_raw"

if ! should_skip "${MULTIQC_RAW_DIR}/multiqc_report.html" "raw Illumina QC"; then
  log_info "Running FastQC on ${NUM_ILLUMINA_LIBS} library pair(s)..."
  
  # Collect all R1 and R2 files for FastQC
  local all_illumina_files=()
  for i in "${!R1_ARRAY[@]}"; do
    all_illumina_files+=("${R1_ARRAY[$i]}" "${R2_ARRAY[$i]}")
  done
  
  mamba run -n "${ENV_FASTQC}" \
    fastqc -t "${THREADS}" \
    -o "${FASTQC_RAW_DIR}" \
    "${all_illumina_files[@]}" \
    > "${LOG_DIR}/qc/${SAMPLE_ID}_fastqc_raw.stdout.log" \
    2> "${LOG_DIR}/qc/${SAMPLE_ID}_fastqc_raw.stderr.log"
  
  log_info "Aggregating FastQC reports with MultiQC..."
  mamba run -n "${ENV_MULTIQC}" \
    multiqc "${FASTQC_RAW_DIR}" \
    -o "${MULTIQC_RAW_DIR}" \
    -n "${SAMPLE_ID}_multiqc_raw" \
    --force \
    > "${LOG_DIR}/qc/${SAMPLE_ID}_multiqc_raw.stdout.log" \
    2> "${LOG_DIR}/qc/${SAMPLE_ID}_multiqc_raw.stderr.log"
  
  log_info "Raw Illumina QC complete. Report: ${MULTIQC_RAW_DIR}/${SAMPLE_ID}_multiqc_raw.html"
fi

#==============================================================================
# ILLUMINA READ TRIMMING WITH FASTP
#==============================================================================

log_step "Trimming Illumina Reads (fastp)"

# Biological rationale: Remove adapter sequences and low-quality bases
# to improve mapping accuracy during polishing. Poor-quality bases can
# introduce false corrections in Pilon.

TRIMMED_R1_ARRAY=()
TRIMMED_R2_ARRAY=()

for i in "${!R1_ARRAY[@]}"; do
  lib_num=$((i + 1))
  log_info "Trimming library ${lib_num}/${NUM_ILLUMINA_LIBS}..."
  
  lib_trim_dir="${TRIM_DIR}/lib${lib_num}"
  mkdir -p "${lib_trim_dir}"
  
  r1_trimmed="${lib_trim_dir}/${SAMPLE_ID}_lib${lib_num}_R1.trimmed.fastq.gz"
  r2_trimmed="${lib_trim_dir}/${SAMPLE_ID}_lib${lib_num}_R2.trimmed.fastq.gz"
  
  TRIMMED_R1_ARRAY+=("${r1_trimmed}")
  TRIMMED_R2_ARRAY+=("${r2_trimmed}")
  
  if should_skip "${r1_trimmed}" "trimming library ${lib_num}"; then
    continue
  fi
  
  mamba run -n "${ENV_FASTP}" \
    fastp \
      -w "${THREADS}" \
      --detect_adapter_for_pe \
      --qualified_quality_phred "${FASTP_QUALIFIED_QUALITY_PHRED}" \
      --length_required "${FASTP_LENGTH_REQUIRED}" \
      --correction \
      -i "${R1_ARRAY[$i]}" \
      -I "${R2_ARRAY[$i]}" \
      -o "${r1_trimmed}" \
      -O "${r2_trimmed}" \
      -h "${lib_trim_dir}/${SAMPLE_ID}_lib${lib_num}_fastp.html" \
      -j "${lib_trim_dir}/${SAMPLE_ID}_lib${lib_num}_fastp.json" \
    > "${LOG_DIR}/trimming/${SAMPLE_ID}_lib${lib_num}_fastp.stdout.log" \
    2> "${LOG_DIR}/trimming/${SAMPLE_ID}_lib${lib_num}_fastp.stderr.log"
  
  check_file "${r1_trimmed}" "Trimmed R1 for library ${lib_num}"
  check_file "${r2_trimmed}" "Trimmed R2 for library ${lib_num}"
done

log_info "All Illumina libraries trimmed successfully"

#==============================================================================
# QUALITY CONTROL: TRIMMED ILLUMINA READS
#==============================================================================

log_step "QC: Trimmed Illumina Reads (FastQC + MultiQC)"

FASTQC_TRIMMED_DIR="${QC_DIR}/illumina/fastqc_trimmed"
MULTIQC_TRIMMED_DIR="${QC_DIR}/illumina/multiqc_trimmed"

if ! should_skip "${MULTIQC_TRIMMED_DIR}/multiqc_report.html" "trimmed Illumina QC"; then
  log_info "Running FastQC on trimmed reads..."
  
  # Collect all trimmed files
  local all_trimmed_files=()
  for i in "${!TRIMMED_R1_ARRAY[@]}"; do
    all_trimmed_files+=("${TRIMMED_R1_ARRAY[$i]}" "${TRIMMED_R2_ARRAY[$i]}")
  done
  
  mamba run -n "${ENV_FASTQC}" \
    fastqc -t "${THREADS}" \
    -o "${FASTQC_TRIMMED_DIR}" \
    "${all_trimmed_files[@]}" \
    > "${LOG_DIR}/qc/${SAMPLE_ID}_fastqc_trimmed.stdout.log" \
    2> "${LOG_DIR}/qc/${SAMPLE_ID}_fastqc_trimmed.stderr.log"
  
  log_info "Aggregating trimmed FastQC reports with MultiQC..."
  mamba run -n "${ENV_MULTIQC}" \
    multiqc "${FASTQC_TRIMMED_DIR}" \
    -o "${MULTIQC_TRIMMED_DIR}" \
    -n "${SAMPLE_ID}_multiqc_trimmed" \
    --force \
    > "${LOG_DIR}/qc/${SAMPLE_ID}_multiqc_trimmed.stdout.log" \
    2> "${LOG_DIR}/qc/${SAMPLE_ID}_multiqc_trimmed.stderr.log"
  
  log_info "Trimmed Illumina QC complete. Report: ${MULTIQC_TRIMMED_DIR}/${SAMPLE_ID}_multiqc_trimmed.html"
fi

#==============================================================================
# DE NOVO ASSEMBLY: FLYE
#==============================================================================

if [[ "${RUN_FLYE}" == "true" ]]; then
  log_step "Assembly: Flye (PacBio Raw Mode)"
  
  # Biological rationale: Flye uses a repeat graph approach that effectively
  # handles long tandem repeats and heterozygous regions. The --pacbio-raw
  # mode is designed for PacBio RS2/Sequel data with ~13-15% error rate.
  
  FLYE_ASM="${ASM_DIR}/flye/assembly.fasta"
  
  if ! should_skip "${FLYE_ASM}" "Flye assembly"; then
    log_info "Running Flye assembler (this may take several hours)..."
    
    run_with_version_log "${ENV_FLYE}" flye \
      flye \
        --pacbio-raw "${FILTERED_READS}" \
        --genome-size "${GENOME_SIZE}" \
        --out-dir "${ASM_DIR}/flye" \
        --threads "${THREADS}" \
      > "${LOG_DIR}/flye/${SAMPLE_ID}_flye.stdout.log" \
      2> "${LOG_DIR}/flye/${SAMPLE_ID}_flye.stderr.log"
    
    check_file "${FLYE_ASM}" "Flye assembly"
    log_info "Flye assembly complete: ${FLYE_ASM}"
  fi
else
  log_info "Skipping Flye assembly (disabled)"
fi

#==============================================================================
# DE NOVO ASSEMBLY: RAVEN
#==============================================================================

if [[ "${RUN_RAVEN}" == "true" ]]; then
  log_step "Assembly: Raven"
  
  # Biological rationale: Raven is a lightweight, fast assembler using
  # overlap-layout-consensus. Often produces more contiguous assemblies
  # than Flye but may struggle with highly repetitive genomes.
  
  RAVEN_ASM="${ASM_DIR}/raven/${SAMPLE_ID}_raven.fasta"
  
  if ! should_skip "${RAVEN_ASM}" "Raven assembly"; then
    log_info "Running Raven assembler..."
    
    run_with_version_log "${ENV_RAVEN}" raven \
      raven \
        --threads "${THREADS}" \
        "${FILTERED_READS}" \
      > "${RAVEN_ASM}" \
      2> "${LOG_DIR}/raven/${SAMPLE_ID}_raven.stderr.log"
    
    check_file "${RAVEN_ASM}" "Raven assembly"
    log_info "Raven assembly complete: ${RAVEN_ASM}"
  fi
else
  log_info "Skipping Raven assembly (disabled)"
fi

#==============================================================================
# DE NOVO ASSEMBLY: WTDBG2
#==============================================================================

if [[ "${RUN_WTDBG2}" == "true" ]]; then
  log_step "Assembly: wtdbg2 + wtpoa-cns"
  
  # Biological rationale: wtdbg2 uses a fuzzy de Bruijn graph approach that
  # is extremely fast and memory-efficient. Best for large genomes or when
  # computational resources are limited. The two-stage process (layout + consensus)
  # allows fine-tuning of consensus quality.
  
  WTDBG2_PREFIX="${ASM_DIR}/wtdbg2/${SAMPLE_ID}_wtdbg2"
  WTDBG2_ASM="${WTDBG2_PREFIX}.ctg.fa"
  
  if ! should_skip "${WTDBG2_ASM}" "wtdbg2 assembly"; then
    log_info "Running wtdbg2 layout stage (preset: ${WTDBG2_PRESET})..."
    
    run_with_version_log "${ENV_WTDBG2}" wtdbg2 \
      wtdbg2 \
        -x "${WTDBG2_PRESET}" \
        -g "${GENOME_SIZE}" \
        -t "${THREADS}" \
        -i "${FILTERED_READS}" \
        -fo "${WTDBG2_PREFIX}" \
      > "${LOG_DIR}/wtdbg2/${SAMPLE_ID}_wtdbg2.stdout.log" \
      2> "${LOG_DIR}/wtdbg2/${SAMPLE_ID}_wtdbg2.stderr.log"
    
    log_info "Running wtpoa-cns consensus stage..."
    mamba run -n "${ENV_WTDBG2}" \
      wtpoa-cns \
        -t "${THREADS}" \
        -i "${WTDBG2_PREFIX}.ctg.lay.gz" \
        -fo "${WTDBG2_ASM}" \
      > "${LOG_DIR}/wtdbg2/${SAMPLE_ID}_wtpoa.stdout.log" \
      2> "${LOG_DIR}/wtdbg2/${SAMPLE_ID}_wtpoa.stderr.log"
    
    check_file "${WTDBG2_ASM}" "wtdbg2 assembly"
    log_info "wtdbg2 assembly complete: ${WTDBG2_ASM}"
  fi
else
  log_info "Skipping wtdbg2 assembly (disabled)"
fi

#==============================================================================
# ILLUMINA POLISHING FUNCTION
#==============================================================================

# Polish a single assembly with Illumina reads using Pilon
# Arguments:
#   1) Assembly name (flye, raven, wtdbg2)
#   2) Input assembly FASTA path
polish_assembly() {
  local asm_name="$1"
  local asm_input="$2"
  
  check_file "${asm_input}" "${asm_name} assembly for polishing"
  
  local asm_polish_dir="${POLISH_DIR}/${asm_name}"
  local final_polished="${asm_polish_dir}/${SAMPLE_ID}_${asm_name}.pilon_polished.fasta"
  
  if should_skip "${final_polished}" "${asm_name} polishing"; then
    return 0
  fi
  
  log_info "Starting ${PILON_ROUNDS}-round Pilon polishing for ${asm_name} assembly"
  
  # Biological rationale: PacBio assemblies have high indel error rates (~13%).
  # Pilon uses high-accuracy Illumina reads to correct these errors through
  # iterative mapping and variant calling. Typically 2 rounds reduce error
  # rate to <0.01%, suitable for gene annotation and comparative genomics.
  
  local current_fasta="${asm_input}"
  
  for round in $(seq 1 "${PILON_ROUNDS}"); do
    log_info "[${asm_name}] Polishing round ${round}/${PILON_ROUNDS}"
    
    local round_dir="${asm_polish_dir}/round_${round}"
    mkdir -p "${round_dir}"
    
    local ref_fa="${round_dir}/assembly.fasta"
    local merged_bam="${round_dir}/illumina.sorted.bam"
    local pilon_output="${SAMPLE_ID}_${asm_name}_round${round}"
    local pilon_fasta="${round_dir}/${pilon_output}.fasta"
    
    # Copy input assembly to round directory
    cp "${current_fasta}" "${ref_fa}"
    
    # Index assembly for BWA
    log_info "[${asm_name}] Round ${round}: Indexing assembly with BWA..."
    mamba run -n "${ENV_BWA}" \
      bwa index "${ref_fa}" \
      > "${LOG_DIR}/polishing/${SAMPLE_ID}_${asm_name}_round${round}_bwa_index.stdout.log" \
      2> "${LOG_DIR}/polishing/${SAMPLE_ID}_${asm_name}_round${round}_bwa_index.stderr.log"
    
    # Map all Illumina libraries and collect BAMs
    local lib_bams=()
    
    for i in "${!TRIMMED_R1_ARRAY[@]}"; do
      lib_num=$((i + 1))
      log_info "[${asm_name}] Round ${round}: Mapping library ${lib_num}/${NUM_ILLUMINA_LIBS}..."
      
      local r1="${TRIMMED_R1_ARRAY[$i]}"
      local r2="${TRIMMED_R2_ARRAY[$i]}"
      local sam="${round_dir}/lib${lib_num}.sam"
      local bam="${round_dir}/lib${lib_num}.sorted.bam"
      
      # Map with BWA-MEM (designed for 70-100bp Illumina reads)
      mamba run -n "${ENV_BWA}" \
        bwa mem -t "${THREADS}" "${ref_fa}" "${r1}" "${r2}" \
        > "${sam}" \
        2> "${LOG_DIR}/polishing/${SAMPLE_ID}_${asm_name}_round${round}_lib${lib_num}_bwa_mem.stderr.log"
      
      # Verify SAM has header
      if ! grep -q '^@SQ' "${sam}"; then
        die "[${asm_name}] BWA-MEM failed for library ${lib_num} (no SAM header)"
      fi
      
      # Sort SAM to BAM
      mamba run -n "${ENV_SAMTOOLS}" \
        samtools sort -@ "${THREADS}" -o "${bam}" "${sam}" \
        2> "${LOG_DIR}/polishing/${SAMPLE_ID}_${asm_name}_round${round}_lib${lib_num}_samtools_sort.stderr.log"
      
      # Cleanup SAM if requested
      if [[ "${CLEANUP_SAMS}" == "true" ]]; then
        rm -f "${sam}"
      fi
      
      check_file "${bam}" "Sorted BAM for library ${lib_num}"
      lib_bams+=("${bam}")
      
      # Check mapping rate for quality control
      local mapping_stats
      mapping_stats=$(mamba run -n "${ENV_SAMTOOLS}" samtools flagstat "${bam}" 2>/dev/null)
      local mapped_reads
      mapped_reads=$(echo "${mapping_stats}" | grep "mapped (" | head -n1 | awk '{print $1}')
      local total_reads
      total_reads=$(echo "${mapping_stats}" | grep "in total" | awk '{print $1}')
      
      if [[ -n "${mapped_reads}" ]] && [[ -n "${total_reads}" ]] && [[ "${total_reads}" -gt 0 ]]; then
        local mapping_rate
        mapping_rate=$(echo "scale=2; 100 * ${mapped_reads} / ${total_reads}" | bc)
        log_info "[${asm_name}] Library ${lib_num} mapping rate: ${mapping_rate}%"
        
        # Warn if mapping rate is suspiciously low
        if (( $(echo "${mapping_rate} < 50" | bc -l) )); then
          log_error "[${asm_name}] WARNING: Low mapping rate (${mapping_rate}%) for library ${lib_num}"
        fi
      fi
    done
    
    # Merge all library BAMs
    log_info "[${asm_name}] Round ${round}: Merging ${NUM_ILLUMINA_LIBS} BAM file(s)..."
    
    if [[ ${#lib_bams[@]} -eq 1 ]]; then
      # Only one library, just copy it
      cp "${lib_bams[0]}" "${merged_bam}"
    else
      # Multiple libraries, merge them
      mamba run -n "${ENV_SAMTOOLS}" \
        samtools merge -@ "${THREADS}" "${merged_bam}" "${lib_bams[@]}" \
        2> "${LOG_DIR}/polishing/${SAMPLE_ID}_${asm_name}_round${round}_samtools_merge.stderr.log"
    fi
    
    # Index merged BAM
    mamba run -n "${ENV_SAMTOOLS}" \
      samtools index "${merged_bam}" \
      2> "${LOG_DIR}/polishing/${SAMPLE_ID}_${asm_name}_round${round}_samtools_index.stderr.log"
    
    check_file "${merged_bam}" "Merged BAM for polishing"
    
    # Cleanup individual library BAMs if requested
    if [[ "${CLEANUP_BAMS}" == "true" ]]; then
      rm -f "${lib_bams[@]}"
    fi
    
    # Run Pilon to correct assembly errors
    log_info "[${asm_name}] Round ${round}: Running Pilon error correction..."
    
    # Build Pilon command with appropriate flags
    local pilon_args=(
      "--genome" "${ref_fa}"
      "--frags" "${merged_bam}"
      "--output" "${pilon_output}"
      "--outdir" "${round_dir}"
      "--changes"
      "--vcf"
      "--fix" "all"
      "--mindepth" "5"
    )
    
    # Add diploid flag if needed
    if [[ "${PLOIDY}" == "diploid" ]]; then
      pilon_args+=("--diploid")
    fi
    
    # Set Java heap size via environment variable
    # This is the correct way to pass memory settings to Pilon
    export _JAVA_OPTIONS="-Xmx${JAVA_MEM}"
    
    mamba run -n "${ENV_PILON}" \
      pilon "${pilon_args[@]}" \
      > "${LOG_DIR}/polishing/${SAMPLE_ID}_${asm_name}_round${round}_pilon.stdout.log" \
      2> "${LOG_DIR}/polishing/${SAMPLE_ID}_${asm_name}_round${round}_pilon.stderr.log"
    
    # Unset to avoid affecting other Java tools
    unset _JAVA_OPTIONS
    
    check_file "${pilon_fasta}" "Pilon output for round ${round}"
    
    # Report number of changes made
    local changes_file="${round_dir}/${pilon_output}.changes"
    if [[ -f "${changes_file}" ]]; then
      local num_changes
      num_changes=$(wc -l < "${changes_file}")
      log_info "[${asm_name}] Round ${round}: Pilon made ${num_changes} corrections"
      
      # If very few changes in this round, polishing may have converged
      if [[ "${num_changes}" -lt 10 ]] && [[ "${round}" -lt "${PILON_ROUNDS}" ]]; then
        log_info "[${asm_name}] Round ${round}: Very few changes detected - assembly may be converged"
      fi
    fi
    
    # Use Pilon output as input for next round
    current_fasta="${pilon_fasta}"
  done
  
  # Copy final polished assembly to output directory
  cp "${current_fasta}" "${final_polished}"
  log_info "[${asm_name}] Polishing complete: ${final_polished}"
}

#==============================================================================
# ILLUMINA POLISHING: ALL ASSEMBLIES
#==============================================================================

log_step "Illumina Polishing with Pilon"

# Polish each assembly that was generated
if [[ "${RUN_FLYE}" == "true" ]] && [[ -f "${FLYE_ASM}" ]]; then
  polish_assembly "flye" "${FLYE_ASM}"
fi

if [[ "${RUN_RAVEN}" == "true" ]] && [[ -f "${RAVEN_ASM}" ]]; then
  polish_assembly "raven" "${RAVEN_ASM}"
fi

if [[ "${RUN_WTDBG2}" == "true" ]] && [[ -f "${WTDBG2_ASM}" ]]; then
  polish_assembly "wtdbg2" "${WTDBG2_ASM}"
fi

log_info "All assemblies polished successfully"

#==============================================================================
# GENOME SCAFFOLDING FUNCTION (RagTag)
#==============================================================================

# Scaffold a polished assembly using a reference genome
# Arguments:
#   1) Assembly name (flye, raven, wtdbg2)
#   2) Input assembly FASTA path (polished)
scaffold_assembly() {
  local asm_name="$1"
  local asm_input="$2"

  check_file "${asm_input}" "${asm_name} polished assembly for scaffolding"
  check_file "${SCAFFOLD_REF}" "Scaffolding reference genome"

  local asm_scaffold_dir="${SCAFFOLD_DIR}/${asm_name}"
  mkdir -p "${asm_scaffold_dir}"

  # RagTag output is a directory; final FASTA is ragtag.scaffold.fasta
  local ragtag_out="${asm_scaffold_dir}/ragtag_output"
  local final_scaffolded="${ragtag_out}/ragtag.scaffold.fasta"

  if should_skip "${final_scaffolded}" "${asm_name} scaffolding"; then
    return 0
  fi

  log_info "[${asm_name}] Starting reference-guided scaffolding with RagTag"
  log_info "[${asm_name}] Reference: ${SCAFFOLD_REF}"

  # Clean output directory if present (RagTag may fail otherwise)
  rm -rf "${ragtag_out}"
  mkdir -p "${ragtag_out}"

  # RagTag syntax:
  # ragtag.py scaffold <ref.fa> <query.fa> -o outdir -t threads
  mamba run -n "${ENV_RAGTAG}" \
    ragtag.py scaffold \
      "${SCAFFOLD_REF}" \
      "${asm_input}" \
      -o "${ragtag_out}" \
      -t "${THREADS}" \
    > "${LOG_DIR}/scaffolding/${SAMPLE_ID}_${asm_name}_ragtag.stdout.log" \
    2> "${LOG_DIR}/scaffolding/${SAMPLE_ID}_${asm_name}_ragtag.stderr.log"

  check_file "${final_scaffolded}" "Scaffolded assembly for ${asm_name}"
  log_info "[${asm_name}] Scaffolding complete: ${final_scaffolded}"
}

#==============================================================================
# GENOME SCAFFOLDING: RAGTAG (OPTIONAL)
#==============================================================================

if [[ "${RUN_SCAFFOLDING}" == "true" ]]; then
  log_step "Genome Scaffolding: RagTag (Reference-guided)"

  # Scaffold each polished assembly
  if [[ "${RUN_FLYE}" == "true" ]]; then
    flye_polished="${POLISH_DIR}/flye/${SAMPLE_ID}_flye.pilon_polished.fasta"
    if [[ -f "${flye_polished}" ]]; then
      scaffold_assembly "flye" "${flye_polished}"
    fi
  fi

  if [[ "${RUN_RAVEN}" == "true" ]]; then
    raven_polished="${POLISH_DIR}/raven/${SAMPLE_ID}_raven.pilon_polished.fasta"
    if [[ -f "${raven_polished}" ]]; then
      scaffold_assembly "raven" "${raven_polished}"
    fi
  fi

  if [[ "${RUN_WTDBG2}" == "true" ]]; then
    wtdbg2_polished="${POLISH_DIR}/wtdbg2/${SAMPLE_ID}_wtdbg2.pilon_polished.fasta"
    if [[ -f "${wtdbg2_polished}" ]]; then
      scaffold_assembly "wtdbg2" "${wtdbg2_polished}"
    fi
  fi

  log_info "Scaffolding complete for all available assemblies"
else
  log_info "Skipping scaffolding (disabled)"
fi

#==============================================================================
# ASSEMBLY QUALITY ASSESSMENT: QUAST
#==============================================================================

if [[ "${RUN_QUAST}" == "true" ]]; then
  log_step "Assembly Quality Assessment: QUAST"

  assemblies_to_assess=()
  asm_labels=()

  #-------------------------------
  # Flye
  #-------------------------------
  if [[ "${RUN_FLYE}" == "true" ]]; then
    flye_polished="${POLISH_DIR}/flye/${SAMPLE_ID}_flye.pilon_polished.fasta"
    flye_scaffolded="${SCAFFOLD_DIR}/flye/ragtag_output/ragtag.scaffold.fasta"

    if [[ -f "${flye_scaffolded}" ]] && [[ -s "${flye_scaffolded}" ]]; then
      assemblies_to_assess+=("${flye_scaffolded}")
      asm_labels+=("Flye_scaffolded")
      log_info "[QUAST] Using Flye scaffolded assembly: ${flye_scaffolded}"
    elif [[ -f "${flye_polished}" ]] && [[ -s "${flye_polished}" ]]; then
      assemblies_to_assess+=("${flye_polished}")
      asm_labels+=("Flye_polished")
      log_info "[QUAST] Using Flye polished assembly: ${flye_polished}"
    else
      log_info "[QUAST] Flye output not found (skipping)"
    fi
  fi

  #-------------------------------
  # Raven
  #-------------------------------
  if [[ "${RUN_RAVEN}" == "true" ]]; then
    raven_polished="${POLISH_DIR}/raven/${SAMPLE_ID}_raven.pilon_polished.fasta"
    raven_scaffolded="${SCAFFOLD_DIR}/raven/ragtag_output/ragtag.scaffold.fasta"

    if [[ -f "${raven_scaffolded}" ]] && [[ -s "${raven_scaffolded}" ]]; then
      assemblies_to_assess+=("${raven_scaffolded}")
      asm_labels+=("Raven_scaffolded")
      log_info "[QUAST] Using Raven scaffolded assembly: ${raven_scaffolded}"
    elif [[ -f "${raven_polished}" ]] && [[ -s "${raven_polished}" ]]; then
      assemblies_to_assess+=("${raven_polished}")
      asm_labels+=("Raven_polished")
      log_info "[QUAST] Using Raven polished assembly: ${raven_polished}"
    else
      log_info "[QUAST] Raven output not found (skipping)"
    fi
  fi

  #-------------------------------
  # wtdbg2
  #-------------------------------
  if [[ "${RUN_WTDBG2}" == "true" ]]; then
    wtdbg2_polished="${POLISH_DIR}/wtdbg2/${SAMPLE_ID}_wtdbg2.pilon_polished.fasta"
    wtdbg2_scaffolded="${SCAFFOLD_DIR}/wtdbg2/ragtag_output/ragtag.scaffold.fasta"

    if [[ -f "${wtdbg2_scaffolded}" ]] && [[ -s "${wtdbg2_scaffolded}" ]]; then
      assemblies_to_assess+=("${wtdbg2_scaffolded}")
      asm_labels+=("wtdbg2_scaffolded")
      log_info "[QUAST] Using wtdbg2 scaffolded assembly: ${wtdbg2_scaffolded}"
    elif [[ -f "${wtdbg2_polished}" ]] && [[ -s "${wtdbg2_polished}" ]]; then
      assemblies_to_assess+=("${wtdbg2_polished}")
      asm_labels+=("wtdbg2_polished")
      log_info "[QUAST] Using wtdbg2 polished assembly: ${wtdbg2_polished}"
    else
      log_info "[QUAST] wtdbg2 output not found (skipping)"
    fi
  fi

  #-------------------------------
  # Run QUAST
  #-------------------------------
  if [[ ${#assemblies_to_assess[@]} -eq 0 ]]; then
    log_info "[QUAST] No assemblies found for assessment (skipping)"
  else
    quast_report="${QUAST_DIR}/report.txt"

    if should_skip "${quast_report}" "QUAST assessment"; then
      log_info "[QUAST] Skipping QUAST (report already exists)"
    else
      log_info "[QUAST] Running QUAST on ${#assemblies_to_assess[@]} assembly(ies)"
      log_info "[QUAST] Assemblies: ${asm_labels[*]}"

      mamba run -n "${ENV_QUAST}" \
        quast.py \
          --threads "${THREADS}" \
          --output-dir "${QUAST_DIR}" \
          --labels "$(IFS=,; echo "${asm_labels[*]}")" \
          "${assemblies_to_assess[@]}" \
        > "${LOG_DIR}/quast/${SAMPLE_ID}_quast.stdout.log" \
        2> "${LOG_DIR}/quast/${SAMPLE_ID}_quast.stderr.log"

      check_file "${quast_report}" "QUAST report"
      log_info "[QUAST] QUAST complete: ${quast_report}"
    fi
  fi
else
  log_info "Skipping QUAST (disabled)"
fi

#==============================================================================
# ASSEMBLY COMPLETENESS ASSESSMENT: BUSCO
#==============================================================================

if [[ "${RUN_BUSCO}" == "true" ]]; then
  log_step "Assembly Completeness Assessment: BUSCO"

  # Set BUSCO_FORCE=true to overwrite existing runs
  BUSCO_FORCE="${BUSCO_FORCE:-false}"

  mkdir -p "${LOG_DIR}/busco"
  mkdir -p "${BUSCO_DIR}"

  assemblies_to_assess=()
  asm_names=()

  #-------------------------------
  # Flye
  #-------------------------------
  if [[ "${RUN_FLYE}" == "true" ]]; then
    flye_polished="${POLISH_DIR}/flye/${SAMPLE_ID}_flye.pilon_polished.fasta"
    flye_scaffolded="${SCAFFOLD_DIR}/flye/ragtag_output/ragtag.scaffold.fasta"

    if [[ -s "${flye_scaffolded}" ]]; then
      assemblies_to_assess+=("${flye_scaffolded}")
      asm_names+=("flye_scaffolded")
      log_info "[BUSCO] Using Flye scaffolded assembly: ${flye_scaffolded}"
    elif [[ -s "${flye_polished}" ]]; then
      assemblies_to_assess+=("${flye_polished}")
      asm_names+=("flye_polished")
      log_info "[BUSCO] Using Flye polished assembly: ${flye_polished}"
    else
      log_info "[BUSCO] Flye output not found (skipping)"
    fi
  fi

  #-------------------------------
  # Raven
  #-------------------------------
  if [[ "${RUN_RAVEN}" == "true" ]]; then
    raven_polished="${POLISH_DIR}/raven/${SAMPLE_ID}_raven.pilon_polished.fasta"
    raven_scaffolded="${SCAFFOLD_DIR}/raven/ragtag_output/ragtag.scaffold.fasta"

    if [[ -s "${raven_scaffolded}" ]]; then
      assemblies_to_assess+=("${raven_scaffolded}")
      asm_names+=("raven_scaffolded")
      log_info "[BUSCO] Using Raven scaffolded assembly: ${raven_scaffolded}"
    elif [[ -s "${raven_polished}" ]]; then
      assemblies_to_assess+=("${raven_polished}")
      asm_names+=("raven_polished")
      log_info "[BUSCO] Using Raven polished assembly: ${raven_polished}"
    else
      log_info "[BUSCO] Raven output not found (skipping)"
    fi
  fi

  #-------------------------------
  # wtdbg2
  #-------------------------------
  if [[ "${RUN_WTDBG2}" == "true" ]]; then
    wtdbg2_polished="${POLISH_DIR}/wtdbg2/${SAMPLE_ID}_wtdbg2.pilon_polished.fasta"
    wtdbg2_scaffolded="${SCAFFOLD_DIR}/wtdbg2/ragtag_output/ragtag.scaffold.fasta"

    if [[ -s "${wtdbg2_scaffolded}" ]]; then
      assemblies_to_assess+=("${wtdbg2_scaffolded}")
      asm_names+=("wtdbg2_scaffolded")
      log_info "[BUSCO] Using wtdbg2 scaffolded assembly: ${wtdbg2_scaffolded}"
    elif [[ -s "${wtdbg2_polished}" ]]; then
      assemblies_to_assess+=("${wtdbg2_polished}")
      asm_names+=("wtdbg2_polished")
      log_info "[BUSCO] Using wtdbg2 polished assembly: ${wtdbg2_polished}"
    else
      log_info "[BUSCO] wtdbg2 output not found (skipping)"
    fi
  fi

  #-------------------------------
  # Run BUSCO
  #-------------------------------
  if [[ ${#assemblies_to_assess[@]} -eq 0 ]]; then
    log_info "[BUSCO] No assemblies found for assessment (skipping)"
  else

    for i in "${!assemblies_to_assess[@]}"; do
      asm="${assemblies_to_assess[$i]}"
      asm_name="${asm_names[$i]}"

      outdir="${BUSCO_DIR}/${asm_name}"
      done_sentinel="${outdir}/${asm_name}.busco.done"

      stdout_log="${LOG_DIR}/busco/${SAMPLE_ID}_busco_${asm_name}.stdout.log"
      stderr_log="${LOG_DIR}/busco/${SAMPLE_ID}_busco_${asm_name}.stderr.log"

      # Skip if already done (unless forcing)
      if [[ "${BUSCO_FORCE}" != "true" ]]; then
        if [[ -f "${done_sentinel}" ]] || ls "${outdir}"/short_summary*.txt >/dev/null 2>&1; then
          log_info "[BUSCO] Skipping BUSCO for ${asm_name} (already complete)"
          continue
        fi
      fi

      log_info "[BUSCO] Running BUSCO for: ${asm_name}"
      log_info "[BUSCO] Input: ${asm}"
      log_info "[BUSCO] Lineage: ${BUSCO_LINEAGE}"
      log_info "[BUSCO] Output folder: ${outdir}"

      # BUSCO options
      busco_force_flag=""
      if [[ "${BUSCO_FORCE}" == "true" ]]; then
        busco_force_flag="-f"
      fi

      # Run BUSCO
      mamba run -n "${ENV_BUSCO}" \
        busco \
          -i "${asm}" \
          -l "${BUSCO_LINEAGE}" \
          -o "${asm_name}" \
          -m genome \
          -c "${THREADS}" \
          --out_path "${BUSCO_DIR}" \
          ${busco_force_flag} \
        > "${stdout_log}" \
        2> "${stderr_log}"

      # Validate output (BUSCO v6 safe)
      if ls "${outdir}"/short_summary*.txt >/dev/null 2>&1; then
        touch "${done_sentinel}"
        log_info "[BUSCO] BUSCO complete for ${asm_name}"
      else
        log_error "[BUSCO] BUSCO finished but summary was not found for ${asm_name}"
        log_error "[BUSCO] Check logs:"
        log_error "         STDOUT: ${stdout_log}"
        log_error "         STDERR: ${stderr_log}"
        die "[BUSCO] BUSCO failed for ${asm_name}"
      fi

    done
  fi

else
  log_info "Skipping BUSCO (disabled)"
fi

#============================================================================
# BUSCO PLOTTING
#============================================================================
mkdir assembly_output/T_harzianum_TW11/assembly_metrics/busco/busco_plots -p

cp assembly_output/T_harzianum_TW11/assembly_metrics/busco/flye_scaffolded/short_summary.specific.hypocreales_odb10.flye_scaffolded* -t assembly_output/T_harzianum_TW11/assembly_metrics/busco/busco_plots 
cp assembly_output/T_harzianum_TW11/assembly_metrics/busco/raven_scaffolded/short_summary.specific.hypocreales_odb10.raven_scaffolded* -t assembly_output/T_harzianum_TW11/assembly_metrics/busco/busco_plots 
cp assembly_output/T_harzianum_TW11/assembly_metrics/busco/wtdbg2_scaffolded/short_summary.specific.hypocreales_odb10.wtdbg2_scaffolded* -t assembly_output/T_harzianum_TW11/assembly_metrics/busco/busco_plots 

mamba run -n busco busco --plot assembly_output/T_harzianum_TW11/assembly_metrics/busco/busco_plots

mamba run -n busco busco --plot assembly_output/T_harzianum_TW11/assembly_metrics/busco/flye_scaffolded
mamba run -n busco busco --plot assembly_output/T_harzianum_TW11/assembly_metrics/busco/raven_scaffolded
mamba run -n busco busco --plot assembly_output/T_harzianum_TW11/assembly_metrics/busco/wtdbg2_scaffolded

#==============================================================================
# PIPELINE COMPLETION
#==============================================================================

log_step "Pipeline Complete"

# Update metadata file with completion time
cat >> "${METADATA_FILE}" << EOF

=== Pipeline Completion ===
End Time: $(date +'%Y-%m-%d %H:%M:%S')
Status: SUCCESS
EOF

# Summary of outputs
log_info "=========================================="
log_info "PIPELINE SUMMARY"
log_info "=========================================="
log_info "Sample ID: ${SAMPLE_ID}"
log_info "Output Directory: ${OUTDIR}/${SAMPLE_ID}"
log_info ""
log_info "Key Outputs:"
log_info "  - Filtered reads: ${FILTERED_READS}"
log_info "  - QC reports: ${QC_DIR}/"

if [[ "${RUN_FLYE}" == "true" ]]; then
  log_info "  - Flye assembly (polished): ${POLISH_DIR}/flye/${SAMPLE_ID}_flye.pilon_polished.fasta"
fi

if [[ "${RUN_RAVEN}" == "true" ]]; then
  log_info "  - Raven assembly (polished): ${POLISH_DIR}/raven/${SAMPLE_ID}_raven.pilon_polished.fasta"
fi

if [[ "${RUN_WTDBG2}" == "true" ]]; then
  log_info "  - wtdbg2 assembly (polished): ${POLISH_DIR}/wtdbg2/${SAMPLE_ID}_wtdbg2.pilon_polished.fasta"
fi

if [[ "${RUN_SCAFFOLDING}" == "true" ]]; then
  log_info ""
  log_info "Scaffolded Outputs (RagTag):"
  log_info "  - Flye scaffolded: ${SCAFFOLD_DIR}/flye/ragtag_output/ragtag.scaffold.fasta"
  log_info "  - Raven scaffolded: ${SCAFFOLD_DIR}/raven/ragtag_output/ragtag.scaffold.fasta"
  log_info "  - wtdbg2 scaffolded: ${SCAFFOLD_DIR}/wtdbg2/ragtag_output/ragtag.scaffold.fasta"
fi

if [[ "${RUN_QUAST}" == "true" ]]; then
  log_info "  - QUAST report: ${METRICS_DIR}/quast/report.html"
fi

if [[ "${RUN_BUSCO}" == "true" ]]; then
  log_info "  - BUSCO results: ${METRICS_DIR}/busco/"
fi

log_info "  - Pipeline metadata: ${METADATA_FILE}"
log_info ""
log_info "Next Steps:"
log_info "  1. Review QC reports (MultiQC, QUAST, BUSCO)"
log_info "  2. Compare assembly metrics to select best assembly"
log_info "  3. Proceed to gene annotation or comparative genomics"
log_info "=========================================="

exit 0