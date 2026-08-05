#!/usr/bin/env bash
# Stages 2 and 3 (raw reads, then re-run on filtered reads): read-level
# QC. PacBio and Illumina use different tools because SRA-archived PacBio
# reads carry no quality scores (stripped during archiving) - length
# distribution, not per-base quality, is the meaningful signal there.

# Fast per-file summary stats (incl. N50, Q-scores) for any FASTQ library.
qc::seqkit_stats() {
    local fastq="$1" out="$2" stderr_log="$3"
    io::run_tool seqkit seqkit stats -a "$fastq" > "$out" 2> "$stderr_log"
}

# Read-length histogram + N50 distribution for PacBio reads. FastQC's
# per-base quality plots aren't meaningful here (see module note above),
# so this is the primary long-read QC signal.
qc::longreadsum() {
    local fastq="$1" out_dir="$2" stdout_log="$3" stderr_log="$4"
    io::run_tool longreadsum longreadsum fq \
        -i "$fastq" -o "$out_dir" \
        > "$stdout_log" 2> "$stderr_log"
}

# FastQC on one or more Illumina FASTQ files.
#
# Args:
#   $1   - output directory
#   $2   - threads
#   $3   - stdout log path
#   $4   - stderr log path
#   $5.. - FASTQ file(s)
qc::fastqc() {
    local out_dir="$1" threads="$2" stdout_log="$3" stderr_log="$4"
    shift 4
    io::run_tool fastqc fastqc -t "$threads" -o "$out_dir" "$@" \
        > "$stdout_log" 2> "$stderr_log"
}

# Aggregate FastQC/fastp reports from $in_dir into one interactive report.
qc::multiqc() {
    local in_dir="$1" out_dir="$2" stdout_log="$3" stderr_log="$4"
    io::run_tool multiqc multiqc "$in_dir" -o "$out_dir" \
        > "$stdout_log" 2> "$stderr_log"
}
