#!/usr/bin/env bash
# Stage 3: read filtering. Long and short reads need different tools -
# fastplong is long-read-aware; standard short-read trimmers (fastp,
# Trimmomatic) are not appropriate for PacBio reads.

# Filter PacBio long reads. Quality filtering is disabled (-Q) because
# quality scores are absent from SRA-archived PacBio data; only a length
# floor is applied (reads too short to help contiguity are dropped).
#
# Args:
#   $1 in.fastq.gz   $2 out.fastq.gz   $3 length_required
#   $4 html report   $5 json report    $6 stdout log   $7 stderr log
filter::pacbio() {
    local in="$1" out="$2" length_required="$3" html="$4" json="$5" \
        stdout_log="$6" stderr_log="$7"
    io::run_tool fastplong fastplong \
        -Q \
        --length_required "$length_required" \
        -i "$in" -o "$out" \
        --html "$html" --json "$json" \
        > "$stdout_log" 2> "$stderr_log"
}

# Filter one Illumina paired-end run: adapter trimming, a quality floor,
# a length floor, and overlap-based error correction in a single pass.
#
# Args:
#   $1 in_r1   $2 in_r2   $3 out_r1   $4 out_r2   $5 threads
#   $6 quality_phred   $7 length_required   $8 html   $9 json
#   $10 stdout log   $11 stderr log
filter::illumina() {
    local in1="$1" in2="$2" out1="$3" out2="$4" threads="$5" \
        quality_phred="$6" length_required="$7" html="$8" json="$9"
    local stdout_log="${10}" stderr_log="${11}"
    io::run_tool fastp fastp \
        -w "$threads" \
        --detect_adapter_for_pe \
        --qualified_quality_phred "$quality_phred" \
        --length_required "$length_required" \
        --correction \
        -i "$in1" -I "$in2" \
        -o "$out1" -O "$out2" \
        --html "$html" --json "$json" \
        > "$stdout_log" 2> "$stderr_log"
}
