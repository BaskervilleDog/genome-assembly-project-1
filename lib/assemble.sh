#!/usr/bin/env bash
# Stage 5: three structurally different assemblers run on the same
# decontaminated PacBio reads, so the best result is picked empirically
# (BUSCO/QUAST, in lib/evaluate.sh) rather than assumed from theory.

# Flye - repeat-graph assembler; --pacbio-corr expects already-corrected
# PacBio reads.
assemble::flye() {
    local reads="$1" genome_size="$2" out_dir="$3" threads="$4" \
        stdout_log="$5" stderr_log="$6"
    io::run_tool flye292 flye \
        --pacbio-corr "$reads" \
        --genome-size "$genome_size" \
        --out-dir "$out_dir" \
        --threads "$threads" \
        > "$stdout_log" 2> "$stderr_log"
}

# Raven - string-graph/OLC assembler with internal Racon polishing.
assemble::raven() {
    local reads="$1" threads="$2" out_fasta="$3" stderr_log="$4"
    io::run_tool raven raven --threads "$threads" "$reads" \
        > "$out_fasta" 2> "$stderr_log"
}

# wtdbg2 - fuzzy de Bruijn graph assembler. Two steps: graph construction,
# then a separate consensus-calling pass (wtpoa-cns). -x rs selects the
# PacBio RS2/Sequel error-tolerance preset.
#
# Args:
#   $1 reads   $2 genome_size   $3 threads   $4 out_prefix (passed to -fo)
#   $5 out_fasta (consensus output)   $6 stdout log   $7 stderr log
assemble::wtdbg2() {
    local reads="$1" genome_size="$2" threads="$3" out_prefix="$4" \
        out_fasta="$5" stdout_log="$6" stderr_log="$7"

    io::run_tool wtdbg2 wtdbg2 \
        -x rs -g "$genome_size" -t "$threads" \
        -i "$reads" -fo "$out_prefix" \
        > "$stdout_log" 2> "$stderr_log"

    io::run_tool wtdbg2 wtpoa-cns \
        -t "$threads" \
        -i "${out_prefix}.ctg.lay.gz" \
        -fo "$out_fasta" \
        >> "$stdout_log" 2>> "$stderr_log"
}
