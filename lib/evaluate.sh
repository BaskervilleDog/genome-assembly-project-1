#!/usr/bin/env bash
# Stage 7, applied at three checkpoints (raw assemblies, polished
# assembly, scaffolded assembly): BUSCO measures gene-space completeness,
# QUAST measures structural contiguity, Merqury measures reference-free
# base accuracy (polished assembly only - it needs the short-read k-mer
# spectrum, not a reference genome).

evaluate::busco() {
    local assembly="$1" lineage="$2" threads="$3" out_dir="$4" \
        stdout_log="$5" stderr_log="$6"
    io::run_tool busco busco \
        -i "$assembly" -f -m genome \
        -l "$lineage" -c "$threads" \
        -o "$out_dir" \
        > "$stdout_log" 2> "$stderr_log"
}

# Compare multiple assemblies (e.g. the three raw assembler outputs) in a
# single QUAST run.
#
# Args:
#   $1 threads   $2 out_dir   $3 stdout log   $4 stderr log
#   $5 labels (comma-separated, matching the assembly order)
#   $6.. assembly FASTA files
evaluate::quast() {
    local threads="$1" out_dir="$2" stdout_log="$3" stderr_log="$4" labels="$5"
    shift 5
    io::run_tool quast quast "$@" \
        --labels "$labels" \
        -t "$threads" \
        -o "$out_dir" \
        > "$stdout_log" 2> "$stderr_log"
}

# QUAST against a reference genome (used for the polished assembly).
# --fragmented avoids false misassembly calls at reference scaffold
# boundaries.
evaluate::quast_reference() {
    local assembly="$1" reference="$2" threads="$3" out_dir="$4" \
        stdout_log="$5" stderr_log="$6"
    io::run_tool quast quast "$assembly" \
        -r "$reference" \
        -t "$threads" \
        --fragmented \
        -o "$out_dir" \
        > "$stdout_log" 2> "$stderr_log"
}

evaluate::meryl_count() {
    local kmer="$1" out_db="$2" stdout_log="$3" stderr_log="$4"
    shift 4
    io::run_tool merqury meryl count k="$kmer" output "$out_db" "$@" \
        > "$stdout_log" 2> "$stderr_log"
}

# Merqury must be run with its own output directory as the cwd.
#
# Args:
#   $1 run_dir (becomes cwd)   $2 reads_meryl db (relative to $1, or absolute)
#   $3 assembly fasta   $4 out_prefix   $5 stdout log   $6 stderr log
evaluate::merqury() {
    local run_dir="$1" reads_meryl="$2" assembly="$3" out_prefix="$4" \
        stdout_log="$5" stderr_log="$6"
    (
        cd "$run_dir"
        io::run_tool merqury merqury.sh "$reads_meryl" "$assembly" "$out_prefix" \
            >> "$stdout_log" 2>> "$stderr_log"
    )
}
