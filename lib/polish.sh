#!/usr/bin/env bash
# Stage 6: long-read structural correction first (Racon), then short-read
# base-level correction (Polypolish, then NextPolish) - structural errors
# left uncorrected would otherwise confuse the short-read alignments the
# later steps depend on.

# One minimap2 + Racon round: align $reads to $assembly_in, then produce
# a corrected $assembly_out. Call twice (second round's assembly_in is
# the first round's assembly_out) - a third round yields diminishing
# returns and isn't applied.
#
# Args:
#   $1 assembly_in   $2 reads   $3 minimap2 preset (e.g. map-pb)
#   $4 threads   $5 paf_out   $6 assembly_out
#   $7 minimap2 stderr log   $8 racon stderr log
polish::racon_round() {
    local assembly_in="$1" reads="$2" preset="$3" threads="$4" \
        paf_out="$5" assembly_out="$6" minimap_log="$7" racon_log="$8"

    io::run_tool minimap2 minimap2 \
        -t "$threads" -x "$preset" \
        "$assembly_in" "$reads" \
        > "$paf_out" 2> "$minimap_log"

    io::run_tool racon1420 racon \
        -t "$threads" \
        "$reads" "$paf_out" "$assembly_in" \
        > "$assembly_out" 2> "$racon_log"
}

# Polypolish: align every short-read library with `bwa mem -a` (ALL
# alignments, not just uniquely-mapping ones, so repetitive regions still
# get polished - this is what Pilon-style unique-alignment polishing
# skips), filter each library, then polish jointly across all of them.
# Cleans up the SAM/BWA-index intermediates it creates.
#
# Args:
#   $1 assembly   $2 threads   $3 work_dir   $4 out_fasta   $5 log_dir
#   $6.. triples of "read1 read2 label" - one triple per short-read library
polish::polypolish() {
    local assembly="$1" threads="$2" work_dir="$3" out_fasta="$4" log_dir="$5"
    shift 5

    io::run_tool bwa bwa index "$assembly"

    local sam_args=()
    while [[ $# -gt 0 ]]; do
        local r1="$1" r2="$2" label="$3"
        shift 3

        io::run_tool bwa bwa mem -t "$threads" -a "$assembly" "$r1" \
            > "$work_dir/alignments_${label}_1.sam" \
            2> "$log_dir/bwa_${label}_1.stderr.log"
        io::run_tool bwa bwa mem -t "$threads" -a "$assembly" "$r2" \
            > "$work_dir/alignments_${label}_2.sam" \
            2> "$log_dir/bwa_${label}_2.stderr.log"

        io::run_tool polypolish polypolish filter \
            --in1  "$work_dir/alignments_${label}_1.sam" \
            --in2  "$work_dir/alignments_${label}_2.sam" \
            --out1 "$work_dir/filtered_${label}_1.sam" \
            --out2 "$work_dir/filtered_${label}_2.sam"

        sam_args+=("$work_dir/filtered_${label}_1.sam" "$work_dir/filtered_${label}_2.sam")
    done

    io::run_tool polypolish polypolish polish "$assembly" "${sam_args[@]}" \
        > "$out_fasta" \
        2> "$log_dir/polypolish.stderr.log"

    rm -f "$work_dir"/alignments_*.sam "$work_dir"/filtered_*.sam \
        "$assembly".{amb,ann,bwt,pac,sa}
}

# One NextPolish round: writes its config file, then runs nextPolish.
# Call twice (round 2's $genome is round 1's
# "$work_dir/genome.nextpolish.fasta") - two rounds were used here given
# the initial QV/BUSCO internal-stop-codon rate indicated residual
# polishing potential after one.
#
# Args:
#   $1 genome fasta   $2 sgs_fofn (paired short-read file-of-filenames)
#   $3 work_dir   $4 cfg_path   $5 stdout log   $6 stderr log
polish::nextpolish_round() {
    local genome="$1" sgs_fofn="$2" work_dir="$3" cfg_path="$4" \
        stdout_log="$5" stderr_log="$6"

    cat > "$cfg_path" << EOF
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 6
multithread_jobs = 3
genome = ${genome}
genome_size = auto
workdir = ${work_dir}
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ${sgs_fofn}
sgs_options = -max_depth 100 -bwa
EOF

    io::run_tool nextpolish39 nextPolish "$cfg_path" \
        > "$stdout_log" 2> "$stderr_log"
}
