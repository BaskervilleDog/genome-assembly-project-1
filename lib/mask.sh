#!/usr/bin/env bash
# Stage 9, the final step before structural annotation: build a de novo
# repeat library from the assembly itself (RepeatModeler2), then
# soft-mask the assembly with it (RepeatMasker). Soft-masking (lowercase),
# not hard-masking (N replacement), is what downstream annotators
# (BRAKER3, Funannotate) require - hard-masking would fragment gene
# models at repeat boundaries.

mask::build_database() {
    local db_name="$1" assembly="$2" stdout_log="$3" stderr_log="$4"
    io::run_tool repeatmodeler BuildDatabase \
        -name "$db_name" "$assembly" \
        > "$stdout_log" 2> "$stderr_log"
}

# RepeatModeler2 must be run with its working directory as the cwd (it
# writes RM_*/ recovery directories relative to cwd - see -recoverDir).
#
# Args:
#   $1 work_dir (becomes cwd)   $2 db_name   $3 parallel_jobs
#   $4 stdout log   $5 stderr log
#
# -LTRStruct activates structural LTR retrotransposon discovery
# (LTR_Harvest + LTR_retriever) in addition to sequence-based detection.
mask::repeatmodeler() {
    local work_dir="$1" db_name="$2" parallel_jobs="$3" \
        stdout_log="$4" stderr_log="$5"
    (
        cd "$work_dir"
        io::run_tool repeatmodeler RepeatModeler \
            -database "$db_name" -pa "$parallel_jobs" -LTRStruct \
            > "$stdout_log" 2> "$stderr_log"
    )
}

# Args:
#   $1 repeat_library   $2 parallel_jobs   $3 out_dir (MUST be absolute -
#   RepeatMasker silently writes to a temporary RM_*/ dir otherwise)
#   $4 assembly   $5 stdout log   $6 stderr log
mask::repeatmasker() {
    local repeat_lib="$1" parallel_jobs="$2" out_dir="$3" assembly="$4" \
        stdout_log="$5" stderr_log="$6"
    io::run_tool repeatmodeler RepeatMasker \
        -lib "$repeat_lib" -pa "$parallel_jobs" -xsmall -gff \
        -dir "$out_dir" "$assembly" \
        > "$stdout_log" 2> "$stderr_log"
}
