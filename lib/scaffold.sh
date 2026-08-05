#!/usr/bin/env bash
# Stage 8: reference-guided scaffolding. Must run after polishing and
# before repeat masking (lib/mask.sh) - RagTag's internal minimap2
# alignment step is inhibited by soft-masked (lowercase) bases.

# Args:
#   $1 reference   $2 assembly   $3 out_dir   $4 threads
#   $5 stdout log   $6 stderr log
#
# -u retains unplaced contigs rather than discarding them, since the
# assembled strain may carry sequence absent from the reference used for
# scaffolding.
scaffold::ragtag() {
    local reference="$1" assembly="$2" out_dir="$3" threads="$4" \
        stdout_log="$5" stderr_log="$6"
    io::run_tool ragtag ragtag.py scaffold \
        "$reference" "$assembly" \
        -o "$out_dir" -t "$threads" -u \
        > "$stdout_log" 2> "$stderr_log"
}
