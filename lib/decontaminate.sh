#!/usr/bin/env bash
# Stage 4: taxonomic decontamination with Kraken2. Soil-derived fungal
# cultures commonly carry bacterial/human DNA from culture media or
# handling; assembling without removing it risks chimeric contigs and
# inflated genome-size estimates. Reads marked UNCLASSIFIED (no database
# hit) are the clean dataset; classified reads are kept separately for
# inspection, not discarded outright.

decontaminate::kraken2_single() {
    local db="$1" threads="$2" fastq="$3" report="$4" \
        unclassified_out="$5" classified_out="$6" kraken_out="$7" stderr_log="$8"
    io::run_tool kraken2 kraken2 \
        --db "$db" --threads "$threads" \
        --report "$report" \
        --unclassified-out "$unclassified_out" \
        --classified-out "$classified_out" \
        "$fastq" \
        > "$kraken_out" 2> "$stderr_log"
}

decontaminate::kraken2_paired() {
    local db="$1" threads="$2" fastq1="$3" fastq2="$4" report="$5" \
        unclassified_out="$6" classified_out="$7" kraken_out="$8" stderr_log="$9"
    io::run_tool kraken2 kraken2 \
        --db "$db" --threads "$threads" --paired \
        --report "$report" \
        --unclassified-out "$unclassified_out" \
        --classified-out "$classified_out" \
        "$fastq1" "$fastq2" \
        > "$kraken_out" 2> "$stderr_log"
}
