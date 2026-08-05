#!/usr/bin/env bats

load 'test_helper'

@test "decontaminate::kraken2_single classifies a single-end library" {
    stub_run_tool
    run decontaminate::kraken2_single kraken2_db 18 pacbio.fastq.gz \
        report.k2report unclassified.fastq classified.fastq out.kraken2 err.log
    [ "$status" -eq 0 ]
    [[ "$(cat "$RUN_TOOL_LOG")" == "kraken2 --db kraken2_db --threads 18 --report report.k2report --unclassified-out unclassified.fastq --classified-out classified.fastq pacbio.fastq.gz" ]]
}

@test "decontaminate::kraken2_paired adds --paired and both mates" {
    stub_run_tool
    run decontaminate::kraken2_paired kraken2_db 18 r1.fastq.gz r2.fastq.gz \
        report.k2report "unclassified#.fastq" "classified#.fastq" out.kraken2 err.log
    [ "$status" -eq 0 ]
    [[ "$(cat "$RUN_TOOL_LOG")" == "kraken2 --db kraken2_db --threads 18 --paired --report report.k2report --unclassified-out unclassified#.fastq --classified-out classified#.fastq r1.fastq.gz r2.fastq.gz" ]]
}
