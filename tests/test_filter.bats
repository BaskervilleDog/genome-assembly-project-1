#!/usr/bin/env bats

load 'test_helper'

@test "filter::pacbio disables quality filtering and sets the length floor" {
    stub_run_tool
    run filter::pacbio in.fastq.gz out.fastq.gz 1000 report.html report.json out.log err.log
    [ "$status" -eq 0 ]
    [[ "$(cat "$RUN_TOOL_LOG")" == "fastplong -Q --length_required 1000 -i in.fastq.gz -o out.fastq.gz --html report.html --json report.json" ]]
}

@test "filter::illumina requests adapter detection, a quality floor, a length floor, and correction" {
    stub_run_tool
    run filter::illumina r1.fastq.gz r2.fastq.gz out1.fastq.gz out2.fastq.gz \
        18 30 150 report.html report.json out.log err.log
    [ "$status" -eq 0 ]
    [[ "$(cat "$RUN_TOOL_LOG")" == "fastp -w 18 --detect_adapter_for_pe --qualified_quality_phred 30 --length_required 150 --correction -i r1.fastq.gz -I r2.fastq.gz -o out1.fastq.gz -O out2.fastq.gz --html report.html --json report.json" ]]
}
