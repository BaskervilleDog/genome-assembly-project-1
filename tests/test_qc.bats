#!/usr/bin/env bats

load 'test_helper'

@test "qc::seqkit_stats runs seqkit stats -a and writes stdout to the output path" {
    stub_run_tool
    local out="$BATS_TEST_TMPDIR/stats.txt" err="$BATS_TEST_TMPDIR/err.log"

    run qc::seqkit_stats "$FIXTURES/sample_reads.fastq" "$out" "$err"
    [ "$status" -eq 0 ]
    [[ "$(cat "$RUN_TOOL_LOG")" == "seqkit stats -a $FIXTURES/sample_reads.fastq" ]]
    [[ "$(cat "$out")" == *"stub-stdout: seqkit stats -a"* ]]
}

@test "qc::longreadsum passes -i/-o to longreadsum fq" {
    stub_run_tool
    run qc::longreadsum "$FIXTURES/sample_reads.fastq" "$BATS_TEST_TMPDIR/qc_out" \
        "$BATS_TEST_TMPDIR/out.log" "$BATS_TEST_TMPDIR/err.log"
    [ "$status" -eq 0 ]
    [[ "$(cat "$RUN_TOOL_LOG")" == "longreadsum fq -i $FIXTURES/sample_reads.fastq -o $BATS_TEST_TMPDIR/qc_out" ]]
}

@test "qc::fastqc forwards every FASTQ file given after the fixed arguments" {
    stub_run_tool
    run qc::fastqc "$BATS_TEST_TMPDIR/out" 18 \
        "$BATS_TEST_TMPDIR/out.log" "$BATS_TEST_TMPDIR/err.log" \
        reads_1.fastq.gz reads_2.fastq.gz
    [ "$status" -eq 0 ]
    [[ "$(cat "$RUN_TOOL_LOG")" == "fastqc -t 18 -o $BATS_TEST_TMPDIR/out reads_1.fastq.gz reads_2.fastq.gz" ]]
}

@test "qc::multiqc points at the input directory and writes to the output directory" {
    stub_run_tool
    run qc::multiqc "$BATS_TEST_TMPDIR/in" "$BATS_TEST_TMPDIR/out" \
        "$BATS_TEST_TMPDIR/out.log" "$BATS_TEST_TMPDIR/err.log"
    [ "$status" -eq 0 ]
    [[ "$(cat "$RUN_TOOL_LOG")" == "multiqc $BATS_TEST_TMPDIR/in -o $BATS_TEST_TMPDIR/out" ]]
}
