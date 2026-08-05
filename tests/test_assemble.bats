#!/usr/bin/env bats

load 'test_helper'

@test "assemble::flye passes --pacbio-corr, genome size, out dir and threads" {
    stub_run_tool
    run assemble::flye clean.fastq.gz 40m out_dir 18 out.log err.log
    [ "$status" -eq 0 ]
    [[ "$(cat "$RUN_TOOL_LOG")" == "flye --pacbio-corr clean.fastq.gz --genome-size 40m --out-dir out_dir --threads 18" ]]
}

@test "assemble::raven writes stdout to the given FASTA path" {
    stub_run_tool
    local out="$BATS_TEST_TMPDIR/raven_assembly.fasta"
    run assemble::raven clean.fastq.gz 18 "$out" "$BATS_TEST_TMPDIR/err.log"
    [ "$status" -eq 0 ]
    [[ "$(cat "$RUN_TOOL_LOG")" == "raven --threads 18 clean.fastq.gz" ]]
    [[ "$(cat "$out")" == *"stub-stdout: raven --threads 18 clean.fastq.gz"* ]]
}

@test "assemble::wtdbg2 runs graph construction then consensus calling, in order" {
    stub_run_tool
    run assemble::wtdbg2 clean.fastq.gz 40m 18 out_prefix out.fasta \
        "$BATS_TEST_TMPDIR/out.log" "$BATS_TEST_TMPDIR/err.log"
    [ "$status" -eq 0 ]
    [ "$(wc -l < "$RUN_TOOL_LOG")" -eq 2 ]
    [ "$(sed -n 1p "$RUN_TOOL_LOG")" = "wtdbg2 -x rs -g 40m -t 18 -i clean.fastq.gz -fo out_prefix" ]
    [ "$(sed -n 2p "$RUN_TOOL_LOG")" = "wtpoa-cns -t 18 -i out_prefix.ctg.lay.gz -fo out.fasta" ]
}
