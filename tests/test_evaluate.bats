#!/usr/bin/env bats

load 'test_helper'

@test "evaluate::busco requests genome mode, the given lineage, and forces overwrite" {
    stub_run_tool
    run evaluate::busco "$FIXTURES/sample_assembly.fasta" hypocreaceae_odb12 18 out_dir out.log err.log
    [ "$status" -eq 0 ]
    [ "$(cat "$RUN_TOOL_LOG")" = "busco -i $FIXTURES/sample_assembly.fasta -f -m genome -l hypocreaceae_odb12 -c 18 -o out_dir" ]
}

@test "evaluate::quast compares every given assembly with matching --labels" {
    stub_run_tool
    run evaluate::quast 18 out_dir out.log err.log "Flye,Raven,wtdbg2" \
        flye.fasta raven.fasta wtdbg2.fasta
    [ "$status" -eq 0 ]
    [ "$(cat "$RUN_TOOL_LOG")" = "quast flye.fasta raven.fasta wtdbg2.fasta --labels Flye,Raven,wtdbg2 -t 18 -o out_dir" ]
}

@test "evaluate::quast_reference passes -r and --fragmented" {
    stub_run_tool
    run evaluate::quast_reference polished.fasta reference.fasta 18 out_dir out.log err.log
    [ "$status" -eq 0 ]
    [ "$(cat "$RUN_TOOL_LOG")" = "quast polished.fasta -r reference.fasta -t 18 --fragmented -o out_dir" ]
}

@test "evaluate::meryl_count builds a k-mer database from every given read file" {
    stub_run_tool
    run evaluate::meryl_count 21 reads.meryl out.log err.log r1.fastq.gz r2.fastq.gz r3.fastq.gz r4.fastq.gz
    [ "$status" -eq 0 ]
    [ "$(cat "$RUN_TOOL_LOG")" = "meryl count k=21 output reads.meryl r1.fastq.gz r2.fastq.gz r3.fastq.gz r4.fastq.gz" ]
}

@test "evaluate::merqury runs merqury.sh with the run directory as cwd" {
    # A custom stub (not stub_run_tool) so the assertion can capture *where*
    # io::run_tool actually ran from, not just what it was called with -
    # this is the one lib/*.sh function that requires a specific cwd.
    RUN_TOOL_LOG="$BATS_TEST_TMPDIR/run_tool.log"
    : > "$RUN_TOOL_LOG"
    io::run_tool() {
        local env="$1"
        shift
        printf '%s\t%s\n' "$(pwd)" "$*" >> "$RUN_TOOL_LOG"
    }
    local run_dir="$BATS_TEST_TMPDIR/merqury"
    mkdir -p "$run_dir"

    run evaluate::merqury "$run_dir" reads.meryl polished.fasta merqury_output out.log err.log
    [ "$status" -eq 0 ]
    local expected
    expected="$(printf '%s\t%s' "$(cd "$run_dir" && pwd)" "merqury.sh reads.meryl polished.fasta merqury_output")"
    [ "$(cat "$RUN_TOOL_LOG")" = "$expected" ]
}
