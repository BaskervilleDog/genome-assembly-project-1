#!/usr/bin/env bats

load 'test_helper'

@test "io::run_tool invokes mamba run -n <env> <command...>" {
    mamba() { printf '%s\n' "$*"; }
    run io::run_tool flye292 flye --version
    [ "$status" -eq 0 ]
    [ "$output" = "run -n flye292 flye --version" ]
}

@test "io::log_step prints a timestamped header containing the message" {
    run io::log_step "Assembly: Flye (repeat graph)"
    [ "$status" -eq 0 ]
    [[ "$output" == *"STEP: Assembly: Flye (repeat graph)"* ]]
}

@test "io::log_info prints a timestamped one-line message" {
    run io::log_info "wrote 00_downloads/SRR10848482_1.fastq.gz"
    [ "$status" -eq 0 ]
    [[ "$output" == *"[INFO] wrote 00_downloads/SRR10848482_1.fastq.gz"* ]]
}

@test "io::require_file fails on a missing file" {
    run io::require_file "$FIXTURES/does_not_exist.fastq"
    [ "$status" -ne 0 ]
    [[ "$output" == *"missing or empty"* ]]
}

@test "io::require_file fails on an empty file" {
    local empty="$BATS_TEST_TMPDIR/empty.fastq"
    : > "$empty"
    run io::require_file "$empty"
    [ "$status" -ne 0 ]
}

@test "io::require_file succeeds on a non-empty file" {
    run io::require_file "$FIXTURES/sample_reads.fastq"
    [ "$status" -eq 0 ]
}
