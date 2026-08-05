#!/usr/bin/env bats

load 'test_helper'

@test "mask::build_database names the BuildDatabase output after the given db name" {
    stub_run_tool
    run mask::build_database trichoderma_harzianum "$FIXTURES/sample_assembly.fasta" out.log err.log
    [ "$status" -eq 0 ]
    [ "$(cat "$RUN_TOOL_LOG")" = "BuildDatabase -name trichoderma_harzianum $FIXTURES/sample_assembly.fasta" ]
}

@test "mask::repeatmodeler runs RepeatModeler with its working directory as cwd" {
    RUN_TOOL_LOG="$BATS_TEST_TMPDIR/run_tool.log"
    : > "$RUN_TOOL_LOG"
    io::run_tool() {
        local env="$1"
        shift
        printf '%s\t%s\n' "$(pwd)" "$*" >> "$RUN_TOOL_LOG"
    }
    local work_dir="$BATS_TEST_TMPDIR/repeatmodeler"
    mkdir -p "$work_dir"

    run mask::repeatmodeler "$work_dir" trichoderma_harzianum 4 out.log err.log
    [ "$status" -eq 0 ]
    local expected
    expected="$(printf '%s\t%s' "$(cd "$work_dir" && pwd)" "RepeatModeler -database trichoderma_harzianum -pa 4 -LTRStruct")"
    [ "$(cat "$RUN_TOOL_LOG")" = "$expected" ]
}

@test "mask::repeatmasker soft-masks (-xsmall) with the custom library and a GFF track" {
    stub_run_tool
    run mask::repeatmasker repeat_lib.fa 4 /abs/out_dir "$FIXTURES/sample_assembly.fasta" out.log err.log
    [ "$status" -eq 0 ]
    [ "$(cat "$RUN_TOOL_LOG")" = "RepeatMasker -lib repeat_lib.fa -pa 4 -xsmall -gff -dir /abs/out_dir $FIXTURES/sample_assembly.fasta" ]
}
