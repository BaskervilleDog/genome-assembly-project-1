#!/usr/bin/env bats
# Integration-ish checks that don't require any of the 20 mamba
# environments or real sequencing data. A real end-to-end run needs both
# (and takes hours-to-days) so isn't something a test suite can exercise;
# what *is* cheap to check on every commit is that nothing is left
# syntactically broken and that the entry-point scripts still wire up to
# every lib/*.sh function they reference.

load 'test_helper'

@test "every lib/*.sh and scripts/*.sh file is syntactically valid bash" {
    for f in "$ROOT"/lib/*.sh "$ROOT"/scripts/*.sh; do
        run bash -n "$f"
        [ "$status" -eq 0 ] || {
            echo "syntax error in $f:"
            echo "$output"
            return 1
        }
    done
}

@test "scripts/check_deps.sh runs to completion (regardless of what's installed here)" {
    run bash "$ROOT/scripts/check_deps.sh"
    # 0 = everything present, 1 = something's missing (expected on a
    # machine without the 20 mamba envs) - anything else means the script
    # itself is broken (syntax error, unbound variable, etc.).
    [ "$status" -eq 0 ] || [ "$status" -eq 1 ]
    [[ "$output" == *"required mamba environments"* ]]
}

@test "01_download.sh references only functions that exist in lib/*.sh" {
    run grep -ohE '\b(io|download|qc)::[a-z0-9_]+' "$ROOT/scripts/01_download.sh"
    [ "$status" -eq 0 ]
    while IFS= read -r fn; do
        [ -z "$fn" ] && continue
        run type -t "$fn"
        [ "$output" = "function" ]
    done <<< "$output"
}

@test "02_assembly.sh references only functions that exist in lib/*.sh" {
    run grep -ohE '\b(io|qc|filter|decontaminate|assemble|polish)::[a-z0-9_]+' "$ROOT/scripts/02_assembly.sh"
    [ "$status" -eq 0 ]
    while IFS= read -r fn; do
        [ -z "$fn" ] && continue
        run type -t "$fn"
        [ "$output" = "function" ]
    done <<< "$output"
}

@test "03_quality_assessment.sh references only functions that exist in lib/*.sh" {
    run grep -ohE '\b(io|evaluate|scaffold|mask)::[a-z0-9_]+' "$ROOT/scripts/03_quality_assessment.sh"
    [ "$status" -eq 0 ]
    while IFS= read -r fn; do
        [ -z "$fn" ] && continue
        run type -t "$fn"
        [ "$output" = "function" ]
    done <<< "$output"
}
