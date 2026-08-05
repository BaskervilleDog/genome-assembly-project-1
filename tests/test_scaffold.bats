#!/usr/bin/env bats

load 'test_helper'

@test "scaffold::ragtag scaffolds against the reference and retains unplaced contigs" {
    stub_run_tool
    run scaffold::ragtag reference.fasta polished.fasta out_dir 18 out.log err.log
    [ "$status" -eq 0 ]
    [ "$(cat "$RUN_TOOL_LOG")" = "ragtag.py scaffold reference.fasta polished.fasta -o out_dir -t 18 -u" ]
}
