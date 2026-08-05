#!/usr/bin/env bats

load 'test_helper'

@test "polish::racon_round aligns with minimap2 then corrects with racon, in order" {
    stub_run_tool
    run polish::racon_round assembly.fasta reads.fastq.gz map-pb 18 \
        out.paf out.fasta minimap.log racon.log
    [ "$status" -eq 0 ]
    [ "$(wc -l < "$RUN_TOOL_LOG")" -eq 2 ]
    [ "$(sed -n 1p "$RUN_TOOL_LOG")" = "minimap2 -t 18 -x map-pb assembly.fasta reads.fastq.gz" ]
    [ "$(sed -n 2p "$RUN_TOOL_LOG")" = "racon -t 18 reads.fastq.gz out.paf assembly.fasta" ]
}

@test "polish::polypolish indexes once, aligns+filters each library, then polishes jointly in library order" {
    # io::run_tool is fully stubbed (see test_helper.bash): it never runs a
    # real `polypolish`/`bwa`, so none of the intermediate SAM files this
    # function references actually get created - the final `rm -f` cleanup
    # is a no-op against them (rm -f does not error on missing files),
    # which is why this test can assert on command construction alone.
    stub_run_tool
    local work="$BATS_TEST_TMPDIR/work" log_dir="$BATS_TEST_TMPDIR/logs"
    mkdir -p "$work" "$log_dir"

    run polish::polypolish assembly.fasta 18 "$work" "$BATS_TEST_TMPDIR/polished.fasta" "$log_dir" \
        libA_1.fastq.gz libA_2.fastq.gz libA \
        libB_1.fastq.gz libB_2.fastq.gz libB
    [ "$status" -eq 0 ]

    run grep -c '^bwa index assembly.fasta$' "$RUN_TOOL_LOG"
    [ "$output" -eq 1 ]

    # The final `polypolish polish` call must list libA's SAMs before
    # libB's, matching the order libraries were passed in.
    run grep '^polypolish polish assembly.fasta' "$RUN_TOOL_LOG"
    [[ "$output" == *"filtered_libA_1.sam filtered_libA_2.sam filtered_libB_1.sam filtered_libB_2.sam"* ]]
}

@test "polish::nextpolish_round writes a config file pointing at the given genome and workdir" {
    stub_run_tool
    local cfg="$BATS_TEST_TMPDIR/run.cfg"
    run polish::nextpolish_round genome.fasta sgs.fofn "$BATS_TEST_TMPDIR/work" "$cfg" \
        "$BATS_TEST_TMPDIR/out.log" "$BATS_TEST_TMPDIR/err.log"
    [ "$status" -eq 0 ]
    [ -f "$cfg" ]
    [[ "$(cat "$cfg")" == *"genome = genome.fasta"* ]]
    [[ "$(cat "$cfg")" == *"workdir = $BATS_TEST_TMPDIR/work"* ]]
    [[ "$(cat "$cfg")" == *"sgs_fofn = sgs.fofn"* ]]
    [[ "$(cat "$RUN_TOOL_LOG")" == "nextPolish $cfg" ]]
}
