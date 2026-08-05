#!/usr/bin/env bats

load 'test_helper'

@test "download::ensure_seqfetcher clones via git when the directory is missing" {
    stub_run_tool
    local dir="$BATS_TEST_TMPDIR/seqfetcher"
    local git_log="$BATS_TEST_TMPDIR/git.log"
    git() {
        echo "$*" >> "$git_log"
        mkdir -p "$dir"
        printf '#!/usr/bin/env bash\n' > "$dir/install.sh"
    }
    run download::ensure_seqfetcher "$dir"
    [ "$status" -eq 0 ]
    [[ "$(cat "$git_log")" == *"clone https://github.com/BaskervilleDog/seqfetcher.git"* ]]
    [ "$(cat "$RUN_TOOL_LOG")" = "./install.sh" ]
}

@test "download::ensure_seqfetcher skips cloning when the directory already exists" {
    stub_run_tool
    local dir="$BATS_TEST_TMPDIR/seqfetcher"
    mkdir -p "$dir"
    printf '#!/usr/bin/env bash\n' > "$dir/install.sh"
    local git_log="$BATS_TEST_TMPDIR/git.log"
    git() { echo "$*" >> "$git_log"; }

    run download::ensure_seqfetcher "$dir"
    [ "$status" -eq 0 ]
    [ ! -f "$git_log" ]
    [ "$(cat "$RUN_TOOL_LOG")" = "./install.sh" ]
}

@test "download::sra_reads writes the accession list and moves downloaded FASTQs" {
    stub_run_tool
    local seqfetcher_dir="$BATS_TEST_TMPDIR/seqfetcher"
    mkdir -p "$seqfetcher_dir/downloads/fastq"
    touch "$seqfetcher_dir/downloads/fastq/SRR10848482_1.fastq.gz"
    local dest="$BATS_TEST_TMPDIR/00_downloads"

    run download::sra_reads "$seqfetcher_dir" "$dest" SRR10848482 SRR10848483
    [ "$status" -eq 0 ]
    [ "$(cat "$seqfetcher_dir/accession_list.txt")" = $'SRR10848482\nSRR10848483' ]
    [ "$(cat "$RUN_TOOL_LOG")" = "seqfetcher download --sra-method parallel --sra-accession-file accession_list.txt" ]
    [ -f "$dest/SRR10848482_1.fastq.gz" ]
}
