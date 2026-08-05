#!/usr/bin/env bash
# Stage 1: pull raw reads from NCBI SRA via SeqFetcher (a wrapper around
# fasterq-dump and related utilities). This is the only lib/ file that
# touches the network or clones a git repo.

# Ensure the SeqFetcher clone exists at $dir and is installed into its
# mamba env. Safe to call every run - it no-ops (except reinstalling) once
# the clone exists.
#
# Args:
#   $1 - directory to clone/install SeqFetcher into
download::ensure_seqfetcher() {
    local dir="$1"
    if [[ ! -d "$dir" ]]; then
        io::log_info "Cloning SeqFetcher..."
        git clone https://github.com/BaskervilleDog/seqfetcher.git "$dir"
    fi
    ( cd "$dir" && chmod +x install.sh && io::run_tool ncbi_tools ./install.sh )
}

# Download a list of SRA accessions into $dest_dir via SeqFetcher.
#
# Args:
#   $1   - SeqFetcher directory (must already exist - see ensure_seqfetcher)
#   $2   - destination directory for the downloaded FASTQs
#   $3.. - one or more SRA accession IDs
download::sra_reads() {
    local seqfetcher_dir="$1" dest_dir="$2"
    shift 2
    local accessions=("$@")

    (
        cd "$seqfetcher_dir"
        printf '%s\n' "${accessions[@]}" > accession_list.txt
        io::run_tool ncbi_tools seqfetcher download \
            --sra-method parallel \
            --sra-accession-file accession_list.txt
    )

    mkdir -p -- "$dest_dir"
    mv "$seqfetcher_dir/downloads/fastq/"*.gz "$dest_dir/"
}
