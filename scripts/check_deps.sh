#!/usr/bin/env bash
# There's no lockfile for a pipeline built on 20 separate mamba
# environments and a handful of system tools - this script is what
# stands in for one. It checks that every tool this pipeline actually
# touches is on PATH, and that every mamba environment 01_download.sh /
# 02_assembly.sh / 03_quality_assessment.sh expect (per envs/*.yml) is
# installed - so a run fails in seconds, not 12 hours into
# RepeatModeler2, because one env was never created on this machine.
#
# Run with:
#   bash scripts/check_deps.sh

set -uo pipefail

ok=0
missing_required=0

check() {
    local name="$1" required="$2" note="$3"
    if command -v "$name" >/dev/null 2>&1; then
        printf '  [x] %-14s %s\n' "$name" "$(command -v "$name")"
        ok=$((ok + 1))
    elif [[ "$required" == "required" ]]; then
        printf '  [ ] %-14s MISSING (required) - %s\n' "$name" "$note"
        missing_required=$((missing_required + 1))
    else
        printf '  [ ] %-14s missing (optional) - %s\n' "$name" "$note"
    fi
}

echo "bash: $BASH_VERSION"
if (( BASH_VERSINFO[0] < 4 )); then
    echo "  warning: bash >= 4 recommended (associative arrays used in scripts/generate_function_graph.sh)"
fi

echo "required (the pipeline scripts won't run without these):"
check mamba required "environment manager - every tool call in lib/*.sh goes through it (io::run_tool)"
check git   required "clones SeqFetcher (lib/download.sh)"
check pigz  required "parallel gzip - compresses decontaminated reads in scripts/02_assembly.sh"
check awk   required "used throughout scripts/*.sh"
check sed   required "used by scripts/generate_function_graph.sh"
check grep  required "used by scripts/generate_function_graph.sh"

echo
echo "required mamba environments (envs/*.yml / specs/*.txt are the source"
echo "of truth for exact tool versions - this only checks the env exists):"

if command -v mamba >/dev/null 2>&1; then
    mamba_envs="$(mamba env list 2>/dev/null | awk '{print $1}')"
    for env in \
        ncbi_tools seqkit longreadsum fastqc multiqc \
        fastplong fastp kraken2 \
        flye292 raven wtdbg2 \
        minimap2 racon1420 bwa polypolish nextpolish39 \
        busco quast merqury ragtag repeatmodeler
    do
        if grep -qx -- "$env" <<< "$mamba_envs"; then
            printf '  [x] %s\n' "$env"
            ok=$((ok + 1))
        else
            printf '  [ ] %-14s MISSING - create from envs/%s.yml (mamba env create -f envs/%s.yml)\n' \
                "$env" "$env" "$env"
            missing_required=$((missing_required + 1))
        fi
    done
else
    echo "  (skipped - mamba itself is missing, see above)"
fi

echo
echo "optional (only needed for development, not for running the pipeline):"
check bats       optional "test runner for tests/*.bats - install: npm i -g bats / brew install bats-core / apt install bats"
check dot        optional "renders graphs/*.svg in scripts/generate_function_graph.sh - without it, the .dot source is still written"
check shellcheck optional "lints lib/*.sh and scripts/*.sh - install: winget install shellcheck / brew install shellcheck / apt install shellcheck"

echo
if (( missing_required > 0 )); then
    echo "$missing_required required item(s) missing - install/create them before running the pipeline."
    exit 1
fi
echo "all required tools and environments present ($ok found)."
