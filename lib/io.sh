#!/usr/bin/env bash
# The one place `mamba run -n <env> ...` appears. Every external tool call
# in lib/*.sh goes through io::run_tool so there is a single seam to mock
# in tests/*.bats (redefine io::run_tool itself, capture the args it was
# called with) instead of needing 20 real mamba environments and a real
# PacBio/Illumina dataset installed just to test that a function builds
# the right command line.
#
# Also holds the timestamped logging helpers every lib/*.sh and
# scripts/*.sh entry point uses, and a small pre-flight guard for catching
# a missing input before an expensive multi-hour step starts.

# Run TOOL inside its mamba environment.
#
# Args:
#   $1   - mamba environment name (see envs/*.yml / specs/*.txt)
#   $2.. - the command (tool + its own flags) to run inside that environment
io::run_tool() {
    local env="$1"
    shift
    mamba run -n "$env" "$@"
}

# Timestamped section header, printed to stdout.
io::log_step() {
    echo ""
    echo "========================================================================"
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] STEP: $*"
    echo "========================================================================"
}

# Timestamped one-line info message, printed to stdout.
io::log_info() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] [INFO] $*"
}

# Fail fast (rather than hours into a run) if an expected input is missing
# or empty.
#
# Args:
#   $1 - path that must exist and be non-empty
io::require_file() {
    local path="$1"
    if [[ ! -s "$path" ]]; then
        echo "io::require_file: missing or empty: $path" >&2
        return 1
    fi
}
