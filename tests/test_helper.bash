# Shared bats setup, loaded from every tests/*.bats file with:
#   load 'test_helper'
#
# Sources every lib/*.sh file so tests can call functions directly - this
# project has no package/install step, so there's no `source`-free way to
# pull them in otherwise.

ROOT="$(cd -- "$(dirname -- "${BATS_TEST_FILENAME}")/.." && pwd)"
FIXTURES="$ROOT/tests/fixtures"

for _f in "$ROOT"/lib/*.sh; do
    # shellcheck source=/dev/null
    source "$_f"
done
unset _f

# Replace io::run_tool - the one seam every external bioinformatics tool
# call in lib/*.sh goes through - with a stub that never touches mamba or
# a real conda environment. Mirrors the real function's `local env="$1";
# shift` so what gets recorded is the tool command exactly as it would run
# inside that env, not the env name. Every call is:
#   1. appended (space-joined args, after the env is shifted off) as one
#      line to $RUN_TOOL_LOG, so a test can assert exactly which
#      tool/flags a lib/*.sh function would have run, and in what order;
#   2. echoed to stdout as `stub-stdout: <args>`, so a test can also
#      assert that a function wired its `> "$out"` redirection at the
#      right call site.
#
# Call this from inside a @test before exercising a lib/*.sh function.
stub_run_tool() {
    RUN_TOOL_LOG="$BATS_TEST_TMPDIR/run_tool.log"
    : > "$RUN_TOOL_LOG"
    io::run_tool() {
        local env="$1"
        shift
        printf '%s\n' "$*" >> "$RUN_TOOL_LOG"
        printf 'stub-stdout: %s\n' "$*"
    }
}
