#!/usr/bin/env bash
# Generate function-call and module-dependency graphs from the source code.
#
# This is the "keep the diagram honest" tool: instead of hand-drawing an
# architecture diagram that drifts out of date (the ASCII pipeline diagram
# in README.md is exactly that kind of diagram), it statically scans
# lib/*.sh and the entry-point scripts with grep/awk and derives two
# Graphviz diagrams straight from the code:
#
#     graphs/call_graph.svg          which function calls which function
#     graphs/dependency_graph.svg    which file sources/invokes what
#
# Run after changing the library's structure:
#
#     bash scripts/generate_function_graph.sh
#
# Requires the Graphviz `dot` CLI on PATH to render SVGs (`winget install
# graphviz` / `brew install graphviz` / `apt install graphviz` / `conda
# install graphviz`). If `dot` isn't found, the .dot source files are
# still written - render them later, or paste them into
# https://dreampuf.github.io/GraphvizOnline/.
#
# There's no `ast`/`parse()` to lean on here - bash has no standard parser
# API - so this is regex/awk best-effort, not a real static-analysis tool:
# it resolves direct calls to functions named `module::function` (the
# convention every lib/*.sh file in this project follows) and skips
# anything it can't confidently resolve.
#
# Note on what this graph does and doesn't show: every external
# bioinformatics tool call in lib/*.sh goes through io::run_tool (`mamba
# run -n <env> <tool> ...`), so the tool name itself is just a string
# argument, not a bare command word - this script will correctly show
# lib.download calling lib.io, but it won't draw a "flye"/"busco"/"quast"
# node for every tool the pipeline invokes. envs/*.yml and specs/*.txt are
# the source of truth for the tool inventory; this graph is about the
# *pipeline's own* structure (which stage calls which, in what order).

set -uo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
LIB_DIR="$ROOT/lib"
OUTPUT_DIR="$ROOT/graphs"
SELF="$(basename -- "${BASH_SOURCE[0]}")"

mkdir -p -- "$OUTPUT_DIR"

# --- 1. Discover source files -----------------------------------------
# lib/*.sh are the "package"; scripts/*.sh (other than this file) are
# entry points, shown as callers/sourcers rather than reusable modules.
mapfile -t FILES < <(
    { find "$LIB_DIR" -maxdepth 1 -name '*.sh' | sort
      find "$ROOT/scripts" -maxdepth 1 -name '*.sh' ! -name "$SELF" | sort
    }
)

dotted_module() {
    # lib/io.sh -> lib.io ; scripts/02_assembly.sh -> scripts.02_assembly
    local rel="${1#"$ROOT"/}"
    rel="${rel%.sh}"
    echo "${rel//\//.}"
}

is_script_module() {
    [[ "$1" == scripts.* ]]
}

# --- 2. Strip heredoc bodies --------------------------------------------
# lib/polish.sh's NextPolish config and 02_assembly.sh's sgs.fofn heredoc
# are data, not bash code, and full of tokens (`genome = ...`,
# `[sgs_option]`) that would otherwise be misread as function/command
# calls. Blank them out, keeping line numbers intact so later line-range
# logic still lines up.
strip_heredocs() {
    awk '
        BEGIN { in_here = 0; delim = "" }
        {
            if (in_here) {
                line = $0
                gsub(/^[ \t]+/, "", line)
                if (line == delim) { in_here = 0; print ""; next }
                print ""
                next
            }
            if (match($0, /<<[-~]?[ \t]*["'"'"']?[A-Za-z_][A-Za-z0-9_]*["'"'"']?/)) {
                tok = substr($0, RSTART, RLENGTH)
                gsub(/<<[-~]?[ \t]*["'"'"']?/, "", tok)
                gsub(/["'"'"']$/, "", tok)
                delim = tok
                in_here = 1
            }
            print
        }
    ' "$1"
}

# --- 2b. Mask single-quoted string contents -----------------------------
# Blank out everything between single quotes (state carries across lines)
# while leaving the quotes themselves and the rest of the line intact.
# Double-quoted content is buffered and only blanked if it does *not*
# contain a `$(` command substitution - an echo/printf message that
# merely *mentions* a function or tool name must not read as a call, but
# `SCAFFOLDED="$(...)"`-style assignments have to stay visible. This is a
# parity approximation, not real parsing.
mask_quotes() {
    # state: 0 = unquoted, 1 = inside '...', 2 = inside "..."
    awk '
        BEGIN { state = 0; buf = ""; has_cmdsub = 0 }
        {
            out = ""
            n = length($0)
            for (i = 1; i <= n; i++) {
                c = substr($0, i, 1)
                if (state == 0) {
                    if (c == "#") { break }  # unquoted "#" -> rest of line is a comment
                    if (c == "\047") { state = 1; out = out c; continue }
                    if (c == "\042") { state = 2; buf = ""; has_cmdsub = 0; out = out c; continue }
                    out = out c
                    continue
                }
                if (state == 1) {
                    if (c == "\047") { state = 0; out = out c } else { out = out " " }
                    continue
                }
                # state == 2: inside "..."
                if (c == "\042") {
                    state = 0
                    if (has_cmdsub) {
                        out = out buf
                    } else {
                        for (k = 1; k <= length(buf); k++) out = out " "
                    }
                    out = out c
                    continue
                }
                buf = buf c
                if (length(buf) >= 2 && substr(buf, length(buf) - 1, 2) == "$(") has_cmdsub = 1
            }
            print out
        }
    ' "$1"
}

# Args: $1 = source file, $2 = destination path for the fully preprocessed copy
preprocess() {
    local src="$1" dest="$2" tmp
    tmp="$(mktemp)"
    strip_heredocs "$src" > "$tmp"
    mask_quotes "$tmp" > "$dest"
    rm -f "$tmp"
}

# --- 3. Bash keywords / builtins / coreutils treated as noise ----------
NOISE=(if then elif else fi for while until do done case esac in function
    select time coproc
    alias bg bind break builtin caller cd command compgen complete continue
    declare dirs disown echo enable eval exec exit export false fc fg
    getopts hash help history jobs kill let local logout mapfile popd printf
    pushd pwd read readarray readonly return set shift shopt source suspend
    test times trap true type typeset ulimit umask unalias unset wait
    awk sed grep egrep fgrep sort uniq cut tr wc head tail cat tac paste
    join comm cp mv rm mkdir rmdir touch chmod chown ls find xargs basename
    dirname date sleep tee printenv env seq expr dd du df realpath readlink
    diff stat main)
declare -A IS_NOISE=()
for w in "${NOISE[@]}"; do IS_NOISE["$w"]=1; done

# --- 4. Pass 1: discover every function definition ----------------------
# Records: file \t start_line \t end_line \t raw_name
FUNC_RECORDS="$(mktemp)"
trap 'rm -f "$FUNC_RECORDS" "$STRIPPED_DIR"/* 2>/dev/null; rmdir "$STRIPPED_DIR" 2>/dev/null' EXIT
STRIPPED_DIR="$(mktemp -d)"

for file in "${FILES[@]}"; do
    stripped="$STRIPPED_DIR/$(basename -- "$file").stripped"
    preprocess "$file" "$stripped"

    awk -v file="$file" '
        function flush_pending() {
            if (name != "") { print file "\t" start "\t" NR-1 "\t" name }
            name = ""; depth = 0
        }
        {
            if (name == "" && match($0, /^[ \t]*(function[ \t]+)?[A-Za-z_][A-Za-z0-9_:]*\(\)[ \t]*\{/)) {
                line = $0
                sub(/^[ \t]*function[ \t]+/, "", line)
                sub(/\(\).*/, "", line)
                gsub(/^[ \t]+|[ \t]+$/, "", line)
                name = line
                start = NR
                depth = gsub(/\{/, "{", $0) - gsub(/\}/, "}", $0)
                if (depth <= 0) { print file "\t" start "\t" NR "\t" name; name = "" }
                next
            }
            if (name != "") {
                depth += gsub(/\{/, "{", $0) - gsub(/\}/, "}", $0)
                if (depth <= 0) { print file "\t" start "\t" NR "\t" name; name = "" }
            }
        }
        END { flush_pending() }
    ' "$stripped" >> "$FUNC_RECORDS"
done

# --- 5. Build the node table: raw_name -> node_id, plus per-file line ranges
declare -A NODE_ID=()          # raw_name -> "dotted.module.short_name"
declare -A FUNC_FILE=()        # raw_name -> file
declare -A FUNC_START=()       # raw_name -> start line
declare -A FUNC_END=()         # raw_name -> end line
declare -a ALL_RAW_NAMES=()

while IFS=$'\t' read -r file start end raw_name; do
    [[ -z "$raw_name" ]] && continue
    mod="$(dotted_module "$file")"
    prefix="$(basename -- "$file" .sh)::"
    if [[ "$raw_name" == "$prefix"* ]]; then
        short="${raw_name#"$prefix"}"
    else
        short="$raw_name"
    fi
    NODE_ID["$raw_name"]="$mod.$short"
    FUNC_FILE["$raw_name"]="$file"
    FUNC_START["$raw_name"]="$start"
    FUNC_END["$raw_name"]="$end"
    ALL_RAW_NAMES+=("$raw_name")
done < "$FUNC_RECORDS"

# --- 6. Pass 2: call edges (best-effort: whole-word match inside a body) --
CALL_EDGES="$(mktemp)"
for caller in "${ALL_RAW_NAMES[@]}"; do
    file="${FUNC_FILE[$caller]}"
    stripped="$STRIPPED_DIR/$(basename -- "$file").stripped"
    start="${FUNC_START[$caller]}"
    end="${FUNC_END[$caller]}"
    body="$(sed -n "$((start + 1)),${end}p" "$stripped")"
    for callee in "${ALL_RAW_NAMES[@]}"; do
        [[ "$callee" == "$caller" ]] && continue
        if grep -qF -- "$callee" <<< "$body"; then
            echo "${NODE_ID[$caller]}	${NODE_ID[$callee]}" >> "$CALL_EDGES"
        fi
    done
done

# Top-level calls in entry-point scripts (e.g. the STEP-by-STEP body of
# 01_download.sh/02_assembly.sh/03_quality_assessment.sh, none of which
# wrap their logic in a `main`): lines outside every discovered function's
# range, attributed to the module itself rather than to one function.
for file in "${FILES[@]}"; do
    mod="$(dotted_module "$file")"
    is_script_module "$mod" || continue
    stripped="$STRIPPED_DIR/$(basename -- "$file").stripped"
    total_lines="$(wc -l < "$stripped")"

    covered="$(mktemp)"
    : > "$covered"
    for name in "${ALL_RAW_NAMES[@]}"; do
        [[ "${FUNC_FILE[$name]}" == "$file" ]] || continue
        seq "${FUNC_START[$name]}" "${FUNC_END[$name]}" >> "$covered"
    done

    top_level="$(mktemp)"
    for ((ln = 1; ln <= total_lines; ln++)); do
        grep -qx "$ln" "$covered" || sed -n "${ln}p" "$stripped" >> "$top_level"
    done

    for callee in "${ALL_RAW_NAMES[@]}"; do
        if grep -qF -- "$callee" "$top_level"; then
            echo "${mod}	${NODE_ID[$callee]}" >> "$CALL_EDGES"
        fi
    done
    rm -f "$covered" "$top_level"
done

sort -u "$CALL_EDGES" -o "$CALL_EDGES"

# --- 7. Dependency graph: internal sourcing + external commands ---------
DEP_EDGES="$(mktemp)"        # solid: sources / invokes
DEP_CHECK_EDGES="$(mktemp)"  # dashed: `command -v x` presence checks

declare -A LIB_MODULES=()
for file in "${FILES[@]}"; do
    mod="$(dotted_module "$file")"
    [[ "$mod" == lib.* ]] && LIB_MODULES["$mod"]=1
done

for file in "${FILES[@]}"; do
    mod="$(dotted_module "$file")"
    stripped="$STRIPPED_DIR/$(basename -- "$file").stripped"

    # A `for f in .../lib/*.sh; do source "$f"; done`-style loop: treat as
    # sourcing every lib module (the pattern every entry-point script uses).
    if grep -qE '/lib/\*\.sh' "$stripped"; then
        for lib_mod in "${!LIB_MODULES[@]}"; do
            echo "${mod}	${lib_mod}" >> "$DEP_EDGES"
        done
    fi
    # A direct `source path/to/x.sh` / `. path/to/x.sh` naming one file.
    while IFS= read -r target; do
        [[ -z "$target" ]] && continue
        target_mod="lib.$(basename -- "$target" .sh)"
        [[ -n "${LIB_MODULES[$target_mod]:-}" ]] && echo "${mod}	${target_mod}" >> "$DEP_EDGES"
    done < <(grep -oE '(source|\.)[ \t]+"?[^"; ]*/([A-Za-z0-9_-]+\.sh)' "$stripped" \
        | grep -oE '[A-Za-z0-9_-]+\.sh$')

    # `command -v NAME` presence checks (what scripts/check_deps.sh does).
    while IFS= read -r tool; do
        echo "${mod}	${tool}" >> "$DEP_CHECK_EDGES"
    done < <(grep -oE 'command -v [A-Za-z0-9_.-]+' "$stripped" | awk '{print $3}')

    # External commands actually invoked: first word of each `;`/`|`/&&/||-
    # separated segment, minus keywords/builtins/coreutils/known functions.
    while IFS= read -r word; do
        [[ -z "$word" ]] && continue
        [[ -n "${IS_NOISE[$word]:-}" ]] && continue
        case "$word" in *::*) continue ;; esac  # a lib/script function call, not external
        [[ -n "${NODE_ID[$word]:-}" ]] && continue
        echo "${mod}	${word}" >> "$DEP_EDGES"
    done < <(sed -E 's/(&&|\|\||[;|])/\n/g' "$stripped" \
        | sed -E 's/^[ \t]+//' \
        | grep -vE '^#' \
        | grep -vE '^[A-Za-z_][A-Za-z0-9_]*\+?=' \
        | grep -oE '^[A-Za-z_][A-Za-z0-9_.:-]*')
done

sort -u "$DEP_EDGES" -o "$DEP_EDGES"
sort -u "$DEP_CHECK_EDGES" -o "$DEP_CHECK_EDGES"

# --- 8. Render -----------------------------------------------------------
render_dot() {
    local dot_source="$1" out_stem="$2"
    printf '%s\n' "$dot_source" > "${out_stem}.dot"
    if command -v dot >/dev/null 2>&1; then
        if dot -Tsvg "${out_stem}.dot" -o "${out_stem}.svg" 2>/tmp/dot_err.$$; then
            echo "wrote ${out_stem#"$ROOT"/}.dot and ${out_stem#"$ROOT"/}.svg"
        else
            echo "wrote ${out_stem#"$ROOT"/}.dot but 'dot' failed:"
            cat /tmp/dot_err.$$ >&2
        fi
        rm -f "/tmp/dot_err.$$"
    else
        echo "wrote ${out_stem#"$ROOT"/}.dot (Graphviz 'dot' not found on PATH - install Graphviz to auto-render SVGs)"
    fi
}

palette=(\#cfe8ff \#d9f2d9 \#ffe8cc \#f4d9ff \#ffe0e0 \#e0f7f5 \#fff6cc)

build_call_graph_dot() {
    echo 'digraph call_graph {'
    echo '  rankdir="LR";'
    echo '  node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=11];'
    echo '  edge [color="#888888"];'

    local i=0
    for file in "${FILES[@]}"; do
        mod="$(dotted_module "$file")"
        color="${palette[$((i % ${#palette[@]}))]}"
        i=$((i + 1))
        echo "  subgraph \"cluster_${mod}\" {"
        echo "    label=\"${mod}\"; style=\"dashed\"; fontname=\"Helvetica\"; fontsize=11;"
        if is_script_module "$mod"; then
            echo "    \"${mod}\" [label=\"${mod}\", fillcolor=\"${color}\", shape=component];"
        fi
        for name in "${ALL_RAW_NAMES[@]}"; do
            [[ "${FUNC_FILE[$name]}" == "$file" ]] || continue
            id="${NODE_ID[$name]}"
            short="${id#"$mod".}"
            echo "    \"${id}\" [label=\"${short}\", fillcolor=\"${color}\", shape=box];"
        done
        echo "  }"
    done

    while IFS=$'\t' read -r caller callee; do
        [[ -z "$caller" ]] && continue
        echo "  \"${caller}\" -> \"${callee}\";"
    done < "$CALL_EDGES"
    echo '}'
}

build_dependency_graph_dot() {
    echo 'digraph dependency_graph {'
    echo '  rankdir="LR";'
    echo '  node [fontname="Helvetica", fontsize=11];'
    echo '  edge [color="#888888"];'

    for file in "${FILES[@]}"; do
        mod="$(dotted_module "$file")"
        echo "  \"${mod}\" [label=\"${mod}\", shape=box, style=\"rounded,filled\", fillcolor=\"#cfe8ff\"];"
    done

    local -A external=()
    while IFS=$'\t' read -r src dst; do
        [[ -z "$src" ]] && continue
        if [[ "$dst" == lib.* || "$dst" == scripts.* ]]; then
            :  # internal node already declared above
        else
            external["$dst"]=1
        fi
    done < <(cat "$DEP_EDGES" "$DEP_CHECK_EDGES")
    for tool in "${!external[@]}"; do
        echo "  \"${tool}\" [label=\"${tool}\", shape=diamond, style=filled, fillcolor=\"#ffe8cc\"];"
    done

    while IFS=$'\t' read -r src dst; do
        [[ -z "$src" ]] && continue
        echo "  \"${src}\" -> \"${dst}\";"
    done < "$DEP_EDGES"
    while IFS=$'\t' read -r src dst; do
        [[ -z "$src" ]] && continue
        echo "  \"${src}\" -> \"${dst}\" [style=dashed, label=\"checks for\", fontsize=9];"
    done < "$DEP_CHECK_EDGES"

    echo '}'
}

render_dot "$(build_call_graph_dot)" "$OUTPUT_DIR/call_graph"
render_dot "$(build_dependency_graph_dot)" "$OUTPUT_DIR/dependency_graph"

rm -f "$CALL_EDGES" "$DEP_EDGES" "$DEP_CHECK_EDGES"
