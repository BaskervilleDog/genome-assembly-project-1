#!/usr/bin/env bash
#==============================================================================
# GENOME ANNOTATION PIPELINE (BRAKER3 + functional tools)
#==============================================================================
# Author: Gianlucca de Urzêda Alves
# Date: 03/02/2026
#
# Purpose:
#   Establish a pipeline for ab-initio and evidence-based genome annotation
#   using BRAKER3 as the authoritative gene set, followed by functional
#   annotation tools (DIAMOND+Swiss-Prot, InterProScan, dbCAN, antiSMASH),
#   producing standardized outputs ready for integration.
#
# Notes:
#   - This script does NOT install tools or databases.
#   - It assumes conda/mamba environments already exist.
#   - BRAKER3 outputs define the authoritative gene set.
#   - All downstream analyses use the SAME predicted protein FASTA.
#
#==============================================================================

#set -euo pipefail

#==============================================================================
# CONFIGURATION (edit these)
#==============================================================================

# --- Inputs ---
GENOME_FASTA="./flye_assembly/ragtag_output/ragtag.scaffold.clean.fasta"
SPECIES_TAG="Tharzianum_T11W"

# --- Threads ---
THREADS="12"

# --- Output root ---
OUTDIR="annotation_pipeline"

# --- Conda environments (must exist) ---
ENV_BRAKER="braker"
ENV_AGAT="agat"
ENV_BUSCO="busco"
ENV_SEQKIT="seqkit"
ENV_DIAMOND="diamond"
ENV_INTERPROSCAN="interproscan"
ENV_DBCAN="dbcan"
ENV_ANTISMASH="antismash"

# --- External databases / resources ---
# Swiss-Prot FASTA (downloaded once and reused)
SWISSPROT_FASTA="${OUTDIR}/databases/uniprot_sprot.fasta.gz"
DIAMOND_DB_PREFIX="${OUTDIR}/databases/swissprot.dmnd"

# dbCAN database directory (downloaded once and reused)
DBCAN_DB_DIR="${OUTDIR}/databases/dbcan_db"

# antiSMASH databases (downloaded once per machine/user)
# The antismash tool manages this internally via download-antismash-databases.

# BUSCO lineage
BUSCO_LINEAGE="hypocreaceae_odb12"

#==============================================================================
# INTERNAL PATHS (do not edit unless you know why)
#==============================================================================

BRAKER_DIR="${OUTDIR}/01_braker"
QC_DIR="${OUTDIR}/02_qc"
PROTEINS_DIR="${OUTDIR}/03_proteins"
DIAMOND_DIR="${OUTDIR}/04_diamond"
INTERPRO_DIR="${OUTDIR}/05_interproscan"
DBCAN_DIR="${OUTDIR}/06_dbcan"
ANTISMASH_DIR="${OUTDIR}/07_antismash"
DB_DIR="${OUTDIR}/databases"
LOG_DIR="${OUTDIR}/logs"

# Canonical outputs (authoritative gene set)
BRAKER_GFF3="${BRAKER_DIR}/braker.gff3"
BRAKER_AA="${BRAKER_DIR}/braker.aa"

# Canonical protein set used for ALL functional tools
PROT_ONELINE="${PROTEINS_DIR}/braker.oneline.faa"
PROT_LONGEST="${PROTEINS_DIR}/braker.longest.faa"
PROT_LONGEST_CLEAN="${PROTEINS_DIR}/braker.longest.clean.faa"

# Longest-isoform-only GFF3 (for antiSMASH mapping)
GFF_LONGEST="${PROTEINS_DIR}/braker.longest.gff3"
GFF_LONGEST_CLEAN="${PROTEINS_DIR}/braker.longest.clean.gff3"

#==============================================================================
# HELPER FUNCTIONS
#==============================================================================

log() {
  echo "[$(date -u +'%Y-%m-%dT%H:%M:%SZ')] $*"
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

require_file() {
  local f="$1"
  [[ -s "$f" ]] || die "Required file not found or empty: $f"
}

require_cmd() {
  local cmd="$1"
  command -v "$cmd" >/dev/null 2>&1 || die "Required command not found in PATH: $cmd"
}

run_logged() {
  # Usage: run_logged LOGFILE command args...
  local logfile="$1"
  shift
  mkdir -p "$(dirname "$logfile")"
  log "Running: $*"
  # shellcheck disable=SC2090
  "$@" 2>&1 | tee "$logfile"
}

ensure_dirs() {
  mkdir -p \
    "$BRAKER_DIR" "$QC_DIR" "$PROTEINS_DIR" "$DIAMOND_DIR" \
    "$INTERPRO_DIR" "$DBCAN_DIR" "$ANTISMASH_DIR" \
    "$DB_DIR" "$LOG_DIR"
}

#==============================================================================
# SANITY CHECKS
#==============================================================================

ensure_dirs
require_file "$GENOME_FASTA"

# mamba must exist for mamba run
require_cmd "mamba"
require_cmd "awk"
require_cmd "grep"
require_cmd "sed"
require_cmd "cut"
require_cmd "sort"
require_cmd "uniq"
require_cmd "wc"

log "Genome: $GENOME_FASTA"
log "Species tag: $SPECIES_TAG"
log "Threads: $THREADS"
log "Output directory: $OUTDIR"


### ================================
### PREFLIGHT: BRAKER3 requirements
### ================================

preflight_braker() {
  log "Preflight: checking BRAKER3 dependencies and inputs..."

  # ---- Input genome ----
  if [[ ! -s "$GENOME" ]]; then
    die "Genome FASTA not found or empty: $GENOME"
  fi

  # ---- Conda env sanity ----
  if ! mamba run -n braker braker.pl --help >/dev/null 2>&1; then
    die "Cannot run braker.pl inside env 'braker'. Is the environment installed correctly?"
  fi

  # ---- GeneMark ----
  if [[ -z "${GENEMARK_PATH:-}" ]]; then
    die "GENEMARK_PATH is not set. Example: export GENEMARK_PATH=/path/to/gmes_linux_64_4"
  fi

  if [[ ! -d "$GENEMARK_PATH" ]]; then
    die "GENEMARK_PATH does not exist or is not a directory: $GENEMARK_PATH"
  fi

  if [[ ! -x "$GENEMARK_PATH/gmes_petap.pl" ]]; then
    die "GeneMark executable not found or not executable: $GENEMARK_PATH/gmes_petap.pl"
  fi

  # GeneMark key (BRAKER/GeneMark requirement)
  if [[ ! -f "$HOME/.gm_key" ]]; then
    die "GeneMark key not found: $HOME/.gm_key (BRAKER requires this)"
  fi

  # Warn if key permissions are unsafe (GeneMark is picky)
  gm_perm="$(stat -c "%a" "$HOME/.gm_key" 2>/dev/null || true)"
  if [[ "$gm_perm" != "600" ]]; then
    log "WARNING: ~/.gm_key permissions are $gm_perm (recommended: 600). Fix with: chmod 600 ~/.gm_key"
  fi

  # ---- Augustus ----
  if ! mamba run -n braker which augustus >/dev/null 2>&1; then
    die "augustus not found inside env 'braker'. Install with: mamba install -n braker -c bioconda augustus"
  fi

  augustus_path="$(mamba run -n braker which augustus | tr -d '\r')"
  if [[ ! -x "$augustus_path" ]]; then
    die "augustus exists but is not executable: $augustus_path"
  fi

  # ---- Environment variables BRAKER expects ----
  # (not strictly required if conda layout is correct, but useful)
  if [[ -z "${AUGUSTUS_CONFIG_PATH:-}" ]]; then
    log "WARNING: AUGUSTUS_CONFIG_PATH is not set. BRAKER will try to guess it."
  fi

  log "Preflight OK: BRAKER3 inputs and dependencies look good."
}

preflight_braker

#==============================================================================
# LAYER 1 — GENE PREDICTION (BRAKER3 is authoritative)
#==============================================================================

log "=== Layer 1: Gene prediction (BRAKER3) ==="

if [[ ! -s "$BRAKER_GFF3" || ! -s "$BRAKER_AA" ]]; then
  run_logged "${LOG_DIR}/braker3.log" \
    mamba run -n "$ENV_BRAKER" braker.pl \
      --genome="$GENOME_FASTA" \
      --species="$SPECIES_TAG" \
      --esmode \
      --fungus \
      --softmasking \
      --threads="$THREADS" \
      --gff3 \
      --workingdir="$BRAKER_DIR"
else
  log "BRAKER outputs already exist. Skipping BRAKER3."
fi

require_file "$BRAKER_GFF3"
require_file "$BRAKER_AA"

# Quick sanity: count gene features
log "BRAKER gene feature count:"
grep -c -w "gene" "$BRAKER_GFF3" || true

#==============================================================================
# QC / BASIC STATS (non-destructive)
#==============================================================================

log "=== QC: AGAT stats + BUSCO ==="

# AGAT statistics
if [[ ! -s "${QC_DIR}/annotation_stats.txt" ]]; then
  run_logged "${LOG_DIR}/agat_stats.log" \
    mamba run -n "$ENV_AGAT" agat_sp_statistics.pl \
      --gff "$BRAKER_GFF3" \
      -o "${QC_DIR}/annotation_stats.txt"
else
  log "AGAT stats already exist. Skipping."
fi

# BUSCO on proteins
BUSCO_OUT="${QC_DIR}/busco"
if [[ ! -s "${BUSCO_OUT}/short_summary.txt" && ! -s "${BUSCO_OUT}/short_summary.specific.${BUSCO_LINEAGE}.txt" ]]; then
  mkdir -p "$BUSCO_OUT"
  run_logged "${LOG_DIR}/busco.log" \
    mamba run -n "$ENV_BUSCO" busco \
      -i "$BRAKER_AA" \
      -o "busco" \
      -m protein \
      -l "$BUSCO_LINEAGE" \
      -c "$THREADS" \
      --offline -f \
      --out_path "$QC_DIR"
else
  log "BUSCO output already exists. Skipping."
fi

#==============================================================================
# CANONICAL PROTEIN SET (Longest isoform per gene)
#==============================================================================

log "=== Canonical protein set: longest isoform per gene ==="
log "This protein FASTA will be used by ALL functional tools."

# Convert to one-line FASTA (stable parsing)
if [[ ! -s "$PROT_ONELINE" ]]; then
  run_logged "${LOG_DIR}/seqkit_oneline.log" \
    mamba run -n "$ENV_SEQKIT" seqkit seq -w 0 "$BRAKER_AA" > "$PROT_ONELINE"
else
  log "Oneline FASTA already exists. Skipping."
fi
require_file "$PROT_ONELINE"

# Extract longest isoform per gene (expects IDs like g1.t1, g1.t2, etc.)
BEST_IDS="${PROTEINS_DIR}/best_ids.txt"
if [[ ! -s "$BEST_IDS" ]]; then
  run_logged "${LOG_DIR}/select_longest_isoforms.log" \
    bash -c "
      mamba run -n '$ENV_SEQKIT' seqkit fx2tab -n -l '$PROT_ONELINE' \
      | awk '
          {
            id=\$1; len=\$2;
            gene=id;
            sub(/\\.t[0-9]+\$/, \"\", gene);

            if (!(gene in best) || len > best[gene]) {
              best[gene]=len;
              bestid[gene]=id;
            }
          }
          END {
            for (g in bestid) print bestid[g];
          }
        ' > '$BEST_IDS'
    "
else
  log "best_ids.txt already exists. Skipping."
fi
require_file "$BEST_IDS"

# Create longest isoform FASTA
if [[ ! -s "$PROT_LONGEST" ]]; then
  run_logged "${LOG_DIR}/seqkit_grep_longest.log" \
    mamba run -n "$ENV_SEQKIT" seqkit grep -f "$BEST_IDS" "$PROT_ONELINE" > "$PROT_LONGEST"
else
  log "Longest isoform FASTA already exists. Skipping."
fi
require_file "$PROT_LONGEST"

# Clean stop codons (*) for tools that dislike them
if [[ ! -s "$PROT_LONGEST_CLEAN" ]]; then
  run_logged "${LOG_DIR}/clean_stop_codons.log" \
    bash -c "sed 's/\\*//g' '$PROT_LONGEST' > '$PROT_LONGEST_CLEAN'"
else
  log "Clean protein FASTA already exists. Skipping."
fi
require_file "$PROT_LONGEST_CLEAN"

# Sanity counts
log "Protein counts:"
log "  BRAKER total proteins:   $(grep -c '^>' "$BRAKER_AA" || true)"
log "  Longest isoforms:        $(grep -c '^>' "$PROT_LONGEST" || true)"
log "  Longest isoforms (clean):$(grep -c '^>' "$PROT_LONGEST_CLEAN" || true)"

#==============================================================================
# LAYER 2 — FUNCTIONAL ANNOTATION (all tools annotate the same proteins)
#==============================================================================

#------------------------------------------------------------------------------
# DIAMOND + Swiss-Prot
#------------------------------------------------------------------------------
log "=== Layer 2: Functional annotation — DIAMOND + Swiss-Prot ==="

mkdir -p "$DB_DIR"

# Download Swiss-Prot once (if missing)
if [[ ! -s "$SWISSPROT_FASTA" ]]; then
  log "Swiss-Prot FASTA not found. Downloading..."
  run_logged "${LOG_DIR}/download_swissprot.log" \
    bash -c "wget -O '$SWISSPROT_FASTA' https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
else
  log "Swiss-Prot FASTA already exists. Skipping download."
fi
require_file "$SWISSPROT_FASTA"

# Build DIAMOND database (once)
if [[ ! -s "${DIAMOND_DB_PREFIX}.dmnd" ]]; then
  run_logged "${LOG_DIR}/diamond_makedb.log" \
    mamba run -n "$ENV_DIAMOND" diamond makedb \
      --in "$SWISSPROT_FASTA" \
      -d "$DIAMOND_DB_PREFIX"
else
  log "DIAMOND DB already exists. Skipping makedb."
fi
require_file "${DIAMOND_DB_PREFIX}.dmnd"

# Run blastp
DIAMOND_OUT="${DIAMOND_DIR}/swissprot_hits.tsv"
if [[ ! -s "$DIAMOND_OUT" ]]; then
  run_logged "${LOG_DIR}/diamond_blastp.log" \
    mamba run -n "$ENV_DIAMOND" diamond blastp \
      -d "$DIAMOND_DB_PREFIX" \
      -q "$PROT_LONGEST_CLEAN" \
      -o "$DIAMOND_OUT" \
      --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
      --max-target-seqs 1 \
      --evalue 1e-5 \
      --threads "$THREADS" \
      --sensitive
else
  log "DIAMOND output already exists. Skipping."
fi
require_file "$DIAMOND_OUT"

#------------------------------------------------------------------------------
# InterProScan
#------------------------------------------------------------------------------
log "=== Layer 2: Functional annotation — InterProScan ==="

INTERPRO_TSV="${INTERPRO_DIR}/interproscan.tsv"
INTERPRO_GFF3="${INTERPRO_DIR}/interproscan.gff3"

if [[ ! -s "$INTERPRO_TSV" ]]; then
  run_logged "${LOG_DIR}/interproscan.log" \
    mamba run -n "$ENV_INTERPROSCAN" interproscan.sh \
      -i "$PROT_LONGEST_CLEAN" \
      -appl TIGRFAM,FunFam,SFLD,SUPERFAMILY,PANTHER,Gene3D,Hamap,ProSiteProfiles,Coils,SMART,PRINTS,PIRSR,AntiFam,Pfam,MobiDBLite,PIRSF \
      -f TSV,GFF3 \
      -goterms \
      -iprlookup \
      -pa \
      -cpu "$THREADS" \
      -d "$INTERPRO_DIR"

  # InterProScan names output files based on input filename; normalize them.
  # This makes the pipeline deterministic for downstream integration.
  if [[ -s "${INTERPRO_DIR}/$(basename "$PROT_LONGEST_CLEAN").tsv" ]]; then
    mv -f "${INTERPRO_DIR}/$(basename "$PROT_LONGEST_CLEAN").tsv" "$INTERPRO_TSV"
  fi
  if [[ -s "${INTERPRO_DIR}/$(basename "$PROT_LONGEST_CLEAN").gff3" ]]; then
    mv -f "${INTERPRO_DIR}/$(basename "$PROT_LONGEST_CLEAN").gff3" "$INTERPRO_GFF3"
  fi
else
  log "InterProScan TSV already exists. Skipping."
fi
require_file "$INTERPRO_TSV"

# Basic sanity: how many unique proteins appear
log "InterProScan unique query IDs:"
cut -f1 "$INTERPRO_TSV" | sort -u | wc -l || true

#------------------------------------------------------------------------------
# dbCAN (CAZyme annotation)
#------------------------------------------------------------------------------
log "=== Layer 2: Functional annotation — dbCAN ==="

mkdir -p "$DBCAN_DIR" "$DBCAN_DB_DIR"

# Download dbCAN database once (if missing)
# Note: run_dbcan "database" step populates --db_dir.
if [[ ! -d "$DBCAN_DB_DIR" || -z "$(ls -A "$DBCAN_DB_DIR" 2>/dev/null || true)" ]]; then
  run_logged "${LOG_DIR}/dbcan_download_db.log" \
    mamba run -n "$ENV_DBCAN" run_dbcan database \
      --db_dir "$DBCAN_DB_DIR" \
      --aws_s3
else
  log "dbCAN DB directory already populated. Skipping download."
fi

# Run dbCAN protein mode
# dbCAN produces multiple output files in its output_dir.
DBCAN_DONE_FLAG="${DBCAN_DIR}/run_dbcan.done"
if [[ ! -s "$DBCAN_DONE_FLAG" ]]; then
  run_logged "${LOG_DIR}/dbcan_run.log" \
    mamba run -n "$ENV_DBCAN" run_dbcan CAZyme_annotation \
      --input_raw_data "$PROT_LONGEST_CLEAN" \
      --output_dir "$DBCAN_DIR" \
      --db_dir "$DBCAN_DB_DIR" \
      --mode protein \
      --threads "$THREADS"

  echo "done" > "$DBCAN_DONE_FLAG"
else
  log "dbCAN appears complete. Skipping."
fi

#------------------------------------------------------------------------------
# antiSMASH (secondary metabolite gene clusters)
#------------------------------------------------------------------------------
log "=== Layer 2: Functional annotation — antiSMASH ==="
log "antiSMASH is run with BRAKER genes (no gene finding)."

mkdir -p "$ANTISMASH_DIR"

# Ensure antiSMASH databases exist
ANTISMASH_DB_FLAG="${OUTDIR}/databases/antismash_db_downloaded.flag"
if [[ ! -s "$ANTISMASH_DB_FLAG" ]]; then
  run_logged "${LOG_DIR}/antismash_download_db.log" \
    mamba run -n "$ENV_ANTISMASH" download-antismash-databases
  echo "done" > "$ANTISMASH_DB_FLAG"
else
  log "antiSMASH databases already downloaded. Skipping."
fi

# Build longest-isoform-only GFF3 for antiSMASH mapping
# We take the protein IDs from the longest FASTA and keep only matching GFF3 lines.
LONGEST_IDS="${PROTEINS_DIR}/longest_isoform_ids.txt"

if [[ ! -s "$LONGEST_IDS" ]]; then
  run_logged "${LOG_DIR}/extract_longest_ids.log" \
    bash -c "grep '^>' '$PROT_LONGEST_CLEAN' | sed 's/^>//' | sed 's/ .*//' > '$LONGEST_IDS'"
else
  log "Longest isoform IDs already exist. Skipping."
fi
require_file "$LONGEST_IDS"

if [[ ! -s "$GFF_LONGEST" ]]; then
  run_logged "${LOG_DIR}/filter_gff_longest.log" \
    bash -c "grep -F -f '$LONGEST_IDS' '$BRAKER_GFF3' > '$GFF_LONGEST'"
else
  log "Longest-only GFF3 already exists. Skipping."
fi
require_file "$GFF_LONGEST"

if [[ ! -s "$GFF_LONGEST_CLEAN" ]]; then
  run_logged "${LOG_DIR}/write_clean_gff3.log" \
    bash -c "printf '##gff-version 3\n' > '$GFF_LONGEST_CLEAN' && cat '$GFF_LONGEST' >> '$GFF_LONGEST_CLEAN'"
else
  log "Clean longest-only GFF3 already exists. Skipping."
fi
require_file "$GFF_LONGEST_CLEAN"

# Run antiSMASH
ANTISMASH_DONE_FLAG="${ANTISMASH_DIR}/antismash.done"
if [[ ! -s "$ANTISMASH_DONE_FLAG" ]]; then
  run_logged "${LOG_DIR}/antismash.log" \
    mamba run -n "$ENV_ANTISMASH" antismash \
      --taxon fungi \
      --output-dir "$ANTISMASH_DIR" \
      --genefinding-tool none \
      --genefinding-gff3 "$GFF_LONGEST_CLEAN" \
      --cpus "$THREADS" \
      "$GENOME_FASTA" \
      --fullhmmer \
      --clusterhmmer \
      --tigrfam \
      --asf \
      --cb-general \
      --pfam2go \
      --tfbs \
      --smcog-trees

  echo "done" > "$ANTISMASH_DONE_FLAG"
else
  log "antiSMASH appears complete. Skipping."
fi

#==============================================================================
# LAYER 3 — ANNOTATION INTEGRATION (outputs prepared for merging)
#==============================================================================

log "=== Layer 3: Annotation integration (preparation) ==="
log "This script prepares standardized outputs for merging by gene/protein ID."
log "Integration itself (parsing+merging into a master table) is intentionally"
log "left as a separate step/module, so the pipeline remains clean and modular."

# Write a manifest of canonical files for downstream integration
MANIFEST="${OUTDIR}/annotation_manifest.tsv"
cat > "$MANIFEST" <<EOF
#key	path	description
genome_fasta	${GENOME_FASTA}	Genome assembly FASTA used for BRAKER3 and antiSMASH
braker_gff3	${BRAKER_GFF3}	Authoritative gene models (GFF3)
braker_proteins_all	${BRAKER_AA}	All predicted proteins from BRAKER3
proteins_longest	${PROT_LONGEST}	Longest isoform per gene (canonical protein set)
proteins_longest_clean	${PROT_LONGEST_CLEAN}	Longest isoform proteins with '*' removed
gff_longest_clean	${GFF_LONGEST_CLEAN}	GFF3 filtered to longest isoforms (for antiSMASH mapping)
diamond_swissprot_hits	${DIAMOND_OUT}	DIAMOND blastp vs Swiss-Prot (TSV outfmt6)
interproscan_tsv	${INTERPRO_TSV}	InterProScan results (TSV)
interproscan_gff3	${INTERPRO_GFF3}	InterProScan annotations (GFF3)
dbcan_dir	${DBCAN_DIR}	dbCAN output directory (multiple files)
antismash_dir	${ANTISMASH_DIR}	antiSMASH output directory (clusters and gene mapping)
EOF

log "Manifest written: $MANIFEST"
log "Pipeline complete."
