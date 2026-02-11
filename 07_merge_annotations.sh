#!/usr/bin/env bash
#==============================================================================
# ANNOTATION INTEGRATION MODULE
#==============================================================================
# Purpose:
#   Merge functional annotations into a unified master table and a GFF3
#   attribute-augmented file using BRAKER3 as authoritative gene set.
#
# Inputs (from pipeline):
#   - BRAKER GFF3 (authoritative coordinates and IDs)
#   - Longest isoform protein FASTA (canonical protein set)
#   - DIAMOND Swiss-Prot hits (TSV)
#   - InterProScan TSV
#   - dbCAN outputs
#   - antiSMASH outputs (genecluster mapping)
#
# Outputs:
#   - master_annotation.tsv
#   - master_annotation.gff3
#
#==============================================================================

set -euo pipefail

#==============================================================================
# CONFIGURATION
#==============================================================================

OUTDIR="annotation_pipeline"

BRAKER_GFF3="${OUTDIR}/01_braker/braker.gff3"
PROT_LONGEST_CLEAN="${OUTDIR}/03_proteins/braker.longest.clean.faa"

DIAMOND_TSV="${OUTDIR}/04_diamond/swissprot_hits.tsv"
INTERPRO_TSV="${OUTDIR}/05_interproscan/interproscan.tsv"
DBCAN_DIR="${OUTDIR}/06_dbcan"
ANTISMASH_DIR="${OUTDIR}/07_antismash"

INTEGRATION_DIR="${OUTDIR}/08_integration"
LOG_DIR="${OUTDIR}/logs"

THREADS="12"  # reserved for future parallel parsing

# Final outputs
MASTER_TSV="${INTEGRATION_DIR}/master_annotation.tsv"
MASTER_GFF3="${INTEGRATION_DIR}/master_annotation.gff3"

#==============================================================================
# HELPERS
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

#==============================================================================
# SANITY CHECKS
#==============================================================================

mkdir -p "$INTEGRATION_DIR" "$LOG_DIR"

require_file "$BRAKER_GFF3"
require_file "$PROT_LONGEST_CLEAN"
require_file "$DIAMOND_TSV"
require_file "$INTERPRO_TSV"

# dbCAN: expected core outputs (depends on dbCAN version; these are common)
DBCAN_OVERVIEW="${DBCAN_DIR}/overview.txt"
DBCAN_HMMOUT="${DBCAN_DIR}/hmm.out"
DBCAN_DIAMONDOUT="${DBCAN_DIR}/diamond.out"
DBCAN_HOTPEP="${DBCAN_DIR}/Hotpep.out"

if [[ ! -s "$DBCAN_OVERVIEW" ]]; then
  log "WARNING: dbCAN overview.txt not found. Will try to parse other dbCAN outputs."
fi

# antiSMASH: there are multiple possible files depending on version.
# We will parse genecluster mappings if present.
ANTISMASH_GENECLUSTERS="${ANTISMASH_DIR}/geneclusters.txt"
if [[ ! -s "$ANTISMASH_GENECLUSTERS" ]]; then
  log "WARNING: antiSMASH geneclusters.txt not found. antiSMASH cluster mapping will be empty."
fi

log "Integration inputs validated."

#==============================================================================
# STEP 1 — BUILD THE AUTHORITATIVE GENE/PROTEIN BACKBONE TABLE
#==============================================================================

log "Step 1: building backbone table from BRAKER GFF3 + longest proteins"

# Extract canonical protein IDs from FASTA
# This defines the authoritative row set for the master table.
PROT_IDS="${INTEGRATION_DIR}/protein_ids.txt"
grep '^>' "$PROT_LONGEST_CLEAN" | sed 's/^>//' | sed 's/ .*//' | sort -u > "$PROT_IDS"
require_file "$PROT_IDS"

# Build a coordinate table from BRAKER GFF3
# We want a per-mRNA or per-protein entry, but BRAKER GFF3 can vary.
# Strategy:
#   - Use CDS features: derive min(start), max(end) per Parent (transcript)
#   - Keep seqid + strand
#
# This is robust across different GFF3 flavors.
COORD_TSV="${INTEGRATION_DIR}/coords.tsv"

awk -F'\t' '
  BEGIN { OFS="\t" }
  $0 ~ /^#/ { next }
  NF < 9 { next }
  $3 != "CDS" { next }

  {
    seqid=$1; start=$4; end=$5; strand=$7; attr=$9;

    parent="";
    if (match(attr, /Parent=[^;]+/)) {
      parent=substr(attr, RSTART+7, RLENGTH-7);
      # If multiple parents, take first
      sub(/,.*/, "", parent);
    } else {
      next;
    }

    # accumulate min/max coords per transcript/protein ID
    if (!(parent in minS) || start < minS[parent]) minS[parent]=start;
    if (!(parent in maxE) || end > maxE[parent]) maxE[parent]=end;

    seq[parent]=seqid;
    str[parent]=strand;
  }
  END {
    for (p in minS) {
      print p, seq[p], minS[p], maxE[p], str[p];
    }
  }
' "$BRAKER_GFF3" | sort -k1,1 > "$COORD_TSV"

# Filter coords to only longest proteins (canonical set)
COORD_LONGEST_TSV="${INTEGRATION_DIR}/coords.longest.tsv"
join -t $'\t' -1 1 -2 1 \
  <(sort -k1,1 "$PROT_IDS") \
  <(sort -k1,1 "$COORD_TSV") \
  > "$COORD_LONGEST_TSV" || true

# If join produced nothing, we fall back to empty coords but keep pipeline alive.
if [[ ! -s "$COORD_LONGEST_TSV" ]]; then
  log "WARNING: Could not map protein IDs to CDS Parent IDs in GFF3."
  log "         Coordinates will be empty in the master table."
fi

#==============================================================================
# STEP 2 — PARSE DIAMOND (Swiss-Prot) INTO NORMALIZED TABLE
#==============================================================================

log "Step 2: parsing DIAMOND Swiss-Prot hits"

DIAMOND_PARSED="${INTEGRATION_DIR}/diamond.parsed.tsv"

# Input format:
# qseqid sseqid pident length evalue bitscore stitle
#
# Output:
# gene_id swissprot_accession swissprot_pident swissprot_alnlen swissprot_evalue swissprot_bitscore swissprot_title
awk -F'\t' 'BEGIN{OFS="\t"}
  NF < 7 { next }
  {
    q=$1; s=$2; pident=$3; alnlen=$4; e=$5; bits=$6;

    # stitle may contain tabs in some edge cases; reconstruct safely
    title=$7;
    if (NF > 7) {
      for (i=8; i<=NF; i++) title=title "\t" $i;
    }

    print q, s, pident, alnlen, e, bits, title;
  }
' "$DIAMOND_TSV" | sort -k1,1 > "$DIAMOND_PARSED"

#==============================================================================
# STEP 3 — PARSE INTERPROSCAN TSV INTO PER-PROTEIN FIELDS
#==============================================================================

log "Step 3: parsing InterProScan TSV (domains + GO)"

INTERPRO_PARSED="${INTEGRATION_DIR}/interpro.parsed.tsv"

# InterProScan TSV columns (common):
# 1 protein_accession
# 4 signature_accession
# 5 signature_description
# 9 InterPro accession
# 10 InterPro description
# 14 GO terms
#
# We aggregate per protein:
#   - interpro_accessions (unique; ;)
#   - interpro_descriptions (unique; ;)
#   - signatures (unique; ;)
#   - go_terms (unique; ;)
#
# Note: Many rows have empty GO/IPR fields.
awk -F'\t' '
  BEGIN { OFS="\t" }
  {
    prot=$1;

    sig_acc=$4;
    sig_desc=$5;

    ipr_acc=$12;  # In many versions InterPro accession is col 12
    ipr_desc=$13; # InterPro description col 13
    go=$14;

    # Some installations shift columns depending on options.
    # Fallback heuristics: if ipr_acc does not look like IPRxxxxxx, search for it.
    if (ipr_acc !~ /^IPR[0-9]+$/) {
      ipr_acc="";
      ipr_desc="";
      for (i=1; i<=NF; i++) {
        if ($i ~ /^IPR[0-9]+$/) {
          ipr_acc=$i;
          # description likely next column
          if (i+1 <= NF) ipr_desc=$(i+1);
          break;
        }
      }
    }

    # GO terms usually look like GO:0000000;GO:....
    if (go !~ /GO:/) {
      # search last columns for GO
      for (i=NF; i>=1; i--) {
        if ($i ~ /GO:/) { go=$i; break; }
      }
    }

    # store unique sets
    if (sig_acc != "" && sig_acc != "-") sigs[prot][sig_acc]=1;
    if (ipr_acc != "" && ipr_acc != "-") iprs[prot][ipr_acc]=1;
    if (ipr_desc != "" && ipr_desc != "-") iprdescs[prot][ipr_desc]=1;

    if (go ~ /GO:/) {
      n=split(go, arr, /\|/);
      for (j=1; j<=n; j++) {
        if (arr[j] ~ /GO:/) gos[prot][arr[j]]=1;
      }
      n=split(go, arr2, /,/);
      for (j=1; j<=n; j++) {
        if (arr2[j] ~ /GO:/) gos[prot][arr2[j]]=1;
      }
      n=split(go, arr3, /;/);
      for (j=1; j<=n; j++) {
        if (arr3[j] ~ /GO:/) gos[prot][arr3[j]]=1;
      }
    }
  }

  function join_set(set,   k, out, first) {
    first=1;
    out="";
    for (k in set) {
      if (first) { out=k; first=0 }
      else { out=out ";" k }
    }
    return out;
  }

  END {
    # print per protein
    for (p in sigs) {
      s = join_set(sigs[p]);
      i = join_set(iprs[p]);
      d = join_set(iprdescs[p]);
      g = join_set(gos[p]);
      print p, i, d, s, g;
    }
  }
' "$INTERPRO_TSV" | sort -k1,1 > "$INTERPRO_PARSED"

#==============================================================================
# STEP 4 — PARSE dbCAN INTO PER-PROTEIN CAZy FAMILIES
#==============================================================================

log "Step 4: parsing dbCAN (CAZy families)"

DBCAN_PARSED="${INTEGRATION_DIR}/dbcan.parsed.tsv"

# Best dbCAN summary is overview.txt (one line per query)
# If not present, parsing is not reliable.
if [[ -s "$DBCAN_OVERVIEW" ]]; then
  # overview.txt usually contains columns including:
  # Gene ID, HMMER, DIAMOND, Hotpep, ... , CAZy families
  #
  # We extract:
  #   gene_id   cazy_families
  #
  # But dbCAN formats vary. A safe approach:
  #   - first field = gene ID
  #   - find tokens that look like CAZy families (GHxx, GTxx, CE, AA, PL, CBM)
  #
  awk '
    BEGIN { OFS="\t" }
    NR==1 { next } # skip header
    {
      gene=$1;
      fams="";

      for (i=1; i<=NF; i++) {
        if ($i ~ /^(GH|GT|CE|AA|PL|CBM)[0-9]+/ || $i ~ /^(GH|GT|CE|AA|PL|CBM)$/) {
          if (fams=="") fams=$i;
          else fams=fams ";" $i;
        }
      }

      # clean duplicates by set
      n=split(fams, a, /;/);
      delete seen;
      out="";
      for (j=1; j<=n; j++) {
        if (a[j]=="" || a[j]=="-") continue;
        if (!(a[j] in seen)) {
          seen[a[j]]=1;
          if (out=="") out=a[j];
          else out=out ";" a[j];
        }
      }

      print gene, out;
    }
  ' "$DBCAN_OVERVIEW" | sort -k1,1 > "$DBCAN_PARSED"
else
  log "WARNING: dbCAN overview.txt missing. Writing empty dbCAN table."
  : > "$DBCAN_PARSED"
fi

#==============================================================================
# STEP 5 — PARSE antiSMASH GENECLUSTERS INTO PER-GENE CLUSTER FIELDS
#==============================================================================

log "Step 5: parsing antiSMASH genecluster mapping"

ANTISMASH_PARSED="${INTEGRATION_DIR}/antismash.parsed.tsv"

# geneclusters.txt typically maps gene IDs to cluster numbers/types.
# Formats vary, so we use a conservative parser:
#   - If a line contains a gene ID from our protein list, capture cluster info.
#
# Output:
# gene_id antismash_cluster_ids antismash_cluster_types
#
if [[ -s "$ANTISMASH_GENECLUSTERS" ]]; then
  # Build a fast lookup of valid protein IDs
  awk '
    NR==FNR { ids[$1]=1; next }
    {
      line=$0;

      # attempt to extract a gene/protein ID token
      # Many antiSMASH versions include gene IDs explicitly in the file.
      for (i=1; i<=NF; i++) {
        if ($i in ids) {
          gene=$i;
          # heuristic: cluster id and type are often near the start of the line
          # Example patterns: "Cluster 1" or "cluster_1" or "region001"
          cid="";
          ctype="";

          if (match(line, /(region[0-9]+|Region[0-9]+|cluster[ _-]?[0-9]+)/)) {
            cid=substr(line, RSTART, RLENGTH);
          }
          if (match(line, /(NRPS|PKS|T1PKS|T2PKS|T3PKS|terpene|RiPP|betalactone|siderophore|indole|fungal)/i)) {
            ctype=substr(line, RSTART, RLENGTH);
          }

          # store sets
          if (cid != "") clid[gene][cid]=1;
          if (ctype != "") cltype[gene][ctype]=1;
        }
      }
    }

    function join_set(set,   k, out, first) {
      first=1;
      out="";
      for (k in set) {
        if (first) { out=k; first=0 }
        else { out=out ";" k }
      }
      return out;
    }

    END {
      for (g in ids) {
        # only print those that were seen in antiSMASH
        if (g in clid || g in cltype) {
          print g, join_set(clid[g]), join_set(cltype[g]);
        }
      }
    }
  ' "$PROT_IDS" "$ANTISMASH_GENECLUSTERS" | sort -k1,1 > "$ANTISMASH_PARSED"
else
  log "WARNING: antiSMASH geneclusters.txt missing. Writing empty antiSMASH table."
  : > "$ANTISMASH_PARSED"
fi

#==============================================================================
# STEP 6 — MERGE ALL TABLES INTO MASTER ANNOTATION TSV
#==============================================================================

log "Step 6: merging into master_annotation.tsv"

# Master schema:
# protein_id  contig  start  end  strand
# swissprot_accession swissprot_pident swissprot_alnlen swissprot_evalue swissprot_bitscore swissprot_title
# interpro_accessions interpro_descriptions signatures go_terms
# cazy_families
# antismash_cluster_ids antismash_cluster_types
#
# We perform a controlled merge using awk with hash maps.
#
# IMPORTANT:
#   - Backbone = canonical protein IDs from longest FASTA
#   - Missing values = empty

awk -F'\t' -v OFS='\t' \
  -v coord_file="$COORD_LONGEST_TSV" \
  -v diamond_file="$DIAMOND_PARSED" \
  -v interpro_file="$INTERPRO_PARSED" \
  -v dbcan_file="$DBCAN_PARSED" \
  -v antismash_file="$ANTISMASH_PARSED" \
  '
  function load_coords(file,   line,a,id) {
    while ((getline line < file) > 0) {
      split(line, a, "\t");
      id=a[1];
      coord_contig[id]=a[2];
      coord_start[id]=a[3];
      coord_end[id]=a[4];
      coord_strand[id]=a[5];
    }
    close(file);
  }

  function load_diamond(file,   line,a,id) {
    while ((getline line < file) > 0) {
      split(line, a, "\t");
      id=a[1];
      dia_acc[id]=a[2];
      dia_pident[id]=a[3];
      dia_alnlen[id]=a[4];
      dia_eval[id]=a[5];
      dia_bits[id]=a[6];

      # title may contain tabs; reconstruct from original line:
      # In our parsed file, it should be 7 columns, but be defensive.
      dia_title[id]=a[7];
      if (length(a) > 7) {
        for (i=8; i in a; i++) dia_title[id]=dia_title[id] "\t" a[i];
      }
    }
    close(file);
  }

  function load_interpro(file,   line,a,id) {
    while ((getline line < file) > 0) {
      split(line, a, "\t");
      id=a[1];
      ipr_acc[id]=a[2];
      ipr_desc[id]=a[3];
      ipr_sig[id]=a[4];
      ipr_go[id]=a[5];
    }
    close(file);
  }

  function load_dbcan(file,   line,a,id) {
    while ((getline line < file) > 0) {
      split(line, a, "\t");
      id=a[1];
      cazy[id]=a[2];
    }
    close(file);
  }

  function load_antismash(file,   line,a,id) {
    while ((getline line < file) > 0) {
      split(line, a, "\t");
      id=a[1];
      as_cid[id]=a[2];
      as_type[id]=a[3];
    }
    close(file);
  }

  BEGIN {
    load_coords(coord_file);
    load_diamond(diamond_file);
    load_interpro(interpro_file);
    load_dbcan(dbcan_file);
    load_antismash(antismash_file);

    print \
      "protein_id","contig","start","end","strand", \
      "swissprot_accession","swissprot_pident","swissprot_alnlen","swissprot_evalue","swissprot_bitscore","swissprot_title", \
      "interpro_accessions","interpro_descriptions","signatures","go_terms", \
      "cazy_families", \
      "antismash_cluster_ids","antismash_cluster_types";
  }

  # Read canonical protein IDs from stdin
  {
    id=$1;

    print \
      id, \
      coord_contig[id], coord_start[id], coord_end[id], coord_strand[id], \
      dia_acc[id], dia_pident[id], dia_alnlen[id], dia_eval[id], dia_bits[id], dia_title[id], \
      ipr_acc[id], ipr_desc[id], ipr_sig[id], ipr_go[id], \
      cazy[id], \
      as_cid[id], as_type[id];
  }
' "$PROT_IDS" > "$MASTER_TSV"

require_file "$MASTER_TSV"
log "Master table written: $MASTER_TSV"

#==============================================================================
# STEP 7 — WRITE FUNCTIONALIZED GFF3 (BRAKER + attributes)
#==============================================================================

log "Step 7: writing master_annotation.gff3"

# We append attributes only to mRNA/transcript-like lines (Parent features),
# because that is where protein IDs are most consistently present.
#
# We attach:
#   SwissProt=...
#   SwissProtTitle=...
#   InterPro=...
#   GO=...
#   CAZy=...
#   antiSMASH=...
#
# This is conservative and keeps the original BRAKER GFF3 intact otherwise.

awk -F'\t' -v OFS='\t' \
  -v diamond_file="$DIAMOND_PARSED" \
  -v interpro_file="$INTERPRO_PARSED" \
  -v dbcan_file="$DBCAN_PARSED" \
  -v antismash_file="$ANTISMASH_PARSED" \
  '
  function load_simple(file, keycol, valcol, arr,   line,a) {
    while ((getline line < file) > 0) {
      split(line, a, "\t");
      if (a[keycol] != "") arr[a[keycol]] = a[valcol];
    }
    close(file);
  }

  function load_diamond(file,   line,a,id) {
    while ((getline line < file) > 0) {
      split(line, a, "\t");
      id=a[1];
      dia_acc[id]=a[2];
      dia_title[id]=a[7];
    }
    close(file);
  }

  function load_interpro(file,   line,a,id) {
    while ((getline line < file) > 0) {
      split(line, a, "\t");
      id=a[1];
      ipr_acc[id]=a[2];
      ipr_go[id]=a[5];
    }
    close(file);
  }

  function load_dbcan(file,   line,a,id) {
    while ((getline line < file) > 0) {
      split(line, a, "\t");
      id=a[1];
      cazy[id]=a[2];
    }
    close(file);
  }

  function load_antismash(file,   line,a,id) {
    while ((getline line < file) > 0) {
      split(line, a, "\t");
      id=a[1];
      as_cid[id]=a[2];
      as_type[id]=a[3];
    }
    close(file);
  }

  BEGIN {
    load_diamond(diamond_file);
    load_interpro(interpro_file);
    load_dbcan(dbcan_file);
    load_antismash(antismash_file);
  }

  $0 ~ /^#/ { print; next }
  NF < 9 { print; next }

  {
    attr=$9;

    # Identify an ID to attach to.
    # Prefer ID=..., otherwise Parent=...
    id="";

    if (match(attr, /ID=[^;]+/)) {
      id=substr(attr, RSTART+3, RLENGTH-3);
    } else if (match(attr, /Parent=[^;]+/)) {
      id=substr(attr, RSTART+7, RLENGTH-7);
      sub(/,.*/, "", id);
    }

    # Only annotate if this ID is present in our maps
    # (typically mRNA IDs match protein IDs for BRAKER, but not always)
    extra="";

    if (id != "") {
      if (id in dia_acc && dia_acc[id] != "") {
        extra=extra ";SwissProt=" dia_acc[id];
      }
      if (id in dia_title && dia_title[id] != "") {
        t=dia_title[id];
        gsub(/;/, ",", t);
        extra=extra ";SwissProtTitle=" t;
      }
      if (id in ipr_acc && ipr_acc[id] != "") {
        extra=extra ";InterPro=" ipr_acc[id];
      }
      if (id in ipr_go && ipr_go[id] != "") {
        extra=extra ";GO=" ipr_go[id];
      }
      if (id in cazy && cazy[id] != "") {
        extra=extra ";CAZy=" cazy[id];
      }
      if (id in as_cid && as_cid[id] != "") {
        extra=extra ";antiSMASHCluster=" as_cid[id];
      }
      if (id in as_type && as_type[id] != "") {
        extra=extra ";antiSMASHType=" as_type[id];
      }
    }

    $9 = attr extra;
    print;
  }
' "$BRAKER_GFF3" > "$MASTER_GFF3"

require_file "$MASTER_GFF3"
log "Master GFF3 written: $MASTER_GFF3"

#==============================================================================
# FINAL REPORT
#==============================================================================

log "Integration complete."
log "Outputs:"
log "  - $MASTER_TSV"
log "  - $MASTER_GFF3"
