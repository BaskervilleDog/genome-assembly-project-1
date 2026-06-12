# 🧬 Hybrid Genome Assembly Pipeline — *Trichoderma harzianum* TW11

A reproducible, multi-assembler pipeline for *de novo* genome assembly using PacBio long reads and Illumina short-read polishing, developed for the fungal *Trichoderma harzianum* TW11.

---

## Overview

This project implements a complete hybrid genome assembly workflow that:

1. **Downloads** raw sequencing data from SRA (NCBI) via SeqFetcher
2. **Quality controls** both PacBio (SeqKit, LongReadSum) and Illumina (FastQC, MultiQC) reads
3. **Filters** long reads with fastplong and short reads with fastp
4. **Decontaminates** reads with Kraken2 (taxonomic k-mer classification)
5. **Assembles** the genome with three independent assemblers (Flye, Raven, wtdbg2)
6. **Evaluates** raw assemblies with BUSCO + QUAST, then selects the best (Flye)
7. **Polishes** the selected assembly sequentially: Racon (×2 long-read) → Polypolish (short-read) → NextPolish (×2 short-read)
8. **Scaffolds** the polished assembly against a chromosome-level reference with RagTag
9. **Masks** repeats with RepeatModeler2 + RepeatMasker, producing a soft-masked assembly ready for structural annotation

Using multiple assemblers mitigates algorithm-specific biases and allows objective selection based on BUSCO completeness and QUAST contiguity metrics.

---

## Biological Context

| Feature | Detail |
|---|---|
| **Organism** | *Trichoderma harzianum* TW11 |
| **Application** | Fungal biocontrol agent genomics |
| **Long reads** | PacBio RS2/Sequel (SRR10848482) — 7.1 Gb, ~178× coverage |
| **Short reads** | Illumina HiSeq 4000 PE (SRR10848483 + SRR10848484) — ~61× combined |
| **Estimated genome size** | ~40 Mb |
| **BUSCO lineage** | `hypocreaceae_odb12` |

---

## Final Assembly Statistics

| Metric | Value |
|---|---|
| Total length | 38.95 Mb |
| Number of contigs | 70 |
| Number of scaffolds | 17 |
| Scaffold N50 | 1 Mbp |
| Gap content | 0.000% |
| BUSCO completeness | 99.8% (`hypocreaceae_odb12`) |
| Merqury QV | 34.4 |

---

## Repository Structure

```
project/
├── 00_downloads/               # Raw FASTQ files from SRA
├── 01_qc_raw/                  # SeqKit stats, LongReadSum, FastQC, MultiQC on raw reads
├── 02_filtered/
│   ├── pacbio/                 # fastplong-filtered long reads
│   └── illumina/               # fastp-filtered short reads
├── 03_qc_filtered/             # SeqKit stats, FastQC, MultiQC on filtered reads
├── 04_contaminants/            # Kraken2 classification reports and classified reads
├── 05_decontaminated/          # Unclassified (host-free) reads
├── 06_assemblies/              # Raw assembler outputs (flye/, raven/, wtdbg2/)
├── 07_assembly_evaluation/
│   ├── raw_assembly/           # BUSCO, QUAST, Merqury on pre-polishing assemblies
│   └── polished_assembly/      # BUSCO, QUAST, Merqury on final assembly
├── 08_polishing/
│   └── longreads/              # PAF alignments, Racon rounds, Polypolish SAMs, NextPolish
├── 09_ragtag/                  # RagTag scaffold output and evaluation
├── 10_repeat_masking/
│   ├── repeatmodeler/          # De novo repeat library
│   └── repeatmasker/           # Soft-masked assembly
└── logs/                       # All stdout/stderr logs
```

---

## Pipeline Stages

### Stage 1 — Data Acquisition

SeqFetcher wraps NCBI's `fasterq-dump` and related utilities to download all three SRA accessions in a single command.

```bash
cat > accession_list.txt << EOF
SRR10848482
SRR10848483
SRR10848484
EOF

mamba run -n ncbi_tools seqfetcher download \
    --sra-method parallel \
    --sra-accession-file accession_list.txt
```

### Stage 2 — Quality Control (Raw Reads)

Quality control differs by sequencing technology:

- **SeqKit stats** (`-a`): fast per-file summary statistics for all libraries
- **LongReadSum**: read-length histograms and N50 distribution plots for PacBio (quality scores absent from SRA)
- **FastQC**: per-base quality, GC content, and adapter detection for Illumina reads
- **MultiQC**: aggregates all FastQC and fastp reports into a single interactive HTML report

> **Note:** SeqKit reports Q=0 for PacBio reads from SRA — quality scores are stripped during SRA archiving. Read-length distribution (LongReadSum) is the primary QC metric for long reads.

### Stage 3 — Read Filtering

| Technology | Tool | Key parameters |
|---|---|---|
| PacBio | fastplong | `-Q` (disable quality filtering), `--length_required 1000` |
| Illumina | fastp | `--detect_adapter_for_pe`, `--qualified_quality_phred 30`, `--length_required 150`, `--correction` |

Quality filtering is disabled for PacBio because quality scores are absent. Length filtering (≥1 kb) removes reads too short to contribute to assembly contiguity. fastp's `--correction` flag uses read-pair overlaps to correct Illumina errors before polishing.

### Stage 4 — Decontamination

Kraken2 classifies reads against the standard database (bacteria, archaea, viruses, human). Reads marked **unclassified** are retained as the clean dataset; classified reads are saved separately for inspection. This step prevents chimeric contigs and inflated genome size estimates from bacterial contamination common in soil-derived fungal cultures.

### Stage 5 — Genome Assembly

Three assemblers with fundamentally different graph approaches are run on the same decontaminated PacBio reads:

```
Decontaminated PacBio reads
        │
        ├──────────────────────────────────────┐
        │                                      │
     [Flye]                   [Raven]      [wtdbg2]
   Repeat graph           String graph    Fuzzy dBG
        │                      │              │
        └──────────────┬────────┘──────────────┘
                       ▼
           BUSCO + QUAST evaluation
                       ▼
             Select best assembly (Flye)
```

| Assembler | Algorithm | Strength |
|---|---|---|
| Flye | Repeat graph | Handles tandem repeats; designed for high PacBio coverage |
| Raven | String/OLC graph | Fast; internal Racon polishing; high contiguity |
| wtdbg2 | Fuzzy de Bruijn graph | Structurally distinct approach; good comparison baseline |

> **Why three assemblers?** Each uses different overlap detection strategies and repeat handling heuristics, meaning they produce different assemblies from identical input. Empirical BUSCO/QUAST evaluation — not theoretical assumptions — determines the best result.

### Stage 6 — Assembly Polishing

Long-read assemblers produce structurally accurate contigs but with residual indel and substitution errors. Polishing proceeds in stages: long-read correction first (structural errors), then short-read correction (base-level errors).

```
Flye assembly
      │
      ▼
minimap2 + Racon × 2   (long-read consensus correction)
      │
      ▼
BWA-MEM (-a) + Polypolish   (all-alignment short-read polishing)
      │
      ▼
NextPolish × 2   (final SNP/indel correction, homopolymer-aware)
```

**Racon (×2):** Two rounds of consensus polishing across all aligned reads correct the largest errors first, then residual errors from round 1.

**Polypolish:** Uses *all* alignments (`bwa mem -a`) rather than uniquely mapping reads, enabling polishing of repetitive regions that Pilon skips. Both Illumina runs (SRR10848483 and SRR10848484) are aligned separately and jointly applied.

**NextPolish (×2):** More accurate than Pilon for indel correction in homopolymer runs — a known PacBio error source. Applied as final refinement with `task = best`.

### Stage 7 — Assembly Evaluation

Three complementary tools are applied at two stages (raw assembler outputs and final polished assembly):

| Tool | What it measures | When applied |
|---|---|---|
| BUSCO (`hypocreaceae_odb12`, 4,323 genes) | Gene space completeness | Raw assemblies + polished |
| QUAST | Structural contiguity (N50, gaps, misassemblies) | Raw assemblies + polished + reference |
| Merqury | Reference-free base accuracy (QV) | Polished assembly only |

**Final polished assembly BUSCO result:**
```
C:99.8% [S:99.5%, D:0.3%], F:0.0%, M:0.1%, n:4323, E:1.6%
```

**Merqury QV = 34.4** (~1 error per 2,800 bp). The 1.6% internal stop codon rate (E:1.6%) reflects residual frameshifts from indel errors, consistent with QV 34.4. This likely represents the ceiling achievable with the available Illumina coverage depth rather than a fixable error.

### Stage 8 — Reference-Guided Scaffolding

RagTag orders, orients, and joins the 70 polished contigs into pseudochromosomes using a chromosome-level *T. harzianum* reference ([GCA_019097725.1, SYAU_Tha_1.0](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_019097725.1/)).

The `-u` flag retains all unplaced contigs in the output, preserving any TW11-specific sequences absent from the reference strain.

> **Scaffolding vs. polishing order:** Scaffolding is performed *after* polishing and *before* repeat masking. RagTag's internal minimap2 alignment is inhibited by soft-masked bases, so masking must always be the final step before annotation.

**Key outputs:**
- `ragtag.scaffold.fasta` — scaffolded pseudochromosome assembly
- `ragtag.scaffold.agp` — contig placement and orientation record
- `ragtag.scaffold.stats` — placed vs. unplaced contig summary

### Stage 9 — Repeat Masking

Repeat masking is a mandatory prerequisite for structural gene annotation. BRAKER3 and Funannotate require **soft-masked** sequence (lowercase repeat bases) to avoid predicting gene models inside transposable elements.

**RepeatModeler2** constructs a custom repeat library from the assembly using RECON and RepeatScout, then classifies families by homology to known repeats. The `-LTRStruct` flag activates structural LTR retrotransposon detection via LTR_Harvest and LTR_retriever. Using a custom library is critical for *Trichoderma*, which has undergone significant TE expansion not fully represented in generic Dfam databases.

**RepeatMasker** applies the custom library to soft-mask (`-xsmall`) all repeat instances in the scaffolded assembly.

> **Note:** RepeatMasker requires an absolute path for `-dir`. Use `$(pwd)` or `${PROJECT_DIR}` to construct the path, or RepeatMasker will silently write output to a temporary `RM_*/` directory.

Expected repeat content: 5–12% for *T. harzianum*.

**Key outputs:**
- `ragtag.scaffold.fasta.masked` — soft-masked assembly; input for all annotation steps
- `ragtag.scaffold.fasta.tbl` — repeat content summary
- `ragtag.scaffold.fasta.gff` — GFF3 repeat annotation track
- `ragtag.scaffold.fasta.out` — full tab-delimited repeat hit table

---

## Pipeline Summary

```
Raw reads (SRA)
    ↓  SeqFetcher
Quality control (SeqKit, LongReadSum, FastQC, MultiQC)
    ↓  fastplong / fastp
Filtered reads
    ↓  Kraken2
Decontaminated reads
    ↓  Flye / Raven / wtdbg2
Raw assemblies → BUSCO + QUAST → select best (Flye)
    ↓  minimap2 + Racon (×2) → Polypolish → NextPolish (×2)
Polished assembly → BUSCO + QUAST + Merqury
    ↓  RagTag scaffold
Scaffolded assembly
    ↓  RepeatModeler2 → RepeatMasker
Soft-masked assembly ← ready for annotation
```

---

## Dependencies

All tools are managed via [mamba](https://mamba.readthedocs.io/) environments.

| Tool | Environment | Purpose |
|---|---|---|
| SeqFetcher | `ncbi_tools` | SRA download |
| SeqKit | `seqkit` | Read statistics |
| LongReadSum | `longreadsum` | PacBio read-length profiling |
| FastQC | `fastqc` | Illumina per-base QC |
| MultiQC | `multiqc` | Aggregate QC reports |
| fastplong | `fastplong` | PacBio long-read filtering |
| fastp | `fastp` | Illumina adapter trimming + error correction |
| Kraken2 | `kraken2` | Taxonomic decontamination |
| Flye | `flye292` | Long-read assembler (repeat graph) |
| Raven | `raven` | Long-read assembler (string graph) |
| wtdbg2 / wtpoa-cns | `wtdbg2` | Long-read assembler (fuzzy dBG) |
| minimap2 | `minimap2` | Long-read alignment for Racon |
| Racon | `racon1420` | Long-read consensus polishing |
| BWA | `bwa` | Short-read alignment for Polypolish |
| Polypolish | `polypolish` | All-alignment short-read polishing |
| NextPolish | `nextpolish39` | Final SNP/indel correction |
| BUSCO | `busco` | Gene completeness assessment |
| QUAST | `quast` | Assembly contiguity statistics |
| Merqury | `merqury` | Reference-free base accuracy (QV) |
| RagTag | `ragtag` | Reference-guided scaffolding |
| RepeatModeler2 | `repeatmodeler` | De novo repeat library construction |
| RepeatMasker | `repeatmodeler` | Genome soft-masking |

---

## Data Availability

| Accession | Technology | Description | Total Bases | Coverage |
|---|---|---|---|---|
| [SRR10848482](https://www.ncbi.nlm.nih.gov/sra/SRR10848482) | PacBio RS2/Sequel | Long reads | 7.1 Gb | ~178× |
| [SRR10848483](https://www.ncbi.nlm.nih.gov/sra/SRR10848483) | Illumina HiSeq 4000 | Short reads (run 1) | 1.5 Gb | ~38× |
| [SRR10848484](https://www.ncbi.nlm.nih.gov/sra/SRR10848484) | Illumina HiSeq 4000 | Short reads (run 2) | 0.95 Gb | ~23× |

Scaffolding reference: [GCA_019097725.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_019097725.1/) (*T. harzianum* SYAU_Tha_1.0)

BioProject: [PRJNA596042](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA596042)

---

## Recommended Next Steps

The soft-masked assembly is ready for structural annotation:

1. **BRAKER3** — combines RNA-seq evidence and protein homology via Augustus and GeneMark; most accurate for fungi when RNA-seq data is available
2. **Funannotate** — end-to-end fungal annotation pipeline with good defaults for *Trichoderma*

Functional annotation:

- **InterProScan** — domain and protein family annotation (Pfam, PANTHER, TIGRFAM)
- **eggNOG-mapper** — COG/KEGG/GO term assignment via orthology
- **antiSMASH** (fungal mode) — secondary metabolite biosynthetic gene cluster identification, highly relevant for *T. harzianum* biocontrol research
- **dbCAN / MEROPS** — CAZyme and protease annotation

---

## Author

**Gianlucca de Urzêda Alves**
Version 3.0 — June 2026
