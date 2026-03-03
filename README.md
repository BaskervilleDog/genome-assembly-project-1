# 🧬 Hybrid Genome Assembly Pipeline — *Trichoderma harzianum* TW11

A reproducible, multi-assembler pipeline for *de novo* genome assembly using PacBio long reads and Illumina short-read polishing, developed for the fungal *Trichoderma harzianum* TW11.

---

## Overview

This project implements a complete hybrid genome assembly workflow that:

1. **Downloads** raw sequencing data from SRA (NCBI/ENA)
2. **Quality controls** both PacBio and Illumina reads
3. **Filters** long reads with FILTlong
4. **Assembles** the genome with three independent assemblers (Flye, Raven, wtdbg2)
5. **Polishes** each assembly with Illumina reads via Pilon (2 rounds)
6. **Scaffolds** polished assemblies against a reference genome with RagTag
7. **Evaluates** assembly quality with QUAST and BUSCO

Using multiple assemblers mitigates algorithm-specific biases and allows selection of the best final assembly based on objective metrics.

---

## Biological Context

| Feature | Detail |
|---|---|
| **Organism** | *Trichoderma harzianum* TW11 |
| **Application** | Fungal biocontrol agent genomics |
| **Long reads** | PacBio RS2/Sequel (SRR10848482) — contiguity across repeats |
| **Short reads** | Illumina PE (SRR10848483, SRR10848484) — base-level accuracy |
| **Estimated genome size** | ~40 Mb |
| **BUSCO lineage** | `hypocreales_odb10` |

---

## Repository Structure

```
.
├── 01_download.sh              # SRA data download via seqfetcher
├── 02_assembly_pipeline.sh     # Main hybrid assembly pipeline
├── 03_read_mapping.sh          # Post-assembly read mapping & QC
├── seqfetcher/                 # Submodule: SRA download tool
└── assembly_output/
    └── T_harzianum_TW11/
        ├── qc/                 # FastQC, MultiQC, seqkit reports
        ├── filtered_reads/     # FILTlong-filtered PacBio reads
        ├── assemblies/         # Raw assemblies (flye/, raven/, wtdbg2/)
        ├── polished_assemblies/# Pilon-polished FASTAs
        ├── scaffolded_assemblies/ # RagTag-scaffolded FASTAs
        ├── assembly_metrics/   # QUAST reports, BUSCO results
        ├── logs/               # Per-tool log files
        └── pipeline_metadata.txt
```

---

## Pipeline Stages

### Stage 1 — Data Acquisition (`01_download.sh`)

- Clones and installs [seqfetcher](https://github.com/BaskervilleDog/seqfetcher)
- Downloads three SRA accessions via ENA
- Outputs raw FASTQ files to `seqfetcher/downloads/fastq/`

### Stage 2 — Assembly Pipeline (`02_assembly_pipeline.sh`)

```
Raw PacBio reads
      │
      ▼
  [seqkit QC]
      │
      ▼
  [FILTlong] ──── min length: 1000 bp
      │
      ▼
  [seqkit QC on filtered reads]
      │
      ├──────────────────────────────┐
      │                              │
  [Flye]                         [Raven]    [wtdbg2]
  --pacbio-raw                   OLC         fuzzy dBG
      │                              │           │
      └──────────────┬───────────────┘           │
                     │◄──────────────────────────┘
                     ▼
              [Pilon Polishing]
              (2 rounds × 3 assemblies)
              BWA-MEM + samtools + Pilon
                     │
                     ▼
              [RagTag Scaffolding]
              (reference-guided)
                     │
                     ▼
           [QUAST + BUSCO Assessment]
```

**Assembler rationale:**

| Assembler | Algorithm | Strength |
|---|---|---|
| Flye | Repeat graph | Handles tandem repeats and heterozygosity |
| Raven | OLC | Fast; often more contiguous assemblies |
| wtdbg2 | Fuzzy de Bruijn graph | Memory-efficient; good for large genomes |

**Polishing rationale:** PacBio RS2 reads have ~13–15% indel error rate. Two rounds of Pilon with high-accuracy Illumina reads reduce the error rate to <0.01%, suitable for gene annotation and comparative genomics.

### Stage 3 — Read Mapping QC (`03_read_mapping.sh`)

- Maps PacBio reads back to the selected assembly with `bwa mem -x pacbio`
- Converts, sorts, and indexes the BAM file
- Runs `bamtools stats` and `qualimap bamqc` to validate assembly completeness

**Result (Flye assembly):** 97.3% of PacBio reads mapped back successfully.

---

## Dependencies

All tools are managed via [mamba](https://mamba.readthedocs.io/) environments.

| Tool | Environment | Purpose |
|---|---|---|
| seqfetcher | `ncbi_tools` | SRA download |
| seqkit | `seqkit` | Read statistics |
| FILTlong | `filtlong` | Long-read filtering |
| FastQC / MultiQC | `fastqc` / `multiqc` | Illumina QC |
| fastp | `fastp` | Illumina adapter trimming |
| Flye | `flye` | Long-read assembler |
| Raven | `raven` | Long-read assembler |
| wtdbg2 / wtpoa-cns | `wtdbg2` | Long-read assembler |
| BWA | `bwa` | Short-read mapping |
| samtools | `samtools` | BAM processing |
| Pilon | `pilon` | Illumina polishing |
| RagTag | `ragtag` | Reference-guided scaffolding |
| QUAST | `quast` | Assembly statistics |
| BUSCO | `busco` | Gene completeness |
| bamtools | `bamtools` | BAM statistics |
| qualimap | `qualimap` | Alignment QC |

---

## Usage

### 1. Download data

```bash
bash 01_download.sh
```

### 2. Run assembly pipeline

```bash
bash 02_assembly_pipeline.sh \
    --pacbio seqfetcher/downloads/fastq/SRR10848482_subreads.fastq.gz \
    --illumina-r1 seqfetcher/downloads/fastq/SRR10848483_1.fastq.gz,seqfetcher/downloads/fastq/SRR10848484_1.fastq.gz \
    --illumina-r2 seqfetcher/downloads/fastq/SRR10848483_2.fastq.gz,seqfetcher/downloads/fastq/SRR10848484_2.fastq.gz \
    --genome-size 40m \
    --sample-id T_harzianum_TW11 \
    --busco-lineage hypocreales_odb10 \
    --scaffold-ref downloads/GCF_003025095.1/GCF_003025095.1_Triha_v1.0_genomic.fna \
    --enable-scaffolding \
    --threads 12
```

### 3. Validate assembly

```bash
bash 03_read_mapping.sh
```

### Key options

| Flag | Default | Description |
|---|---|---|
| `--genome-size` | required | e.g. `40m`, `4.5m`, `3.2g` |
| `--pilon-rounds` | 2 | Number of polishing iterations |
| `--skip-flye` / `--skip-raven` / `--skip-wtdbg2` | — | Disable individual assemblers |
| `--busco-lineage` | `bacteria_odb10` | BUSCO dataset name |
| `--target-bases` | 0 (off) | Downsample long reads (e.g. `6000000000`) |
| `--no-resume` | — | Force re-run all steps |
| `--cleanup-all` | — | Delete SAMs and BAMs after use |

---

## Pipeline Features

- **Checkpointing** (`--resume`): skips completed steps based on output file existence, enabling safe re-runs after interruption
- **Multi-library support**: accepts any number of Illumina PE library pairs as comma-separated lists
- **Coverage warnings**: alerts if post-filter coverage drops below 30× or read N50 is unexpectedly low
- **Pilon convergence detection**: reports when the number of corrections is very low, indicating assembly quality has stabilized
- **Structured logging**: timestamped logs per tool per step, stored in `logs/`
- **Metadata tracking**: all parameters, input files, and tool versions written to `pipeline_metadata.txt`

---

## Data Availability

| Accession | Type | Library |
|---|---|---|
| [SRR10848482](https://www.ncbi.nlm.nih.gov/sra/SRR10848482) | PacBio subreads | Long reads |
| [SRR10848483](https://www.ncbi.nlm.nih.gov/sra/SRR10848483) | Illumina PE | Short reads (lib 1) |
| [SRR10848484](https://www.ncbi.nlm.nih.gov/sra/SRR10848484) | Illumina PE | Short reads (lib 2) |

Scaffolding reference: [GCF_003025095.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003025095.1/) (*T. harzianum* Triha v1.0)

---

## Author

**Gianlucca de Urzêda Alves**
Version 2.0 — February 2026