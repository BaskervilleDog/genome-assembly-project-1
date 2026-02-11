#!/bin/bash

# Data: 28/01/2026
# Autor: Gianlucca de Urzêda Alves

# Script construido para obter sequenciamentos de genoma
# Ferramentas usadas: mamba | jellyfish

### Estimar tamanho do genoma usando jellyfish
# Decompress the files
gunzip -c ./filtered_data/SRR10848483_1.fastq.gz > ./filtered_data/SRR10848483_1.fastq
gunzip -c ./filtered_data/SRR10848483_2.fastq.gz > ./filtered_data/SRR10848483_2.fastq

# Run jellyfish on uncompressed files
mamba run -n jellyfish jellyfish count -C -m 21 -s 5G -t 10 -o SRR10848483_reads.jf \
    ./filtered_data/SRR10848483_1.fastq \
    ./filtered_data/SRR10848483_2.fastq

# Clean up after
rm ./filtered_data/SRR10848483_1.fastq ./filtered_data/SRR10848483_2.fastq

mamba run -n jellyfish jellyfish histo -t 10 SRR10848483_reads.jf > SRR10848483_reads.histo

awk '$1 >= 10 {sum += $1 * $2} END {print sum}' SRR10848483_reads.histo

# **Manual calculation:**
# Genome size = Total k-mers / Coverage peak
# Total k-mers = 1204126364
# Coverage peak = 31
# Genome size = 38,842,785.93 (38,8 Mb)

# Decompress the files
gunzip -c ./filtered_data/SRR10848484_1.fastq.gz > ./filtered_data/SRR10848484_1.fastq
gunzip -c ./filtered_data/SRR10848484_2.fastq.gz > ./filtered_data/SRR10848484_2.fastq

# Run jellyfish on uncompressed files
mamba run -n jellyfish jellyfish count -C -m 21 -s 5G -t 10 -o SRR10848484_reads.jf \
    ./filtered_data/SRR10848484_1.fastq \
    ./filtered_data/SRR10848484_2.fastq

# Clean up after
rm ./filtered_data/SRR10848484_1.fastq ./filtered_data/SRR10848484_2.fastq

mamba run -n jellyfish jellyfish histo -t 10 SRR10848484_reads.jf > SRR10848484_reads.histo

awk '$1 >= 10 {sum += $1 * $2} END {print sum}' SRR10848484_reads.histo

# **Manual calculation:**
# Genome size = Total k-mers / Coverage peak
# Total k-mers = 779147778
# Coverage peak = 28
# Genome size = 27,826,706.35 (27.82 Mb)

## Calculo de tamanho com GenomeScope

# Biblioteca Illumina IS1000
mamba run -n genomescope genomescope2 -i SRR10848483_reads.histo \
    -o SRR10848483_genomescope_output \
    -k 21 \
    -p 1 \
    -l 150 \
    --fitted_hist

# Genome Haploid Length         36,522,162 bp --  40,555,625 bp  -- 45,590,594 bp

# Biblioteca Illumina IS270
mamba run -n genomescope genomescope2 -i SRR10848484_reads.histo \
    -o SRR10848484_genomescope_output \
    -k 21 \
    -p 1 \
    -l 150 \
    --fitted_hist

# Genome Haploid Length         4,017,228 bp --  33,616,007 bp -- Inf bp

### Assembly usando Flye para long reads

mamba run -n flye flye --pacbio-raw ./filtered_data/SRR10848482_subreads.fastq.gz \
     --out-dir flye_assembly \
     --genome-size 40m \
     --threads 12

#[2026-01-29 15:49:32] INFO: Assembly statistics: 13:50 - 15:50 
#
#        Total length:   41120800
#        Fragments:      25
#        Fragments N50:  3472118
#        Largest frg:    7660367
#        Scaffolds:      0
#        Mean coverage:  143

### Polindo assembly gerado pelo Flye

### Indexando o assembly usando bwa
mamba run -n bwa bwa index ./flye_assembly/assembly.fasta

### Mapeando as bibliotecas Illumina no assembly

### Biblioteca 1
mamba run -n bwa bwa mem -t 12 ./flye_assembly/assembly.fasta \
    ./filtered_data/SRR10848483_1.fastq.gz ./filtered_data/SRR10848483_2.fastq.gz | \
    mamba run -n samtools samtools view -@ 12 -bS - | \
    mamba run -n samtools samtools sort -@ 12 -o ./flye_assembly/mapped_lib1.bam -

### Biblioteca 2
mamba run -n bwa bwa mem -t 12 ./flye_assembly/assembly.fasta \
    ./filtered_data/SRR10848484_1.fastq.gz ./filtered_data/SRR10848484_2.fastq.gz | \
    mamba run -n samtools samtools view -@ 12 -bS - | \
    mamba run -n samtools samtools sort -@ 12 -o ./flye_assembly/mapped_lib2.bam -

### Unir os arquivos de alinhamento
mamba run -n samtools samtools merge -@ 12 \
    ./flye_assembly/mapped_reads_merged.bam \
    ./flye_assembly/mapped_lib1.bam \
    ./flye_assembly/mapped_lib2.bam

### Indexar o arquivo bam unido
mamba run -n samtools samtools index ./flye_assembly/mapped_reads_merged.bam

### Executar o Pilon para polir o assembly: 16:00 17:10
java -Xms4G -Xmx32G -jar  /home/gianlucca/miniforge3/envs/pilon/share/pilon-1.24-0/pilon.jar \
    --genome ./flye_assembly/assembly.fasta \
    --frags ./flye_assembly/mapped_reads_merged.bam \
    --output pilon_round1 \
    --outdir ./flye_assembly/pilon_round1 \
    --changes

### Segundo round de polimento

### Indexar o assembly polido no round 1
mamba run -n bwa bwa index ./flye_assembly/pilon_round1/pilon_round1.fasta

### Mapeando as bibliotecas Illumina no assembly polido no round 1

### Biblioteca 1
mamba run -n bwa bwa mem -t 12 ./flye_assembly/pilon_round1/pilon_round1.fasta \
    ./filtered_data/SRR10848483_1.fastq.gz ./filtered_data/SRR10848483_2.fastq.gz | \
    mamba run -n samtools samtools view -@ 12 -bS - | \
    mamba run -n samtools samtools sort -@ 12 -o ./flye_assembly/mapped_lib1_round2.bam -

### Biblioteca 2
mamba run -n bwa bwa mem -t 12 ./flye_assembly/pilon_round1/pilon_round1.fasta \
    ./filtered_data/SRR10848484_1.fastq.gz ./filtered_data/SRR10848484_2.fastq.gz | \
    mamba run -n samtools samtools view -@ 12 -bS - | \
    mamba run -n samtools samtools sort -@ 12 -o ./flye_assembly/mapped_lib2_round2.bam -

### Unir os arquivos de alinhamento
mamba run -n samtools samtools merge -@ 12 \
    ./flye_assembly/mapped_reads_merged_round2.bam \
    ./flye_assembly/mapped_lib1_round2.bam \
    ./flye_assembly/mapped_lib2_round2.bam

### Indexar o arquivo bam unido
mamba run -n samtools samtools index ./flye_assembly/mapped_reads_merged_round2.bam

### Executar o Pilon para polir o assembly
java -Xms4G -Xmx32G -jar  /home/gianlucca/miniforge3/envs/pilon/share/pilon-1.24-0/pilon.jar \
    --genome ./flye_assembly/pilon_round1/pilon_round1.fasta \
    --frags ./flye_assembly/mapped_reads_merged_round2.bam \
    --output pilon_round2 \
    --outdir ./flye_assembly/pilon_round2 \
    --changes

# Original assembly stats
mamba run -n assembly_stats assembly-stats ./flye_assembly/assembly.fasta > ./flye_assembly/original_assembly_stats.txt
cat ./flye_assembly/original_assembly_stats.txt

# Round 1 stats
mamba run -n assembly_stats assembly-stats ./flye_assembly/pilon_round1/pilon_round1.fasta > ./flye_assembly/round1_assembly_stats.txt
cat ./flye_assembly/round1_assembly_stats.txt

# Final assembly stats
mamba run -n assembly_stats assembly-stats ./flye_assembly/pilon_round2/pilon_round2.fasta > ./flye_assembly/final_assembly_stats.txt
cat ./flye_assembly/final_assembly_stats.txt

# Analisar a qualidade dos assemblies com QUAST

mamba run -n quast quast ./flye_assembly/assembly.fasta -o ./flye_assembly/original_assembly_quast

mamba run -n quast quast ./flye_assembly/pilon_round1/pilon_round1.fasta -o ./flye_assembly/pilon_round1_assembly_quast

mamba run -n quast quast ./flye_assembly/pilon_round2/pilon_round2.fasta -o ./flye_assembly/pilon_round2_assembly_quast

### Analisar os assemblies de acordo com o BUSCO

### Encontrar o dataset de linhagem

mamba run -n busco busco --list-datasets

### Baixar a linhagem que corresponde ao organismo de interesse
mamba run -n busco busco --download hypocreaceae_odb12

### Executar o BUSCO para cada assembly

mamba run -n busco busco \
    -i ./flye_assembly/assembly.fasta \
    -o ./flye_assembly/busco_results_original \
    -m genome \
    -l hypocreaceae_odb12 \
    -c 12 \
    --offline

mamba run -n busco busco \
    -i ./flye_assembly/pilon_round1/pilon_round1.fasta \
    -o ./flye_assembly/busco_results_polished1 \
    -m genome \
    -l hypocreaceae_odb12 \
    -c 12 \
    --offline

mamba run -n busco busco \
    -i ./flye_assembly/pilon_round2/pilon_round2.fasta \
    -o ./flye_assembly/busco_results_polished2 \
    -m genome \
    -l hypocreaceae_odb12 \
    -c 12 \
    --offline

### Realizar o Scaffolding
### Pode ser realizado usando Referencia ou reads

### Referencia
# Genoma de referencia para Trichoderma harzianum GCF_003025095.1 (Trichoderma harzianum CBS 226) versão Triha v1.0
mamba run -n ragtag /home/gianlucca/miniforge3/envs/ragtag/bin/ragtag.py scaffold ./downloads/GCF_003025095.1/ncbi_dataset/data/GCF_003025095.1/GCF_003025095.1_Triha_v1.0_genomic.fna ./flye_assembly/pilon_round2/pilon_round2.fasta -o ./flye_assembly/ragtag_output -t 12

### Gap Filling
mamba run -n ragtag /home/gianlucca/miniforge3/envs/ragtag/bin/ragtag.py patch ./flye_assembly/ragtag_output/ragtag.scaffold.fasta ./flye_assembly/pilon_round2/pilon_round2.fasta -o ./flye_assembly/ragtag_output_patched -t 12

### Avaliar com QUAST
mamba run -n quast quast ./flye_assembly/ragtag_output/ragtag.scaffold.fasta -o ./flye_assembly/assembly_scaffold_quast

### Avaliar com BUSCO
mamba run -n busco busco \
    -i ./flye_assembly/ragtag_output/ragtag.scaffold.fasta \
    -o ./flye_assembly/busco_results_scaffold \
    -m genome \
    -l hypocreaceae_odb12 \
    -c 12 \
    --offline


### Realizar assembly usando Raven 11:11 13:40
mkdir ./raven_assembly

mamba run -n raven raven --threads 12 \
      --polishing-rounds 3 \
      -k 15 \
      -w 5 \
      ./filtered_data/SRR10848482_subreads.fastq.gz \
      > ./raven_assembly/raven_assembly.fasta \
      2> ./raven_assembly/raven.log

### Avaliar com QUAST
mamba run -n quast quast ./raven_assembly/raven_assembly.fasta -o ./raven_assembly/assembly_quast

### Avaliar com BUSCO
mamba run -n busco busco \
    -i ./raven_assembly/raven_assembly.fasta \
    -o ./raven_assembly/busco_results \
    -m genome \
    -l hypocreaceae_odb12 \
    -c 12 \
    --offline

### Realizar assembly usando wtdbg2 13:49 13:55
mkdir ./wtdbg2_assembly

mamba run -n wtdbg2 wtdbg2 -x rs \
       -g 40m \
       -i ./filtered_data/SRR10848482_subreads.fastq.gz \
       -t 12 \
       -L 5000 \
       -S 1 \
       -fo ./wtdbg2_assembly/assembly_wtdbg2

# 13:57 14:11 
mamba run -n wtdbg2 wtpoa-cns -t 12 \
          -i ./wtdbg2_assembly/assembly_wtdbg2.ctg.lay.gz \
          -fo ./wtdbg2_assembly/wtdbg2_assembly.fasta

### Avaliar com QUAST
mamba run -n quast quast ./wtdbg2_assembly/wtdbg2_assembly.fasta -o ./wtdbg2_assembly/assembly_quast

### Avaliar com BUSCO
mamba run -n busco busco \
    -i ./wtdbg2_assembly/wtdbg2_assembly.fasta \
    -o ./wtdbg2_assembly/busco_results \
    -m genome \
    -l hypocreaceae_odb12 \
    -c 12 \
    --offline

# Avaliar com Merqury
mkdir ./merqury -p
mkdir -p logs/flye_assembly
mkdir -p logs/raven_assembly

mamba run -n merqury meryl k=21 count \
    ./filtered_data/SRR10848483_1.fastq.gz ./filtered_data/SRR10848483_2.fastq.gz \
    ./filtered_data/SRR10848484_1.fastq.gz ./filtered_data/SRR10848484_2.fastq.gz \
    output ./merqury/illumina.meryl

mamba run -n merqury \
  /home/gianlucca/miniforge3/envs/merqury/bin/merqury.sh \
  ./merqury/illumina.meryl \
  ./flye_assembly/ragtag_output/ragtag.scaffold.fasta \
  ./flye_assembly/merqury

mamba run -n merqury \
  /home/gianlucca/miniforge3/envs/merqury/bin/merqury.sh \
  ./merqury/illumina.meryl \
  ./raven_assembly/raven_assembly.fasta \
  ./raven_assembly/merqury