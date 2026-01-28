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
mamba run -n jellyfish jellyfish count -C -m 21 -s 1G -t 4 -o reads.jf \
    ./filtered_data/SRR10848483_1.fastq \
    ./filtered_data/SRR10848483_2.fastq

# Clean up after
rm ./filtered_data/SRR10848483_1.fastq ./filtered_data/SRR10848483_2.fastq

mamba run -n jellyfish jellyfish histo -t 4 reads.jf > reads.histo

awk '$1 >= 10 {sum += $1 * $2} END {print sum}' reads.histo

# **Manual calculation:**
# Genome size = Total k-mers / Coverage peak
# Total k-mers = 1204126364
# Coverage peak = 31
# Genome size = 38,842,785.93 (38,8 Mb)

## Calculo de tamanho com GenomeScope

mamba run -n genomescope genomescope2 -i reads.histo \
    -o genomescope_output \
    -k 21 \
    -p 1 \
    -l 150 \
    --fitted_hist

# Genome Haploid Length         36,467,499 bp     45,667,458 bp

### Assembly usando Flye para long reads

mamba run -n flye flye --pacbio-raw ./filtered_data/SRR10848482_subreads.fastq.gz \
     --out-dir flye_assembly \
     --genome-size 40m \
     --threads 4

### Ativar ambiente contendo as ferramentas necessárias para download dos dados
#### sra-tools | parallel-fastq-dump | ncbi-datasets-cli | ncbi-entrez-direct