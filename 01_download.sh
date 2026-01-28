#!/bin/bash

# Data: 28/01/2026
# Autor: Gianlucca de Urzêda Alves

# Script construido para obter sequenciamentos de genoma do SRA e realizar controle de qualidade
# Ferramentas usadas: mamba | sra-tools | FASTQC | MultiQC | fastp | nanoplot | longQC | filtlong

### Ativar ambiente contendo as ferramentas necessárias para download dos dados
#### sra-tools | parallel-fastq-dump | ncbi-datasets-cli | ncbi-entrez-direct

# 1. Download de ferramenta para downloads ###############################################

git clone https://github.com/BaskervilleDog/seqfetcher.git
cd seqfetcher

chmod +x install.sh
./install.sh

## Baixar dados do banco SRA

echo "SRR10848482
SRR10848483
SRR10848484" > accession_list.txt

mamba run -n ncbi_tools seqfetcher download --sra-method ena --sra-accession-file accession_list.txt

cd ..

# 2. Controle de Qualidade ###############################################

## Obter estatisticas iniciais usando seqkit

mamba run -n seqkit seqkit stats ./seqfetcher/downloads/fastq/*.fastq.gz

## Controle de qualidade usando FASTQC

### Criar diretório para conter os arquivos do relatório de qualidade
mkdir ./qc_report -p

### Executar programa fastqc
mamba run -n fastqc fastqc ./seqfetcher/downloads/fastq/*.fastq.gz -o ./qc_report --threads 4

### Executar programa multiqc na pasta contendo todos relatorios fastqc
mamba run -n multiqc multiqc ./qc_report 

### QC de long reads realizado por diferentes softwares
#### Nanoplot
mamba run -n nanoplot NanoPlot --fastq ./seqfetcher/downloads/fastq/*subreads.fastq.gz -o ./qc_report/nanoplot_output

#### LongQC
for f in ./seqfetcher/downloads/fastq/*subreads.fastq.gz; do
  sample=$(basename "$f" .subreads.fastq.gz)
  sudo docker run --rm \
    -v "$(pwd)/seqfetcher/downloads/fastq:/input" \
    -v "$(pwd)/qc_report/longqc_output:/output" \
    longqc sampleqc \
      -x pb-rs2 \
      -p 4 \
      -o "/output/$sample" \
      "/input/$(basename "$f")"
done

### Relatório dos dados brutos foram gerados
### Filtragem das reads dos dados brutos

### Criar diretório para conter os dados filtrados
mkdir -p filtered_data

### Executar programa fastp para cada arquivo baixado
# Loop for onde vai iterar sobre todos arquivos na pasta DATA que terminam com "_1.fastq.gz";
# Cada nome de arquivo forward é armazenado na variável r1
# Cria nomes de arquivos substituindo "_1.fastq.gz" com "_2.fastq.gz"
# sample extrai o nome de cada amostra (basename remove o caminho de diretorio do nome;  "_1.fastq.gz" é retirado do nome)
## DATA/sampleA_1.fastq.gz → sampleA
# Executa fastp 
## -i são inputs forward
## -I são inputs reverse
## -o são outputs forward
## -O são outputs reverse
## -q remove todas bases com qualidade Phred < 20 (trimming)
## -l 100 descarta reads com menos de 100 pares de base após o trimming
## -D desativa avaliação de duplicação

for r1 in ./seqfetcher/downloads/fastq/*_1.fastq.gz; do
    r2=${r1/_1.fastq.gz/_2.fastq.gz}
    sample=$(basename "$r1" _1.fastq.gz)

    mamba run -n fastp fastp \
        -i "$r1" \
        -I "$r2" \
        -o "filtered_data/${sample}_1.fastq.gz" \
        -O "filtered_data/${sample}_2.fastq.gz" \
        -q 20 \
        -l 100 \
        -D
done

for f in ./seqfetcher/downloads/fastq/*subreads.fastq.gz; do
  sample=$(basename "$f" .fastq.gz)

  mamba run -n filtlong filtlong \
    --min_length 1000 \
    --keep_percent 90 \
    "$f" | gzip > "./filtered_data/${sample}.fastq.gz"
done

# 3. Reavaliar usando fastqc e multiqc ###############################################

### Criar diretório para conter os arquivos do relatório de qualidade
mkdir ./filtered_data/qc_report

### Executar programa fastqc
mamba run -n fastqc fastqc ./filtered_data/*fastq.gz -o ./filtered_data/qc_report --threads 2

### Executar programa multiqc na pasta contendo todos relatorios fastqc
mamba run -n multiqc multiqc ./filtered_data/qc_report/

### QC de long reads realizado por diferentes softwares
#### Nanoplot
mamba run -n nanoplot NanoPlot --fastq ./filtered_data/*subreads.fastq.gz -o ./filtered_data/qc_report/nanoplot_output

#### LongQC
for f in filtered_data/*subreads.fastq.gz; do
  sample=$(basename "$f" .subreads.fastq.gz)
  sudo docker run --rm \
    -v "$(pwd)/seqfetcher/downloads/fastq:/input" \
    -v "$(pwd)/filtered_data/qc_report/longqc_output:/output" \
    longqc sampleqc \
      -x pb-rs2 \
      -p 4 \
      -o "/output/$sample" \
      "/input/$(basename "$f")"
done

###############################################