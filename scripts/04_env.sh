#!/usr/bin/env bash

# Get environment names and export them to envs/ and specs/

mkdir -p envs specs

for env in \
    ncbi_tools \
    seqkit \
    longreadsum \
    fastqc \
    multiqc \
    fastplong \
    fastp \
    kraken2 \
    flye292 \
    raven \
    wtdbg2 \
    minimap2 \
    racon1420 \
    bwa \
    polypolish \
    nextpolish39 \
    busco \
    quast \
    merqury \
    ragtag \
    repeatmodeler
do
    conda env export -n "$env" > "envs/${env}.yml"

    conda list \
        -n "$env" \
        --explicit \
        > "specs/${env}.txt"
done