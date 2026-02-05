#!/bin/bash
#==============================================================================
# ASSEMBLY ANNOTATION USING Braker3
#==============================================================================
# Author: Gianlucca de Urzêda Alves
# Date: 03/02/2026
# Description:
#   - Stablish a pipeline for ab-initio and evidence based annotations
#==============================================================================


#==============================================================================
# CONFIGURATION
#==============================================================================

# Install braker (if needed)
mamba create -n braker -c bioconda -c conda-forge braker3 perl-hash-merge perl-mce

mkdir braker_annotation

# Install GeneMark ES 
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_zyAHJ/gmes_linux_64_4.tar.gz
# Get software key
wget http://genemark.bme.gatech.edu/tmp/GMtool_iK0YV/gm_key.gz

tar -xvzf gmes_linux_64_4.tar.gz

gunzip gm_key.gz

mv /home/gianlucca/genome-assembly-project-1/gm_key ~/.gm_key
chmod 600 ~/.gm_key

chmod +x /home/gianlucca/genome-assembly-project-1/gmes_linux_64_4/*

export GENEMARK_PATH=/home/gianlucca/genome-assembly-project-1/gmes_linux_64_4

# Install Augustus
git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augustus
make clean
make

mamba activate braker

export AUGUSTUS_BIN_PATH=/home/gianlucca/genome-assembly-project-1/Augustus/bin
export AUGUSTUS_SCRIPTS_PATH=/home/gianlucca/genome-assembly-project-1/Augustus/scripts
export AUGUSTUS_CONFIG_PATH=/home/gianlucca/genome-assembly-project-1/Augustus/config
export PERL5LIB=/home/gianlucca/miniforge3/envs/braker/lib/perl5/site_perl:/home/gianlucca/miniforge3/envs/braker/lib/perl5/vendor_perl

braker.pl \
  --genome=flye_assembly/ragtag_output/ragtag.scaffold.fasta \
  --species=Tharzianum_T11W \
  --esmode \
  --fungus \
  --softmasking \
  --threads=12 \
  --gff3 \
  --workingdir=./braker_annotation

grep -c "gene" ./braker_annotation/braker.gff3

mamba create -n agat -c bioconda agat

mamba run -n agat agat_sp_statistics.pl --gff ./braker_annotation/braker.gff3 -o ./braker_annotation/annotation_stats.txt

mamba run -n busco busco \
        -i "./braker_annotation/braker.aa" \
        -o "./braker_annotation/busco" \
        -m protein \
        -l hypocreaceae_odb12 \
        -c "12" \
        --offline -f \
        2>&1 | tee "./braker_annotation/busco/busco.log"

mamba run -n seqkit seqkit seq -w 0 ./braker_annotation/braker.aa > ./braker_annotation/braker.oneline.faa

mamba run -n seqkit seqkit fx2tab -n -l ./braker_annotation/braker.oneline.faa \
| awk '
  {
    id=$1; len=$2;
    gene=id;
    sub(/\.t[0-9]+$/, "", gene);

    if (!(gene in best) || len > best[gene]) {
      best[gene]=len;
      bestid[gene]=id;
    }
  }
  END {
    for (g in bestid) print bestid[g];
  }
' > ./braker_annotation/best_ids.txt

mamba run -n seqkit seqkit grep -f ./braker_annotation/best_ids.txt ./braker_annotation/braker.oneline.faa > ./braker_annotation/braker.longest.faa

grep -c "^>" ./braker_annotation/braker.aa
grep -c "^>" ./braker_annotation/braker.longest.faa

#mamba create -n eggnog-mapper -c bioconda eggnog-mapper

#mkdir ./eggnog_mapper_annotation

#mamba run -n eggnog_mapper download_eggnog_data.py -y --data_dir ./eggnog_mapper_annotation

#mamba run -n eggnog_mapper emapper.py -i ./braker_annotation/braker.longest.faa \
#  --itype proteins \
#  -m mmseqs \
#  --data_dir ./eggnog_mapper_annotation \
#  -o trichoderma_annotation \
#  --output_dir ./eggnog_mapper_annotation \
#  --cpu 12

mamba create -n diamond -c bioconda diamond

mkdir -p ./diamond_annotation

cd diamond_annotation
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
cd ..

mamba run -n diamond diamond makedb --in ./diamond_annotation/uniprot_sprot.fasta.gz -d swissprot

mamba run -n diamond diamond blastp -d swissprot \
  -q ./braker_annotation/braker.longest.faa \
  -o ./diamond_annotation/trichoderma_swissprot_annotations.txt \
  --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
  --max-target-seqs 1 \
  --evalue 1e-5 \
  --threads 12 \
  --sensitive


mamba create -n interproscan -c bioconda interproscan

mkdir -p interproscan_annotation

sed 's/\*//g' ./braker_annotation/braker.longest.faa > ./braker_annotation/braker.longest.clean.faa
grep -c "\*" braker_annotation/braker.longest.clean.faa

mamba run -n interproscan hmmpress ~/miniforge3/envs/interproscan/share/InterProScan/data/superfamily/1.75/hmmlib_1.75

mamba run -n interproscan interproscan.sh \
  -i braker_annotation/braker.longest.clean.faa \
  -appl TIGRFAM,FunFam,SFLD,SUPERFAMILY,PANTHER,Gene3D,Hamap,ProSiteProfiles,Coils,SMART,PRINTS,PIRSR,AntiFam,Pfam,MobiDBLite,PIRSF \
  -f TSV,GFF3 \
  -goterms \
  -iprlookup \
  -pa \
  -cpu 12 \
  -d interproscan_annotation/


mamba create -n orthofinder -c bioconda orthofinder

mkdir ./comparative_genomics/protein_sequences -p

cp ./braker_annotation/braker.longest.faa ./comparative_genomics/protein_sequences/braker.faa

mamba run -n ncbi_tools seqfetcher search --organism "Trichoderma"

mamba run -n ncbi_tools seqfetcher download --proteome --assembly GCF_020647795.1 
mamba run -n ncbi_tools seqfetcher download --proteome --assembly GCF_020647865.1 
mamba run -n ncbi_tools seqfetcher download --proteome --assembly GCF_028502605.1 
mamba run -n ncbi_tools seqfetcher download --proteome --assembly GCF_000167675.1 
mamba run -n ncbi_tools seqfetcher download --proteome --assembly GCF_000170995.1

cp ./downloads/GCF_003025095.1/ncbi_dataset/data/GCF_003025095.1/protein.faa \
 ./comparative_genomics/protein_sequences/Trichoderma_harzianum.faa

cp ./downloads/proteomes/GCF_000167675.1/ncbi_dataset/data/GCF_000167675.1/protein.faa \
 ./comparative_genomics/protein_sequences/Trichoderma_reesei.faa

cp ./downloads/proteomes/GCF_000170995.1/ncbi_dataset/data/GCF_000170995.1/protein.faa \
 ./comparative_genomics/protein_sequences/Trichoderma_virens.faa

cp ./downloads/proteomes/GCF_020647795.1/ncbi_dataset/data/GCF_020647795.1/protein.faa \
 ./comparative_genomics/protein_sequences/Trichoderma_atroviride.faa

cp ./downloads/proteomes/GCF_020647865.1/ncbi_dataset/data/GCF_020647865.1/protein.faa \
 ./comparative_genomics/protein_sequences/Trichoderma_asperellum.faa

cp ./downloads/proteomes/GCF_028502605.1/ncbi_dataset/data/GCF_028502605.1/protein.faa \
 ./comparative_genomics/protein_sequences/Trichoderma_breve.faa

mamba run -n orthofinder orthofinder \
  -f ./comparative_genomics/protein_sequences/ \
  -t 12 \
  -a 12 \
  -M msa \
  -A mafft \
  -S diamond \
  -o ./comparative_genomics/orthofinder_results 