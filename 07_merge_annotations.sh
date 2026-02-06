#!/bin/bash

# Create master annotation table using bash/awk

echo "Creating master annotation table..."

mkdir ./master_annotation

# Get all unique gene IDs
echo "Extracting gene IDs..."
cat braker_annotation/braker.longest.clean.faa | grep ">" | sed 's/>//' | sed 's/ .*//' | sort -u > ./master_annotation/all_genes.txt

# Create header
echo -e "gene_id\tswissprot_hit\tswissprot_desc\tswissprot_evalue\tdomains\tinterpro_ids\tgo_terms\tcazyme_family" > ./master_annotation/master_annotations.tsv

# Process each gene
total=$(wc -l < ./master_annotation/all_genes.txt)
counter=0

while read gene_id; do
    counter=$((counter + 1))
    if [ $((counter % 1000)) -eq 0 ]; then
        echo "Processing gene $counter/$total..."
    fi
    
    # DIAMOND/SwissProt annotation
    swissprot_line=$(grep "^$gene_id\s" diamond_annotation/trichoderma_swissprot_annotations.txt | head -1)
    if [ -n "$swissprot_line" ]; then
        swissprot_hit=$(echo "$swissprot_line" | cut -f2)
        swissprot_desc=$(echo "$swissprot_line" | cut -f7-)
        swissprot_evalue=$(echo "$swissprot_line" | cut -f5)
    else
        swissprot_hit=""
        swissprot_desc=""
        swissprot_evalue=""
    fi
    
    # InterProScan annotations
    interpro_lines=$(grep "^$gene_id\s" interproscan_annotation/braker.longest.clean.faa.tsv | grep -v "MobiDBLite")
    if [ -n "$interpro_lines" ]; then
        domains=$(echo "$interpro_lines" | cut -f6 | sort -u | grep -v "^$" | tr '\n' ';' | sed 's/;$//')
        interpro_ids=$(echo "$interpro_lines" | cut -f12 | sort -u | grep -v "^$\|^-$" | tr '\n' ';' | sed 's/;$//')
        go_terms=$(echo "$interpro_lines" | cut -f14 | sort -u | grep -v "^$\|^-$" | tr '\n' ';' | sed 's/;$//')
    else
        domains=""
        interpro_ids=""
        go_terms=""
    fi
    
    # dbCAN annotations
    if [ -f "dbcan_annotation/overview.tsv" ]; then
        cazyme=$(grep "^$gene_id\s" dbcan_annotation/overview.tsv | cut -f2 | head -1)
    else
        cazyme=""
    fi
    
    # Output line
    echo -e "$gene_id\t$swissprot_hit\t$swissprot_desc\t$swissprot_evalue\t$domains\t$interpro_ids\t$go_terms\t$cazyme"
    
done < ./master_annotation/all_genes.txt >> ./master_annotation/master_annotations.tsv

echo "Master annotation table created!"

# Print summary statistics
echo ""
echo "=== Annotation Summary ==="
total_genes=$(wc -l < ./master_annotation/all_genes.txt)
swissprot_hits=$(cut -f2 ./master_annotation/master_annotations.tsv | grep -v "^$\|^swissprot_hit$" | wc -l)
domain_hits=$(cut -f5 ./master_annotation/master_annotations.tsv | grep -v "^$\|^domains$" | wc -l)
go_hits=$(cut -f7 ./master_annotation/master_annotations.tsv | grep -v "^$\|^go_terms$" | wc -l)
cazyme_hits=$(cut -f8 ./master_annotation/master_annotations.tsv | grep -v "^$\|^cazyme_family$" | wc -l)

echo "Total genes: $total_genes"
echo "Genes with SwissProt hits: $swissprot_hits"
echo "Genes with InterPro domains: $domain_hits"
echo "Genes with GO terms: $go_hits"
echo "Genes with CAZyme annotation: $cazyme_hits"

echo ""
echo "Done! Output: master_annotations.tsv"


# Create a visualization package directory
mkdir -p genome_package

# Copy/compress genome assembly
cp flye_assembly/ragtag_output/ragtag.scaffold.fasta genome_package/
bgzip genome_package/ragtag.scaffold.fasta
mamba run -n samtools samtools faidx genome_package/ragtag.scaffold.fasta.gz

# Create an enriched GFF3 with all annotations
cat > create_annotated_gff3.sh << 'EOF'
#!/bin/bash

# Add functional annotations to GFF3
awk 'BEGIN{FS=OFS="\t"}
NR==FNR {
    # Load annotations
    if(NR>1) {
        gene=$1
        swissprot=$3
        domains=$5
        go=$7
        cazyme=$8
        annot[gene]["product"] = swissprot
        annot[gene]["domains"] = domains
        annot[gene]["go"] = go
        annot[gene]["cazyme"] = cazyme
    }
    next
}
/^#/ {print; next}
{
    # Parse GFF3 line
    if($9 ~ /ID=/) {
        split($9, attrs, ";")
        for(i in attrs) {
            if(attrs[i] ~ /^ID=/) {
                split(attrs[i], id_parts, "=")
                gene_id = id_parts[2]
                break
            }
        }
        
        # Add annotations if available
        new_attrs = $9
        if(gene_id in annot) {
            if(annot[gene_id]["product"] != "")
                new_attrs = new_attrs ";product=" annot[gene_id]["product"]
            if(annot[gene_id]["go"] != "")
                new_attrs = new_attrs ";Ontology_term=" annot[gene_id]["go"]
            if(annot[gene_id]["cazyme"] != "")
                new_attrs = new_attrs ";Note=CAZyme:" annot[gene_id]["cazyme"]
        }
        $9 = new_attrs
    }
    print
}' master_annotation/master_annotations.tsv braker_annotation/braker.gff3 > genome_package/genome_annotated.gff3

echo "Annotated GFF3 created!"
EOF

chmod +x create_annotated_gff3.sh
./create_annotated_gff3.sh

# Sort and compress GFF3
mamba create -n bedtools -c bioconda bedtools

mamba run -n bedtools bedtools sort -i genome_package/genome_annotated.gff3 > genome_package/genome_annotated.sorted.gff3
bgzip genome_package/genome_annotated.sorted.gff3
tabix -p gff genome_package/genome_annotated.sorted.gff3.gz