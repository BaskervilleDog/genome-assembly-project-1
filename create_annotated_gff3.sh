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
