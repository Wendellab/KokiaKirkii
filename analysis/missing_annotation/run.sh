#!/usr/bin/env bash

cd $PBS_O_WORKDIR

PATH=$PATH:/usr/local/igbb/blast-2.5.0+/bin/
PATH=$PATH:/usr/local/igbb/bedtools-2.25.0/bin/
PATH=$PATH:/usr/local/igbb/EMBOSS-6.6.0/bin/
PATH=$PATH:/usr/local/igbb/exonerate-2.2.0/bin/
PATH=$PATH:/usr/local/igbb/gmap-2017-05-08/bin/

# Create GMAP Database for Kirkii, masking all annotated genes
awk '$3 == "gene"' analysis/annotation/annotations/kirkii.renamed.gff \
                 > analysis/gene_copy_number/kirkii.genes.gff
bedtools maskfasta -fi db/kirkii.renamed.fa \
                   -bed analysis/gene_copy_number/kirkii.genes.gff \
                   -fo analysis/gene_copy_number/kirkii.gene_mask.fa
gmap_build -D analysis/gene_copy_number/db \
           -d kirkii analysis/gene_copy_number/kirkii.gene_mask.fa


# Create GMAP Database for Kokia, masking all annotated genes
awk '$3 == "gene"' analysis/annotation/annotations/kokia.corrected.gff \
                 > analysis/gene_copy_number/kokia.genes.gff
bedtools maskfasta -fi analysis/assemble/abyss/kokia/kokia.renamed.fa \
                   -bed analysis/gene_copy_number/kokia.genes.gff \
                   -fo analysis/gene_copy_number/kokia.gene_mask.fa
gmap_build -D analysis/gene_copy_number/db \
           -d kokia analysis/gene_copy_number/kokia.gene_mask.fa


# Create GMAP Database for G. raimondii, masking all annotated genes
awk '$3 == "gene"' db/G.raimondii_JGI_221_v2.1.transcripts.gff3 \
                 > analysis/gene_copy_number/gorai.genes.gff
bedtools maskfasta -fi db/G.raimondii_JGI_221_v2.0.assembly.fasta \
                   -bed analysis/gene_copy_number/gorai.genes.gff \
                   -fo analysis/gene_copy_number/gorai.gene_mask.fa         
gmap_build -D analysis/gene_copy_number/db \
           -d gorai analysis/gene_copy_number/gorai.gene_mask.fa


#Search all given fastas against all databases, only output genes with at least
# 90% identity over 90% coverage
for db in kokia kirkii gorai ; do
    for fasta in {Kirkii,Kokia}_{loss,gain}; do
        gmap -t 20 \
             -D analysis/gene_copy_number/db/ \
             -d $db \
             -f gff3_gene \
             --min-trimmed-coverage=0.9 \
             --min-identity=0.9 \
             analysis/gene_copy_number/fastas/$fasta.fasta \
             > analysis/gene_copy_number/$fasta.$db.gmap.gff
    done
done
