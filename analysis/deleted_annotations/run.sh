#!/usr/bin/env bash

cd $PBS_O_WORKDIR

PATH=$PATH:/usr/local/igbb/blast-2.5.0+/bin/
PATH=$PATH:/usr/local/igbb/bedtools-2.25.0/bin/
PATH=$PATH:/usr/local/igbb/EMBOSS-6.6.0/bin/
PATH=$PATH:/usr/local/igbb/exonerate-2.2.0/bin/
PATH=$PATH:/usr/local/igbb/gmap-2017-05-08/bin/
PATH=$PATH:/usr/local/igbb/gffutils/gffread-0.9.9/
swsetup bioperl

# setup databases for easy coding
ln -s  $PBS_O_WORKDIR/analysis/assemble/abyss/kokia/kokia.renamed.fa analysis/deleted_annotations/kokia.fa
ln -s  $PBS_O_WORKDIR/analysis/annotation/annotations/kokia.corrected.gff analysis/deleted_annotations/kokia.gff

ln -s  $PBS_O_WORKDIR/db/kirkii.renamed.fa analysis/deleted_annotations/kirkii.fa
ln -s  $PBS_O_WORKDIR/analysis/annotation/annotations/kirkii.renamed.gff analysis/deleted_annotations/kirkii.gff

# Create GMAP Databases
for db in kokia kirkii ; do
    gmap_build -D analysis/deleted_annotations/db \
               -d $db analysis/deleted_annotations/$db.fa
done

# Search fasta against all databases, only output genes with at least
# 90% identity over 90% coverage
for db in kokia kirkii ; do
    gmap -t 20 \
         -D analysis/deleted_annotations/db/ \
         -d $db \
         -f gff3_gene \
         --min-trimmed-coverage=0.90 \
         --min-identity=0.90 \
         --no-chimeras \
         analysis/deleted_annotations/Unclustered_Raimondii.CDS.fasta \
         > analysis/deleted_annotations/Unclustered_Raimondii.$db.gmap.gff
    awk '$3 == "mRNA"' analysis/deleted_annotations/Unclustered_Raimondii.$db.gmap.gff \
        > analysis/deleted_annotations/Unclustered_Raimondii.$db.gmap.mRNA.gff
done


# Create table and lists of Gorai genes
for db in kokia kirkii ; do
    echo $db

    # Number of Gorai genes that map to the database 
    awk '$3 == "mRNA" { print $9 }' analysis/deleted_annotations/Unclustered_Raimondii.$db.gmap.gff |
        sed -e 's/ID=//' -e 's/.mrna[0-9];.*//' |
        sort -u > analysis/deleted_annotations/$db.gmap.hit.lst

    # Number of Genes with at least on mRNA hit with a valid protein (^M[^.]*\.)
    gffread -g analysis/deleted_annotations/$db.fa -y  - \
            analysis/deleted_annotations/Unclustered_Raimondii.$db.gmap.gff | 
        bp_seqconvert.pl --from fasta --to tab |
        awk ' $2 !~ /^M[^.]*\.$/ { print $1 } ' FS="\t" |
        sed 's/.mrna[0-9]//' |
        sort -u > analysis/deleted_annotations/$db.gmap.valid_protein.lst

    # Number of genes with at least on mRNA intersecting a current gene annotation on the same strand
    awk '$3 == "gene"' analysis/deleted_annotations/$db.gff |
        bedtools intersect -a analysis/deleted_annotations/Unclustered_Raimondii.$db.gmap.mRNA.gff -b -  -s |
        sed -e 's/.*ID=//' -e 's/.mrna[0-9].*//' |
        sort -u > analysis/deleted_annotations/$db.gmap.intersect.lst

    # Number of genes with both valid protein and intersecting a currently annotated gene
    comm -12 analysis/deleted_annotations/$db.gmap.valid_protein.lst \
         analysis/deleted_annotations/$db.gmap.intersect.lst \
         > analysis/deleted_annotations/$db.gmap.both.lst

    #Display table
    wc -l  analysis/deleted_annotations/$db.gmap.hit.lst \
       analysis/deleted_annotations/$db.gmap.valid_protein.lst \
       analysis/deleted_annotations/$db.gmap.intersect.lst \
       analysis/deleted_annotations/$db.gmap.both.lst |
        sed -e '/total/d' -e 's/^  *//' -e 's/ .*//'
done |
    paste - - - - -| 
    /usr/local/igbb/abyss-2.0.1/bin/abyss-tabtomd
