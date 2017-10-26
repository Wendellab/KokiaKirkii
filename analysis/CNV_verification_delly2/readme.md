
### Make the gff into bed, because bedtools hates this gff format
`sed 's/<DEL>/0/g' Kokia.DEL.pass.vcf | awk '{gsub("A|T|G|C","1",$5)}1' | sed 's/ /\t/g' | sed '/#/d' | awk '{ OFS="\t" } { print $1, $2, $2-1+length($4)-$5 }' > Kokia.DEL.bed` 

`sed 's/<DEL>/0/g' Gk.DEL.pass.vcf | awk '{gsub("A|T|G|C","1",$5)}1' | sed 's/ /\t/g' | sed '/#/d' | awk '{ OFS="\t" } { print $1, $2, $2-1+length($4)-$5 }' > Gk.DEL.bed`


### intersect the gff of primary transcript exons with the deletions file  
pipe out only exons that contain a deletion  
extract chromosome, gene name  
two warnings were manually checked, and do not overlap genes, i.e.,  

`WARNING: File Gk.DEL.bed has a record where naming convention (leading zero) is inconsistent with other files:  `
> Chr10   37651   37707


`WARNING: File Kokia.DEL.bed has a record where naming convention (leading zero) is inconsistent with other files:  `
> Chr10   37565   37585


# _G. kirkii_

`bedtools intersect -wo -a Dgenome2_13.exons.gff -b Gk.DEL.bed | cut -f1,9,13 | sed '/WARNING/d' | sed 's/__.__.//g' | sed 's/__..__.//g' | awk ' ($3+1)%3!=0 { print } ' | cut -f2 | sort | uniq | wc -l`

### number of unique genes with disruptive deletions = 1011
number of gene deletions suggested by OrthoFinder = 2957 (includes possible whole gene deletions)
number of gene deletions verified by gmap = 1156 (out of 2144 checked)


# _K. drynarioides_

`bedtools intersect -wo -a Dgenome2_13.exons.gff -b Kokia.DEL.bed | cut -f1,9,13 | sed '/WARNING/d' | sed 's/__.__.//g' | sed 's/__..__.//g' | awk ' ($3+1)%3!=0 { print } ' | cut -f2 | sort | uniq | wc -l`

### number of unique genes with deletions = 1426
number of gene deletions suggested by OrthoFinder = 2008 (includes possible whole gene deletions)
number of gene deletions verified by gmap = 944 (out of 1458 checked)
