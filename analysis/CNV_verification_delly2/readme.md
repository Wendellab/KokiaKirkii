



> intersect the gff of primary transcript exons with the deletions file  
> pipe out only exons that contain a deletion  
> extract chromosome, gene name  
> two warnings were manually checked, and do not overlap genes, i.e.,  

WARNING: File Gk.DEL.bed has a record where naming convention (leading zero) is inconsistent with other files:  
Chr10   37651   37707


WARNING: File Kokia.DEL.bed has a record where naming convention (leading zero) is inconsistent with other files:  
 Chr10   37565   37585


# G. kirkii

bedtools intersect -wo -a Dgenome2_13.exons.gff -b Gk.DEL.bed | cut -f1,9,13 | sed '/WARNING/d' | sed 's/__.__.//g' | sed 's/__..__.//g' | awk ' ($3+1)%3!=0 { print } ' | cut -f2 | sort | uniq | wc -l

# number of unique genes with deletions = 1225


# K. drynarioides

bedtools intersect -wo -a Dgenome2_13.exons.gff -b Kokia.DEL.bed | cut -f1,9,13 | sed '/WARNING/d' | sed 's/__.__.//g' | sed 's/__..__.//g' | awk ' ($3+1)%3!=0 { print } ' | cut -f2 | sort | uniq | wc -l

# number of unique genes with deletions = 1887
