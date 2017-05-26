### trim files ###
module load trimmomatic/0.36

for a in *_R1.fq.gz; do trimmomatic PE -threads 20 $a ${a%_R1.fq.gz}_R2.fq.gz -baseout ${a%_R1.fq.gz}.trim.fq.gz ILLUMINACLIP:Adapters.fa:2:30:15 LEADING:28 TRAILING:28 SLIDINGWINDOW:8:28 SLIDINGWINDOW:1:10 MINLEN:65 TOPHRED33; done
rename 1P R1 *
rename 2P R2 *
for a in *1U.fq.gz; do cat $a ${a%_1U.fq.gz}_2U.fq.gz > ${a%_1U.fq.gz}_S.fq.gz; done

### create bwa index and map ###
# genome from ftp://ftp.bioinfo.wsu.edu/species/Gossypium_raimondii/JGI_221_G.raimondii_Dgenome/assembly/G.raimondii_JGI_221_v2.0.assembly.fasta.gz, taking only the 13 Chromosomes
module load bwa/0.7.15
bwa index Dgenome2_13.fasta

ls *_R1.fq.gz | cut -f1 -d '.' | while read line; do bwa mem -R "@RG\\tID:$line\\tPL:ILLUMINA\\tPU:NULL\\tLB:Shotgun_LIB\\tSM:$line-NULL" -t 200 -M Dgenome2_13.fasta $line.trim_R1.fq.gz $line.trim_R2.fq.gz > $line.PE.sam; bwa mem -R "@RG\\tID:$line\\tPL:ILLUMINA\\tPU:NULL\\tLB:Shotgun_LIB\\tSM:$line-NULL" -t 20 -M Dgenome2_13.fasta $line.trim_S.fq.gz > $line.SE.sam; done

### create picard and samtools indexes required for gatk
module load picard/2.9.0
module load samtools/1.2 # avoids the error [W::bam_merge_core2] No @HD tag found

picard CreateSequenceDictionary R=Dgenome2_13.fasta O=Dgenome2_13.dict
samtools faidx Dgenome2_13.fasta

### merge and go sam2bam

for a in *.PE.bam; do samtools merge -@ 200 ${a%.PE.bam}.bam $a ${a%.PE.bam}.SE.bam; done
for a in *.bam; do samtools sort $a ${a%.bam}.sort; done

#### start picard/gatk ####
module load parallel/20160422
module rm java/1.7.0_55
module load gatk/3.6

### remove duplicates, then index 
ls *.sort.bam | parallel --jobs 2 'picard MarkDuplicates I={} O={.}.dedup.bam REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE METRICS_FILE={.}.dedup.stats &> {.}.log'
ls *.dedup.bam | parallel 'samtools index {} &> {.}.log'

### realign intervals
ls *.dedup.bam | parallel --jobs 2 'gatk -T RealignerTargetCreator -R Dgenome2_13.fasta -I {} -o {.}.realign.intervals -nt 20 2> {.}.RTC.err'
ls *.dedup.bam | parallel --jobs 2 'gatk -T IndelRealigner -R Dgenome2_13.fasta -I {} -targetIntervals {.}.realign.intervals -o {.}.realign.bam 2> {.}.indel.err'

### call variants and then joint genotyping

ls *.realign.bam | parallel --jobs 2 'gatk -T HaplotypeCaller -R Dgenome2_13.fasta -I {} -nct 20 -stand_emit_conf 10 -stand_call_conf 30 --emitRefConfidence GVCF -variant_index_type LINEAR -variant_index_parameter 128000 --genotyping_mode DISCOVERY -o {.}.raw_variants.g.vcf 2>{.}.HaplotypeCaller.err'

### join the vcf for each chromosome for all accessions

gatk -T GenotypeGVCFs -R Dgenome2_13.fasta -stand_emit_conf 10 -stand_call_conf 30 --variant Gk.raw.variants.g.vcf --variants Kok.raw.variants.g.vcf --out all.vcf 2> joint.genotype.err

### generate hard filter for SNPs and indels #runs very fast

gatk -T SelectVariants -R Dgenome2_13.fasta -V all.vcf -selectType SNP -o all.raw_snp.vcf 2> all.SelectVariants_snp.err
 
gatk -T SelectVariants -R Dgenome2_13.fasta -V all.vcf -selectType INDEL -o all.raw_indel.vcf 2> all.SelectVariants_indel.err

# apply hard filters to snp 
 
gatk -T VariantFiltration -R Dgenome2_13.fasta -V all.raw_snp.vcf --filterExpression "QD<2.0||FS>60.0||MQ<40.0||SOR>4.0||MQRankSum<-12.5||ReadPosRankSum<-8.0" --filterName all_snp_filter -o all.filtered_snp.vcf 2> all.VariantFiltration_snp.err

#apply hard filters to indel 

gatk -T VariantFiltration -R Dgenome2_13.fasta -V all.raw_indel.vcf --filterExpression "QD<2.0||FS>200.0||SOR>10.0||ReadPosRankSum<-20.0" --filterName TANG_indel_filter -o all.filtered_indel.vcf 2> all.VariantFiltration_indel.err

 
gatk -T SelectVariants -R Dgenome2_13.fasta --variant all.filtered_snp.vcf --excludeFiltered -o all.PASS_snp.vcf 2> all.PASS_snp.err

gatk -T SelectVariants -R Dgenome2_13.fasta --variant all.filtered_indel.vcf --excludeFiltered -o all.PASS_indel.vcf 2> all.PASS_indel.err

# make minimal vcf for indels
gatk -T CombineVariants -R Dgenome2_13.fasta --variant all.PASS_indel.vcf -o indel.minimal.vcf --minimalVCF

# preparing the file for R
sed '/[.][/][.]/d' indel.minimal.vcf > indel.minimal.pruned.vcf 
sed -i.bak '/0[/]1/d' indel.minimal.pruned.vcf
sed -i.bak '/1[/]2/d' indel.minimal.pruned.vcf
sed -i.bak '/0[/]2/d' indel.minimal.pruned.vcf
sed -i.bak '/0[/]3/d' indel.minimal.pruned.vcf
sed -i.bak '/2[/]3/d' indel.minimal.pruned.vcf
sed -i.bak '/1[/]3/d' indel.minimal.pruned.vcf
sed -i.bak '/1[/]4/d' indel.minimal.pruned.vcf
sed -i.bak '/1[/]2/d' indel.minimal.pruned.vcf
cut -f1,2,4,5,10,11 indel.minimal.pruned.vcf | sed '/##/d' > indels.table

