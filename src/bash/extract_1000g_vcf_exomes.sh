#!/bin/bash
#
#

set -x 

chrOfInt=`echo chr{1..21} chrX`
vcfRoot="/s/genomes/human/hg19/1000g/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
outDir="Data/raw_data/variants"
vcfFileSuffix="phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

wes_regions="/s/project/mitoMultiOmics/db_data/enrichment_kit_files/agilent_sureselect_human_all_exon_xt_50mb_v5_plus_utrs_probe_targets.bed"
wes_regions_ensembl="${wes_regions%%.bed}_ensembl.bed"

annoFile="Data/raw_data/annotation.tsv"

# samples of interes
saOfInt=`cat $annoFile | cut -f5 | sed 's/"//g' | tail  -n +2 | tr '\n'  ',' | sed 's/,$//'`

# change to correct chromosome names in ENSEMBL 
cat $wes_regions | sed 's/^chr//' > $wes_regions_ensembl

# extract the regions
allVcfs=""
for chr in $chrOfInt
do
    vcfFile=$vcfRoot/ALL.${chr}.*.genotypes.vcf.gz
    outFile="$outDir/`basename $vcfFile`"
    outFile="${outFile%%.vcf.gz}.exome.vcf.gz"
    allVcfs="${allVcfs} $outFile"

    echo "`date`: Running for vcf File $vcfFile into $outFile ..."

#    bcftools view $vcfFile --force-samples -s $saOfInt -R $wes_regions_ensembl \
#        | bcftools sort -m 10G -T $TMP | bgzip -@10 -c -l 9 \
#        > $outFile
#    tabix $outFile
done

# merge all vcf files into one
bcftools concat $allVcfs -o $outDir/1000G_subset_exome.vcf.gz -O z --threads 10
tabix $outDir/1000G_subset_exome.vcf.gz

