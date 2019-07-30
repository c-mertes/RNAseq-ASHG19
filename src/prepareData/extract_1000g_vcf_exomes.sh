#!/bin/bash
#
#

set -x 

# input files/dirs
outDir="Data/input_data/variants"
injectionsVcf="${outDir}/injections.vcf.gz"
vcfRoot="/s/genomes/human/hg19/1000g/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
wes_regions="/s/project/mitoMultiOmics/db_data/enrichment_kit_files/agilent_sureselect_human_all_exon_xt_50mb_v5_plus_utrs_probe_targets.bed"
vep_cache_dir=/opt/modules/i12g/ensembl-vep/94/cachedir 
vep_plugin_dir=/opt/modules/i12g/ensembl-vep/94/cachedir/Plugins
cadd_dir=/s/genomes/human/hg19/CADD/v1.3

# samples to extract
annoFile="$(dirname ${outDir})/annotation.tsv"

# default variables
chrOfInt=`echo chr{1..21} chrX chrY`
vcfFileSuffix="phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
wes_regions_ensembl="${wes_regions%%.bed}_ensembl.bed"

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

    if [ ! -e "${outFile}.tbi" ]; then
        bcftools view $vcfFile --force-samples -s $saOfInt -R $wes_regions_ensembl -Ou \
            | bcftools sort -m 10G -T $TMP/$RANDOM -Ou \
            | bcftools view --min-ac 1 -m2 -M2 \
            | bgzip -@10 -c \
            > $outFile &&
        tabix $outFile &
    fi
done

# wait for jobs to finish
wait

# merge all vcf files into one
bcftools concat $allVcfs $injectionsVcf -o $outDir/1000G_subset_exome.vcf.gz -O z --threads 10
tabix $outDir/1000G_subset_exome.vcf.gz

# run vep on top
vep --cache --dir_cache $vep_cache_dir --port 3337 --force_overwrite \
    --compress_output bgzip --minimal --allele_number --verbose --fork 10 \
    --assembly GRCh37 --merged --dir_plugins $vep_plugin_dir --buffer_size 10000 \
    --sift b --polyphen s --total_length --numbers --symbol --hgvs --ccds \
    --uniprot --xref_refseq --af --max_af --af_gnomad --pubmed --canonical \
    --plugin CADD,$cadd_dir/whole_genome_SNVs.tsv.gz,$cadd_dir/InDels.tsv.gz \
    --format vcf --vcf --filter_common --biotype --failed 1 \
    --input_file $outDir/1000G_subset_exome.vcf.gz \
    --output_file $outDir/1000G_subset_exome.vep.vcf.gz


