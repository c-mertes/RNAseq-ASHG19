library(data.table)
library(magrittr)
library(ensemblVEP)

devtools::load_all("../tMAE/")
counts_1 <- fread('/s/project/drop-analysis/processed_data/mae/allelic_counts/HG00106--HG00106.4.M_120208_5.csv.gz')
counts_2 <- fread('/s/project/drop-analysis/processed_data/mae/allelic_counts/HG00111--HG00111.2.M_111215_4.csv.gz')
counts_1[position == 120917382 & MAE_ID == 'HG00106--HG00106.4.M_120208_5', c('altCount', 'refCount', 'totalCount', 'rawDepth') := list(33, 0, 33, 33)]
counts_2[position == 97373558 & MAE_ID == 'HG00111--HG00111.2.M_111215_4', c('altCount', 'refCount', 'totalCount', 'rawDepth') := list(193, 5, 198, 198)]

DIR <- '/s/project/drop-analysis/processed_data/mae/allelic_counts/'
mae_counts_all <- lapply(list.files(DIR, full.names = T)[1:10], function(r){
  fread(r, fill = T)
}) %>% rbindlist()
mae_counts_all[position == 120917382 & MAE_ID == 'HG00106--HG00106.4.M_120208_5', c('altCount', 'refCount', 'totalCount', 'rawDepth') := list(33, 0, 33, 33)]
mae_counts_all[position == 97373558 & MAE_ID == 'HG00111--HG00111.2.M_111215_4', c('altCount', 'refCount', 'totalCount', 'rawDepth') := list(193, 5, 198, 198)]

fwrite(mae_counts_all[MAE_ID %in% c('HG00106--HG00106.4.M_120208_5', 'HG00111--HG00111.2.M_111215_4')], '/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/mae/allelic_counts.tsv', sep = '\t', quote = F, row.names = F, na = NA)
setorder(mc, contig)
fwrite(mc[MAE_ID %in% c('HG00106--HG00106.4.M_120208_5', 'HG00111--HG00111.2.M_111215_4'), .(contig, position, variantID, refAllele, altAllele, refCount, altCount, totalCount, lowMAPQDepth, lowBaseQDepth, 
                                                                                             rawDepth, otherBases, improperPairs, MAE_ID)], 
       '/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/mae/allelic_counts.tsv', sep = '\t', quote = F, row.names = F, na = NA)

vcf_file <- '/s/project/ashg19_rnaseq_workschop/Data/input_data/variants/1000G_subset_exome.vep.vcf.gz'
# samples(scanVcfHeader(vcf_file))
vcf <- readVcf(vcf_file, param=ScanVcfParam(samples=c("HG00106", "HG00111")))

vcf <- vcf[geno(vcf)$GT[,1] %in% c("1|0", "0|1") | geno(vcf)$GT[,2] %in% c("1|0", "0|1")]
csq_dt <- parseCSQToGRanges(vcf) %>% as.data.table()

# Subset for relevant columns only
csq_dt <- csq_dt[, .(seqnames, start, end, Existing_variation, Consequence, SYMBOL, Gene, gnomAD_AF, MAX_AF)]
# Remove duplicated rows
csq_dt <- unique(csq_dt)
# Subset for variants that are annotated to a gene
csq_dt <- csq_dt[grep('^ENSG', Gene)]
# Remove double consequences
csq_dt[, Consequence := gsub('&.+', '', Consequence)]
csq_dt <- unique(csq_dt)

ranked_consequences <- c("splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", 
                         "stop_lost", "start_lost", "inframe_insertion", "missense_variant", "protein_altering_variant",
                         "incomplete_terminal_codon_variant", "stop_retained_variant", "synonymous_variant", "coding_sequence_variant",
                         "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant", 
                         "intron_variant", "non_coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant")

csq_dt[, Consequence := factor(Consequence, levels = ranked_consequences)]
setorderv(csq_dt, c('seqnames', 'start', 'Consequence'))
csq_dt <- csq_dt[, .SD[1], by = .(seqnames, start)]


fwrite(csq_dt, '/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/mae/annot_info.tsv', sep = '\t', quote = F, row.names = F, na = NA)

mc <- mc[, 1:14]
resMAE <- DESeq4MAE(mc[MAE_ID %in% c('HG00106--HG00106.4.M_120208_5', 'HG00111--HG00111.2.M_111215_4')], minCoverage = 10)
resMAE[, signif := padj < .05]
resMAE[, signif_ALT := signif == TRUE & altRatio >= .8]
resMAE[signif == TRUE, .N, by = MAE_ID]
resMAE[signif_ALT == TRUE, .N, by = MAE_ID]

resAnnot <- merge(resMAE, csq_dt, by.x = 'variantID', by.y = 'Existing_variation')
resAnnot[, MAX_AF := as.numeric(MAX_AF)]
csq_dt
ee
sample1 <- 'HG00106--HG00106.4.M_120208_5'
plotMA(resAnnot[MAE_ID == sample1], rare_column = 'rare', title = sample1)
sample2 <- 'HG00111--HG00111.2.M_111215_4'
plotMA(resAnnot[MAE_ID == sample2], rare_column = 'rare', title = sample2) 


ggplot(mc, aes(refCount, altCount)) + geom_point() + scale_y_log10() + scale_x_log10() + theme_bw()
