library(data.table)
library(magrittr)
DIR <- '/s/public_webshare/public/genetic_diagnosis_demo/processed_data/mae/allelic_counts'

r <- list.files(DIR, full.names = T)[1]
mae_reads <- lapply(list.files(DIR, full.names = T), function(r){
    sample <- unlist(strsplit(unlist(strsplit(r, split = '/'))[9], split = ".csv.gz"))[1]
    mr <- fread(r)
    mr[, sample := sample]
    return(mr)
}) %>% rbindlist()

fwrite(mae_reads, 'Data/input_data/variants/chr1_allellic_counts.tsv', sep = '\t', row.names = F, quote = F)


mae21 <- fread('Data/input_data/variants/chr21_allellic_counts.tsv')
mae1 <- fread('Data/input_data/variants/chr1_allellic_counts.tsv')

# Insert corrupted counts
x1 <-  data.table('chr3', 119234712, 'rs781525096', 'A', 'G', 3, 92, 95, 0, 0, 95, 0, 0, "HG00150--HG00150.4.M_120208_7")
x2 <-  data.table('chr22', 36876814, 'rs754022333', 'C', 'T', 1, 87, 88, 0, 0, 88, 0, 0, "HG00178--HG00178.4.M_120208_8")
x3 <-  data.table('chr10', 120917382, 'rs367932369', 'C', 'T', 0, 31, 31, 0, 0, 31, 0, 0, "HG00132--HG00132.2.M_111215_4")
x4 <-  data.table('chr10', 97373558, 'rs759267778', 'G', 'A', 73, 2, 75, 0, 1, 75, 0, 0, "HG00132--HG00132.2.M_111215_4")


mae_all <- rbind(mae21, mae1, x1, x2, x3, x4, use.names = F)

# Change other counts
mae_all[40, c('refCount', 'altCount', 'totalCount') := as.list(c(2, 101, 103))]
mae_all[140, c('refCount', 'altCount', 'totalCount') := as.list(c(5, 72, 77))]
mae_all[240, c('refCount', 'altCount', 'totalCount') := as.list(c(0, 20, 20))]
mae_all[340, c('refCount', 'altCount', 'totalCount') := as.list(c(4, 49, 53))]

setorder(mae_all, contig)

fwrite(mae_all, '/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/mae/allellic_counts.tsv', sep = '\t', row.names = F, quote = F)


###
# Create gnomad allele frequency table
###

gnomad21 <- readRDS('/s/genomes/human/hg19/gnomAD/allele_freq_gnomad_ukbb/gnomad_ukbb_21.Rds')
gnomad21 <- gnomad21[pos %in% mae21$position]

gnomad1 <- readRDS('/s/genomes/human/hg19/gnomAD/allele_freq_gnomad_ukbb/gnomad_ukbb_1.Rds')
gnomad1 <- gnomad1[pos %in% mae1$position]

gnomad_all <- rbind(gnomad1, gnomad21)
gnomad_all[, chr := paste0('chr', chr)]

y1 <- data.table(x1[, .(V1, V2, V4, V5)], 0, 0.00009558, 0, 0, 0, 0.0001946, 0, 0, 0.0001946, 0.0001946 )
y2 <- data.table(x2[, .(V1, V2, V4, V5)], 0, 0.000003995, 0, 0, 0, 0.000008857, 0, 0, 0.000008857, 0.000008857 )
y3 <- data.table(x3[, .(V1, V2, V4, V5)], 0, 0.000003978, 0, 0, 0, 0.000008798, 0, 0, 0.000008798, 0.000008798 )
y4 <- data.table(x4[, .(V1, V2, V4, V5)], 0, 0.000007955, 0, 0, 0, 0.000008793, 0, 0, 0.00003266, 0.00003266 )
gnomad_all <- rbind(gnomad_all, y1, y2, y3, y4, use.names = F)

setorder(gnomad_all, chr)

fwrite(gnomad_all, '/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/mae/AF_gnomAD.tsv', sep = '\t', row.names = F, quote = F)