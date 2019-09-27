# 
# This script takes the raw OUTRIDER counts and manipulates them to incorporate 
# the corresponding outliers needed for the workshop
# 
library(OUTRIDER)

# In and output
input        <- "/s/project/ashg19_rnaseq_workschop/processed_data/v29/counts/all/total_counts.Rds"
outCounts    <- "/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/outrider/raw_counts.tsv.gz"
outResPval   <- "/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/outrider/results_pvalue.tsv.gz"
outResZscore <- "/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/outrider/results_zscore.tsv.gz"


# load the outrider object
ods <- OutriderDataSet(readRDS(input))

# Insert corrupted counts
cts <- counts(ods)
map <- mcols(ods)[,"gene_id"]
names(map) <- mcols(ods)[,"gene_name"]

# 3 outliers to detect
# Obvious:  NA18873.4.M_120208_7 (TIMMDC1)
# Many:     NA11918.1.M_111124_3 (MCOLN1)
# Exercise: HG00103.4.M_120208_3 (TXN2)
# 
cts[map["TIMMDC1"], "NA18873.4.M_120208_7"] <- 341
cts[map["MCOLN1"],  "NA11918.1.M_111124_3"] <- 16
cts[map["TXN2"],    "HG00103.4.M_120208_3"] <- 852


#'
#' Create OUTRIDER result file
#' 
txdb <- makeTxDbFromGFF("Data/input_data/annotations/gencode.v29lift37.annotation.gtf.gz")
odsInj <- OutriderDataSet(countData=cts)
odsInj <- filterExpression(odsInj, gtfFile=txdb, filterGenes=TRUE)

# visualize it
odsInj <- plotCountGeneSampleHeatmap(odsInj, normalized=FALSE)
odsInj <- plotCountCorHeatmap(odsInj, normalized=FALSE)

# run outrider stick to the colab configurations!:
# https://colab.research.google.com/drive/1_U_kK69Zh2_Yl2Ncggw0k1B2jbdmptme#scrollTo=JGNob8EEdvk1&line=2&uniqifier=1
odsInj <- OUTRIDER(odsInj, q=15, iterations=2)


#' 
#' Save output
#' 
write_tsv_gz <- function(x, file, ...){
    gzf <- gzfile(file, "w")
    write.table(x, gzf, sep="\t", ...)
    close(gzf)
}

write_tsv_gz(cts, outCounts)
write_tsv_gz(results(odsInj), outResPval)
write_tsv_gz(results(odsInj, padjCut=1, zScoreCut=2), outResZscore)

resInj <- results(odsInj)
resInj[geneID == map["TIMMDC1"]]
resInj[geneID == map["MCOLN1"]]
resInj[geneID == map["TXN2"]]

