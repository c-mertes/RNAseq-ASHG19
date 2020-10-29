# 
# This script takes the raw FraseR counts and manipulates them to incorporate 
# the corresponding outliers needed for the workshop
# 
# BiocManager::install("FRASER")
library(FRASER)
package.version("FRASER")

# 
# Input
# 
fds_raw_obj  <- "/s/project/drop-analysis/processed_data/aberrant_splicing/datasets/savedObjects/raw-all/fds-object.RDS"

#
# Output
# 
outJCounts   <- "/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/splicing/raw_junction_counts.tsv.gz"
outSCounts   <- "/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/splicing/raw_site_counts.tsv.gz"
realAnno     <- "/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/annotation.tsv"

outResPval   <- "/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/splicing/results_pvalue.tsv.gz"
outResZscore <- "/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/splicing/results_zscore.tsv.gz"
outfitFds    <- "/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/splicing/fitted_fds.RDS"

# load global annotation
anno <- fread(realAnno)[,1:7]
anno

# FRASER object
fds_raw <- loadFraserDataSet(file=fds_raw_obj)
dim(fds_raw)

# adjust to global annotation
colnames(fds_raw) <- gsub("\\..*", "", colnames(fds_raw))
colData(fds_raw)[,"sampleID"] <- colnames(fds_raw)
colnames(nonSplicedReads(fds_raw)) <- colnames(fds_raw)
fds_raw[,anno[,INDIVIDUAL]]

# extract count data from fds object
jcountsdt <- cbind(
        as.data.table(rowRanges(fds_raw, type="psi5")[,c("startID", "endID")]),
        as.data.table(counts(fds_raw, type="psi5")))
sscountsdt <- cbind(
        as.data.table(rowRanges(fds_raw, type="theta")[,c("spliceSiteID", "type")]),
        as.data.table(counts(fds_raw, "theta")))

#'
#' Inject outliers
#' 
#' NA11918.1.M_111124_3
#' TIMMDC1 chr3:119217368-119243937
jcountsdt[seqnames == "chr3" & start == 119232567  & end == 119234706, c("NA11918"):=list(62)]
jcountsdt[seqnames == "chr3" & start == 119234787  & end == 119236051, c("NA11918"):=list(62)]
jcountsdt[seqnames == "chr3" & start == 119232567  & end == 119236051, c("NA11918"):=list(2)]

#' NA20505.1.M_111124_6
#' TAZ chrX:153638798-153651712
jcountsdt[seqnames == "chrX" & start == 153641905  & end == 153642443, c("NA20505"):=list(28)]
jcountsdt[seqnames == "chrX" & start == 153641905  & end == 153642437, c("NA20505"):=list(1)]

#' HG00132.2.M_111215_4
#' KCTD7 chr7:66091868-66110216
jcountsdt[seqnames == "chr7" & start == 66098432 & end == 66100508, c("HG00132"):=list(14)]
jcountsdt[seqnames == "chr7" & start == 66100901 & end == 66103239, c("HG00132"):=list(14)]
jcountsdt[seqnames == "chr7" & start == 66098432 & end == 66103239, c("HG00132"):=list(0)]


#' 
#' Save raw output
#' 
write_tsv_gz <- function(x, file, row.names=FALSE, ...){
    gzf <- gzfile(file, "w", compression=9)
    write.table(x, gzf, sep="\t", row.names=row.names, ...)
    close(gzf)
}

write_tsv_gz(jcountsdt, outJCounts)
write_tsv_gz(sscountsdt, outSCounts)


#' 
#' Load FraseR object and run full analysis
#' 
jc <- fread(outJCounts)
sc <- fread(outSCounts)
sanno <- fread(realAnno)[,1:6]
an <- colData(fds_raw)
colData <- DataFrame(row.names=sanno$INDIVIDUAL,
        merge(sanno, an, by.x="INDIVIDUAL", by.y="sampleID", sort=FALSE))
colData$sampleID <- colData$INDIVIDUAL

fds <- FraserDataSet(colData=colData, junctions=jc, spliceSites=sc)
dontWriteHDF5(fds) <- TRUE
workingDir(fds) <- "/tmp/FRASER_workshop"
fds

#' 
#' Filtering
register(MulticoreParam(15, 30, progre=TRUE))
set.seed(42)
patient_ids <- c("NA11918", "NA20505", "HG00132")
fds <- fds[,unique(c(patient_ids, sample(colnames(fds))))[1:60]]
fds <- calculatePSIValues(fds)

fds <- filterExpressionAndVariability(fds, minDeltaPsi=0.1, filter=FALSE)
plotFilterExpression(fds, bins=100)
fds <- fds[mcols(fds)[["passed"]],]
plotCountCorHeatmap(fds, "psi5", logit=TRUE, sampleClustering=NA, 
        annotation_col=c("SEX", "LAB", "ORIGIN"), topN=20000)

#'
#' Find hyper parameter (encoding dimension)
fds
qs <- integer(3)
names(qs) <- psiTypes
register(MulticoreParam(22, progre=TRUE))
for(i in psiTypes){
    message("\n", date(), ": hyper opt for type: ", i, "\n")
    fds <- optimHyperParams(fds, type=i, q_param=unique(sort(c(2, 4:10, 12, 
            seq(14, min(50, ncol(fds)), by = 3)))))
    plotEncDimSearch(fds, i)
    qs[i] <- bestQ(fds, i)
}

#' 
#' Fit FRASER model
register(MulticoreParam(10, progre=TRUE))
qs
fds <- FRASER(fds, q=qs)
fds <- annotateRangesWithTxDb(fds)

saveRDS(fds, file=outfitFds, compress="gzip")

#'
#' extract results and safe them
resp <- results(fds)
respgene <- as.data.table(resultsByGenes(resp))
write_tsv_gz(respgene, outResPval)


resz <- results(fds, fdrCutoff=NA, zScoreCutoff=2)
reszgene <- as.data.table(resultsByGenes(resz))
write_tsv_gz(reszgene, outResZscore)


#'
#' Visualize results
#' 
plotJunctionCounts(fds, "psi5", basePlot=FALSE, result=reszgene[hgncSymbol == "TIMMDC1"][1])
plotJunctionCounts(fds, "psi5", basePlot=FALSE, result=reszgene[hgncSymbol == "TAZ"][1])
plotJunctionCounts(fds, "psi5", basePlot=FALSE, result=reszgene[hgncSymbol == "KCTD7"][1])

reszgene[sampleID == "HG00132"]
sort(table(respgene[,sampleID]))
