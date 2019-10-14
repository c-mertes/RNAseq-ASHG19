# 
# This script takes the raw FraseR counts and manipulates them to incorporate 
# the corresponding outliers needed for the workshop
# 
devtools::load_all("../FraseR/")
library(BSgenome)
library(data.table)
register(MulticoreParam(15, 30, progressbar=TRUE))
genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg19")

# In and output
wdir         <- "/s/public_webshare/public/genetic_diagnosis_demo/processed_data/aberrant_splicing/datasets"
dname        <- "raw-all"
outJCounts   <- "/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/splicing/raw_junction_counts.tsv.gz"
outSCounts   <- "/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/splicing/raw_site_counts.tsv.gz"
outAnno      <- "/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/splicing/annotation.tsv.gz"

outResPval   <- "/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/splicing/results_pvalue.tsv.gz"
outResZscore <- "/s/public_webshare/public/workshops/RNAseq_ASHG19/input_data/splicing/results_zscore.tsv.gz"



# read annotation
anno <- fread("/s/project/drop-analysis/processed_data/aberrant_splicing/annotations/all.tsv")

# read j counts
jfiles <- list.files("/s/project/drop-analysis/processed_data/aberrant_splicing/datasets/cache/splitCounts/", full.names=TRUE)
names(jfiles) <- gsub(".RDS", "", gsub(".*-", "", jfiles))
jcounts <- mergeCounts(bplapply(jfiles, readRDS), assumeEqual=FALSE, BPPARAM=bpparam())
motif <- GenomicAlignments:::.extract_unoriented_intron_motif(genome, jcounts)
strand(jcounts) <- GenomicAlignments:::.infer_intron_strand(motif)
strand(jcounts)[strand(jcounts) == "*"] <- "+"
jcounts <- annotateSpliceSite(jcounts)


# read ss counts
ssfiles <- list.files("/s/project/drop-analysis/processed_data/aberrant_splicing/datasets/cache/nonSplicedCounts/", full.names=TRUE)
names(ssfiles) <- gsub(".RDS", "", gsub(".*-", "", ssfiles))
sscounts <- mergeCounts(bplapply(ssfiles, readRDS), assumeEqual=TRUE)

# prepare text input 
missingSamples <- anno[!anno[,sampleID] %in% colnames(mcols(sscounts)), sampleID]
imputedData <- mcols(sscounts)[,sample(3:ncol(mcols(sscounts)), length(missingSamples), replace=TRUE)]
colnames(imputedData) <- missingSamples
mcols(sscounts) <- cbind(mcols(sscounts), imputedData)

countedSamples <- colnames(mcols(sscounts))[c(-1,-2)]
newOrder <- unique(c(colnames(as.data.table(granges(jcounts))), "startID", "endID", countedSamples))
jcountsdt  <- as.data.table(jcounts)[,..newOrder]
newOrder <- unique(c(colnames(as.data.table(granges(jcounts))), "spliceSiteID", "type", countedSamples))
sscountsdt <- as.data.table(sscounts)[,..newOrder]

length(missingSamples)

#'
#' Inject outliers
#' 
#' NA11918.1.M_111124_3
#' TIMMDC1 chr3:119217368-119243937
jcountsdt[seqnames == "chr3" & start == 119232567  & end == 119234706, c("NA11918.1.M_111124_3"):=list(62)]
jcountsdt[seqnames == "chr3" & start == 119234787  & end == 119236051, c("NA11918.1.M_111124_3"):=list(62)]
jcountsdt[seqnames == "chr3" & start == 119232567  & end == 119236051, c("NA11918.1.M_111124_3"):=list(2)]

#' NA20505.1.M_111124_6
#' TAZ chrX:153638798-153651712
jcountsdt[seqnames == "chrX" & start == 153641905  & end == 153642443, c("NA20505.1.M_111124_6"):=list(28)]
jcountsdt[seqnames == "chrX" & start == 153641905  & end == 153642437, c("NA20505.1.M_111124_6"):=list(1)]

#' HG00132.2.M_111215_4
#' KCTD7 chr7:66091868-66110216
jcountsdt[seqnames == "chr7" & start == 66098432 & end == 66100508, c("HG00132.2.M_111215_4"):=list(14)]
jcountsdt[seqnames == "chr7" & start == 66100901 & end == 66103239, c("HG00132.2.M_111215_4"):=list(14)]
jcountsdt[seqnames == "chr7" & start == 66098432 & end == 66103239, c("HG00132.2.M_111215_4"):=list(0)]

#' 
#' remove names
colnames(jcountsdt) <- gsub("\\..*", "", colnames(jcountsdt))
colnames(sscountsdt) <- gsub("\\..*", "", colnames(sscountsdt))
anno[,sampleID:=gsub("\\..*", "", sampleID)]

#' 
#' Save raw output
#' 
write_tsv_gz <- function(x, file, ...){
    gzf <- gzfile(file, "w")
    write.table(x, gzf, sep="\t", ...)
    close(gzf)
}

write_tsv_gz(jcountsdt, outJCounts)
write_tsv_gz(sscountsdt, outSCounts)
write_tsv_gz(anno[condition %in% countedSamples], outAnno)


#' 
#' Load FraseR object and run full analysis
#' 
jc <- fread(outJCounts, sep="\t", drop=1)
sc <- fread(outSCounts, sep="\t", drop=1)
an <- fread(outAnno, sep="\t", drop=1)
colData <- DataFrame(an, row.names = an$sampleID)

fds <- FraseRDataSet(colData=colData, junctions=jc, spliceSites=sc, 
        parallel=bpparam())
dontWriteHDF5(fds) <- TRUE
fds

#' 
#' Filtering 
fds <- calculatePSIValues(fds)
fds <- filterExpression(fds, minDeltaPsi=0.00, filter=FALSE)
plotFilterExpression(fds, bins=100)
fds <- fds[mcols(fds)[["passed"]],]

#'
#' Fit autoencoder
correction <- "PCA"
iterations <- 2
parallel(fds) <- bpparam()
fds
for(i in psiTypes){
    message("\n", date(), ": type: ", i, "\n")
    fds <- fit(fds, correction=correction, q=12, weighted=FALSE, 
            iterations=iterations, type=i, BPPARAM=bpparam())
    
    message("\n", date(), ": pvalues: ", i, "\n")
    fds <- calculatePvalues(fds, type=i)
    
    message("\n", date(), ": padjust: ", i, "\n")
    fds <- calculatePadjValues(fds, type=i)
    
    message("\n", date(), ": zscore: ", i, "\n")
    fds <- calculateZScores(fds, type=i)
}

fds <- annotateRanges(fds)

#'
#' extract results and safe them
resp <- results(fds, fdrCutoff=0.05)
respgene <- as.data.table(resultsByGenes(resp))
write_tsv_gz(respgene, outResPval)


resz <- results(fds, fdrCutoff=NA, zScoreCutoff=2)
reszgene <- as.data.table(resultsByGenes(resz))
write_tsv_gz(reszgene, outResZscore)


#'
#' Visualize results
#' 
plotJunctionCounts(fds, "psi5", basePlot=FALSE, result=respgene[hgncSymbol == "TIMMDC1"])
plotJunctionCounts(fds, "psi5", basePlot=FALSE, result=reszgene[hgncSymbol == "TAZ"])
plotJunctionCounts(fds, "psi5", basePlot=FALSE, result=reszgene[hgncSymbol == "KCTD7"])


