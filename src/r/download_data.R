library(data.table)
library(Rsamtools)
library(BiocParallel)

register(MulticoreParam(30, 500, progress=TRUE))

raw_data <- file.path("Data", "raw_data")
raw_rna_data <- file.path(raw_data, "rna-seq")
anno_file_url <- "https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt"
anno_file <- file.path(raw_data, basename(anno_file_url))


# prepare file system/folders
if(!file.exists(raw_data)){
    dir.create(raw_data, recursive=TRUE)
}
if(!file.exists(raw_rna_data)){
    dir.create(raw_rna_data, recursive=TRUE)
}


# get anno
download.file(anno_file_url, anno_file, method="wget", extra="-c")
anno <- fread(anno_file)


# get bam files
res <- bplapply(anno[grepl("_1.fastq", get("Comment[FASTQ_URI]")), `Assay Name`], function(s){
    bam_url <- file.path(dirname(anno_file_url), "processed", paste0(s, ".bam"))
    bam_local <- file.path(raw_rna_data, paste0(s, ".bam"))
    message("download sample: ", s)
    download.file(bam_url, bam_local, method="wget", extra="--continue -nv")
    if(is.na(index(BamFile(bam_local))) || nrow(idxstatsBam(bam_local)) <= 23){
        indexBam(bam_local)
    }
})

# create needed annotation files
# file mapping (ID -> filename)
anno_filemapping <- anno[,.(
        ID=`Assay Name`, 
        FILE=file.path("../RNAseq-ASHG19/Data/raw_data/rna-seq", paste0(`Assay Name`, ".bam")), 
        ASSAY="RNA_ID")]
write.table(anno_filemapping, "Data/raw_data/filemapping.tsv", sep="\t", row.names=FALSE)

# sample annotation (VCF_ID, RNA_ID, co-variats, ...)
anno_final <- anno[,.(RNA_ID=`Assay Name`, OUTRIDER_GROUP="all")]
write.table(anno_final, "Data/raw_data/annotation.tsv", sep="\t", row.names=FALSE)

