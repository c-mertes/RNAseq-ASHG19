library(data.table)
library(Rsamtools)
library(BiocParallel)

register(MulticoreParam(30, 500, progress=TRUE))

raw_data <- file.path("Data", "input_data")
raw_rna_data <- file.path(raw_data, "rna-seq")
raw_gtf_data <- file.path(raw_data, "annotations")
raw_vcf_data <- file.path(raw_data, "variants")
anno_file_url <- "https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt"
anno_file <- file.path(raw_data, basename(anno_file_url))
vcf_file_url <- "https://i12g-gagneurweb.in.tum.de/public/workshops/RNAseq_ASHG19/input_data/variants/1000G_subset_exome.vcf.gz"
gtf_file_url <- "https://i12g-gagneurweb.in.tum.de/public/workshops/RNAseq_ASHG19/input_data/annotations/gencode.v29lift37.annotation.gtf.gz"


# prepare file system/folders
sapply(c(raw_data, raw_rna_data, raw_gtf_data, raw_vcf_data), function(x){
    if(!file.exists(x)){
        dir.create(x, recursive=TRUE)
    }
})

# download gtf and variant data
download.file(vcf_file_url, file.path(raw_vcf_data, basename(vcf_file_url)), method="wget", extra="-c")
download.file(gtf_file_url, file.path(raw_gtf_data, basename(gtf_file_url)), method="wget", extra="-c")

# get anno
download.file(anno_file_url, anno_file, method="wget", extra="-c")
anno <- fread(anno_file)


# create needed annotation files

# random subset of 3 main labs
useSampleIdx <- c(1, 7, 10, 14, 18, 28, 34, 48, 49, 62, 64, 67, 68, 72, 74, 97, 
        102, 105, 109, 114, 117, 120, 129, 134, 137, 146, 148, 150, 151, 153,
        156, 161, 170, 172, 174, 175, 183, 184, 190, 198, 201, 211, 213, 214,
        220, 231, 243, 249, 250, 254, 255, 256, 257, 268, 286, 289, 294, 295,
        300, 306, 307, 311, 314, 316, 317, 322, 323, 340, 345, 349, 352, 358,
        354, 356, 361, 363, 364, 365, 373, 376, 381, 390, 392, 393, 394, 395, 
        402, 403, 407, 409, 413, 415, 427, 431, 433, 439, 442, 447, 448, 454)

# sample annotation (VCF_ID, RNA_ID, co-variats, ...)
anno_final <- unique(anno[,.(RNA_ID=`Assay Name`, SEX=`Characteristics[sex]`,
        ORIGIN=`Characteristics[ancestry category]`, 
        LAB=Performer, INDIVIDUAL=`Factor Value[individual]`,
        phase1TG=`Comment[1000g Phase1 Genotypes]`,
        OUTRIDER_GROUP="all",
        PAIRED_END=TRUE, COUNT_MODE="IntersectionStrict", 
        COUNT_OVERLAP=TRUE, COUNT_AS_STRANDED="no")])[useSampleIdx,]
write.table(anno_final, file.path(raw_data, "annotation.tsv"), sep="\t", row.names=FALSE)

# file mapping (ID -> filename)
anno_filemapping <- anno_final[,.(ID=RNA_ID, ASSAY="RNA_ID",
    FILE=file.path(raw_rna_data, paste0(RNA_ID, ".bam")))]
write.table(anno_filemapping, file.path(raw_data, "filemapping.tsv"), sep="\t", row.names=FALSE)


# get bam files
res <- bplapply(anno_final[,RNA_ID], function(s){
    bam_url <- file.path(dirname(anno_file_url), "processed", paste0(s, ".bam"))
    bam_local <- file.path(raw_rna_data, paste0(s, ".bam"))
    message("download sample: ", s)
    download.file(bam_url, bam_local, method="wget", extra="--continue -nv")
    if(is.na(index(BamFile(bam_local))) || nrow(idxstatsBam(bam_local)) <= 23){
        indexBam(bam_local)
    }
})
