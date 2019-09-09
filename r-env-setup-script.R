#' 
#' This script sets up the R environment in the 
#' Google Colab virtual maschine. Mainly it 
#' downloads the precompiles R packages and 
#' installs needed Ubuntu packages.
#' Furthermore it downloads all needed datasets
#' to run the workshop analysis pipelines.
#' 

#' install needed ubuntu packages
system("apt update")
system("apt install -y libmysqlclient-dev libudunits2-dev libgeos-dev libgdal-dev libcairo2-dev")

#' Function to download bigger files from google drive automatically.
downloadGDriveFile <- function(id, out){
    system(paste0(
        'wget --continue --load-cookies /tmp/cookies.txt ',
            '"https://docs.google.com/uc?export=download&confirm=$(wget ',
            '--quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate ',
            "'https://docs.google.com/uc?export=download&id=", id, "' -O- | ",
            "sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\\1\\n/p')&id=", id, '" ',
            "-O ", out, " && rm -rf /tmp/cookies.txt"), intern=TRUE)
}

#' R package cache (precompiled packages based on colab env)
downloadGDriveFile(out="r_binaries.tar.gz", id="1-VqZb_Cv63AH8ogWwhSp48DUowzv_GQf")

#' Unpack cache locally
system("tar -xzf r_binaries.tar.gz -C /", intern=TRUE)

#' Download count matrix and sample annotation
system("wget --continue https://i12g-gagneurweb.in.tum.de/public/workshops/RNAseq_ASHG19/input_data/annotation.tsv")
system("wget --continue -P outrider https://i12g-gagneurweb.in.tum.de/public/workshops/RNAseq_ASHG19/input_data/outrider/raw_counts.tsv.gz")
system("wget --continue -P annotations https://i12g-gagneurweb.in.tum.de/public/workshops/RNAseq_ASHG19/input_data/annotations/gencode.v29lift37.annotation.txdb")

#' Set correct library path
.libPaths(Sys.getenv("R_LIBS_USER"))
