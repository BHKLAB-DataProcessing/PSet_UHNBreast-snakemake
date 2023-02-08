library(downloader)
options(timeout=1000)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- paste0(args[1], "download")
rna_seq_data_file <- args[2]

# out_dir <- '/Users/minoru/Code/bhklab/DataProcessing/PSet/getUHNBreast_2019/download'
# rna_seq_data_file <- 'Kallisto_0.46.1.tar.gz'

basePath <- "https://orcestradata.blob.core.windows.net/uhn/2019/Drug_Recomputed"
download(file.path(basePath, "UHN_recomputed.RData"), destfile = file.path(out_dir, "UHN_recomputed.RData"))

basePath <- "https://orcestradata.blob.core.windows.net/uhn/2019/RNA-seq"
download(file.path(basePath, rna_seq_data_file), destfile = file.path(out_dir, rna_seq_data_file))
