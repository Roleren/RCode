#
# index.Rsubread <- function() {
#   library(Rsubread)
#   dir.create("/export/valenfs/data/references/Rnor_6.0_rat/Rsubread-index/", recursive = TRUE, showWarnings = FALSE)
#   Rsubread::buildindex(basename = "/export/valenfs/data/references/Rnor_6.0_rat/Rsubread-index/test_index",
#                        reference = "/export/valenfs/data/references/Rnor_6.0_rat/NONCODE_ncRNA_rat.fa")
#   # Took 1.5 minutes
#   dir.create("/export/valenfs/data/processed_data/Ribo-seq/testest1/", recursive = TRUE, showWarnings = FALSE)
#   Rsubread::align(index = "/export/valenfs/data/references/Rnor_6.0_rat/Rsubread-index/test_index",
#                   readfile1 = "/export/valenfs/data/raw_data/Ribo-Seq/Preeti_Jain_2020_Rattus_norvegicus/DHPG_set6_S3_R1_001.fastq.gz",
#                   output_file = "/export/valenfs/data/processed_data/Ribo-seq/testest1/DHPG_set6_S3_R1_001.bam",
#                   nthreads = 64)
# }
#
# df.rfp <- read.experimentl("Preeti_RFP")
# txdb <- loadTxdb(df.rfp)
# lincs <- ORFik:::importGtfFromTxdb(txdb)
# type <- lincs
# part <- "lincRNA"
# tx <- NULL
# valids <- type[grep(x = type$transcript_biotype, pattern = part)]
# if (length(valids) == 0) stop("found no valid transcript of type", part)
# if (is.null(tx)) tx <- loadRegion(path)
#
# tx[unique(valids$transcript_id)]
#
# ncrnas <- ORFik:::loadTranscriptType(txdb, part = "lincRNA")
