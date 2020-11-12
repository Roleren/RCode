#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Analysis data of peaks in Ribo-seq made by Preeti, September 2020
library(ORFikPipeline)

df.rfp <- read.experimentl("Preeti_RFP")
outputLibs(df.rfp, type = "pshifted")

cds <- loadRegion(df.rfp, part = "cds", names.keep = filterTranscripts(df.rfp))

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Non coding RNAs peaks
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
gff.df <- ORFik:::importGtfFromTxdb(df.rfp)
ncNames <- unique(gff.df[grep(gff.df$transcript_biotype, pattern = "lincRNA"),]$transcript_id)
tx <- loadRegion(df.rfp, part = "tx", names.keep = ncNames)

dt_nc_wt1 <- findPeaksPerGene(tx, RFP_WT_r1, type = "maxmedian")
dt_nc_wt2 <- findPeaksPerGene(tx, RFP_WT_r2, type = "maxmedian")

nc_merged <- data.table::merge.data.table(dt_nc_wt1, dt_nc_wt2, by = c("genes", "gene_id"))
nc_merged_eq <- nc_merged[position.x == position.y,]

nrow(nc_merged_eq); nrow(nc_merged_eq) / nrow(nc_merged)
head(nc_merged_eq[order(count.x, decreasing = TRUE), ], 10)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Wild type correlation of peaks
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

dt_wt1 <- findPeaksPerGene(cds, RFP_WT_r1, type = "maxmedian")
dt_wt2 <- findPeaksPerGene(cds, RFP_WT_r2, type = "maxmedian")

# Filter to equal genes
wt_merged <- data.table::merge.data.table(dt_wt1, dt_wt2, by = c("genes", "gene_id"))
wt_merged_eq <- wt_merged[position.x == position.y,]
# Number of genes that have matching peaks.
nrow(wt_merged_eq); nrow(wt_merged_eq) / nrow(wt_merged)
head(wt_merged_eq[order(count.x, decreasing = TRUE), ], 10)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# NMDA and DHPG peaks correlation
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

dt_DHPG1 <- findPeaksPerGene(cds, RFP_DHPG_r1, type = "maxmedian")
dt_DHPG2 <- findPeaksPerGene(cds, RFP_DHPG_r2, type = "maxmedian")
dt_NMDA1 <- findPeaksPerGene(cds, RFP_NMDA_r1, type = "maxmedian")
dt_NMDA2 <- findPeaksPerGene(cds, RFP_NMDA_r2, type = "maxmedian")

DHPG_merged <- data.table::merge.data.table(dt_DHPG1, dt_DHPG2, by = c("genes", "gene_id"))
NMDA_merged <- data.table::merge.data.table(dt_NMDA1, dt_NMDA2, by = c("genes", "gene_id"))

DHPG_merged_eq <- DHPG_merged[position.x == position.y,]
NMDA_merged_eq <- NMDA_merged[position.x == position.y,]

nrow(DHPG_merged_eq); nrow(DHPG_merged_eq) / nrow(DHPG_merged)
nrow(NMDA_merged_eq); nrow(NMDA_merged_eq) / nrow(NMDA_merged)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Merge 
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# WT vs DHPG
nrow(wt_merged_eq); nrow(DHPG_merged_eq)
wt_DHPG_merged <- data.table::merge.data.table(wt_merged_eq, DHPG_merged_eq,
                                           by = c("genes", "gene_id"),
                                           suffixes = c("_wt", "_DHPG"))
nrow(wt_DHPG_merged)
wt_DHPG_merged_eq <- wt_DHPG_merged[position.x_wt == position.x_DHPG,]
nrow(wt_DHPG_merged_eq)

# WT vs NMDA
nrow(wt_merged_eq); nrow(NMDA_merged_eq)
wt_NMDA_merged <- data.table::merge.data.table(wt_merged_eq, NMDA_merged_eq,
                                           by = c("genes", "gene_id"),
                                           suffixes = c("_wt", "_NMDA"))
nrow(wt_NMDA_merged)
wt_NMDA_merged_eq <- wt_NMDA_merged[position.x_wt == position.x_NMDA,]
nrow(wt_NMDA_merged_eq)
# Create result table of matches between WT, DHPG and NMDA
res <- data.table(name = unique(df.rfp$condition),
                  total_genes = rep(length(cds), 3), 
                  valid_peak_genes = c(nrow(wt_merged), nrow(DHPG_merged), nrow(NMDA_merged)), 
                  match_wt_genes = c(nrow(wt_merged), nrow(wt_DHPG_merged), nrow(wt_NMDA_merged)),
                  match_wt_peak_pos = c(nrow(wt_merged), nrow(wt_DHPG_merged_eq), nrow(wt_NMDA_merged_eq)))
