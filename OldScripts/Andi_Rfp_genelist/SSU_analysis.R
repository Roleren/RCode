#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Analysis of clusters from Andi using Chew data, September 2020
# SSU data


library(ORFikPipeline); library(ggplot2)

output_folder <- "/export/valenfs/projects/Hakon/Andi_genes_zebrafish/Feature_tables/"
# create.experimentl("/export/valenfs/data/processed_data/TCP-seq/valen_2019_zebrafish_15_trim3_15nt_tRNAscan/aligned_GRCz10_tidy_tRNA_multi/merged",
#                    exper = "Valen15", viewTemplate = FALSE)
create.experimentl("/export/valenfs/data/processed_data/TCP-seq/valen_2018_zebrafish_7_trim3_15nt_tRNAscan/aligned_GRCz10_tidy_tRNA/merged",
                   exper = "Valen6M", viewTemplate = FALSE)


df.zf.ssu.15 <- read.experimentl(paste0("Valen15")); df.zf.ssu.15
#convertLibs(df.zf.ssu.15); remove.experiments(df.zf.ssu.15)
df.zf.ssu.6 <- read.experimentl("Valen6M")
# convertLibs(df.zf.ssu.6); remove.experiments(df.zf.ssu.6)
df.zf.ssu.6 <- df.zf.ssu.6[df.zf.ssu.6$libtype == "SSU",]

df.zf.ssu <- df.zf.ssu.15[df.zf.ssu.15$libtype == "SSU",]; df.zf.ssu
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load annotation
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# path to danrer10 gtf
gtf <- df.zf.ssu@txdb # Get gtf from experiment
txdb <- loadTxdb(gtf) 
fa <- ORFik:::findFa(df.zf.ssu)
# Get longest isoform per gene, CDS > 30 bp
txNames <- filterTranscripts(txdb, longestPerGene = TRUE) 
loadRegions(txdb, c("mrna", "cds"), names.keep = txNames)
cds_1 <- cds; names(cds_1) <- paste0(names(cds_1), "_1"); cds_1 <- ORFik:::removeMetaCols(cds_1)
geneNames <- ORFik:::txNamesToGeneNames(names(cds), txdb)



#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Cluster info table
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
gl <- "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/aligned_Zv11/aligned/QC_STATS/clusterMZgenes.txt"
gl <- data.table::fread(gl)
gl <- gl[gl$ensembl_gene_id != "", ] # only valid Ensembl ids
gl <- gl[gl$ensembl_gene_id %in% geneNames, ]
gl <- gl[!duplicated(ensembl_gene_id)]; print(nrow(gl))



#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# TIS windows +/- 30 (per cluster per stage) 
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Shield
# update annotation
cluster <- gl[, ensembl_gene_id, by = cluster]
cluster <- split(cluster$ensembl_gene_id, cluster$cluster)
TIS_region <- startRegion(cds, mrna, TRUE, 29, 30); unique(widthPerGroup(TIS_region))
names(TIS_region) <- geneNames

# TIS plot
dt <- data.table()
for (stage in df.zf.ssu$stage) {
  remove.experiments(df.zf.ssu[df.zf.ssu$stage == stage,])
  outputLibs(df.zf.ssu[df.zf.ssu$stage == stage,], type = "ofst") # output wanted cell stage
  clu <- 1
  for (i in cluster) {
    grl <- TIS_region[i]
    temp <- metaWindow(SSU, grl, scoring = NULL, zeroPosition = 30,
                       feature = as.character(clu), fraction = stage)
    dt <- rbindlist(list(dt, temp))
    # zscore gives shape, a good starting plot
    clu <- clu + 1
  }
}

# All stages
loadRegions(txdb, c("mrna", "cds"), names.keep = txNames)
cluster <- gl[, ensembl_gene_id, by = cluster]
cluster <- split(cluster$ensembl_gene_id, cluster$cluster)
TIS_region <- startRegion(cds, mrna, TRUE, 29, 30); unique(widthPerGroup(TIS_region))
names(TIS_region) <- geneNames

# TIS plot
dt_new <- data.table()
for (stage in df.zf.ssu.6$stage) {
  remove.experiments(df.zf.ssu.6[df.zf.ssu.6$stage == stage,])
  outputLibs(df.zf.ssu.6[df.zf.ssu.6$stage == stage,], type = "ofst") # output wanted cell stage
  clu <- 1
  for (i in cluster) {
    grl <- TIS_region[i]
    temp <- metaWindow(SSU, grl, scoring = NULL, zeroPosition = 30,
                       feature = as.character(clu), fraction = stage)
    dt_new <- rbindlist(list(dt_new, temp))
    # zscore gives shape, a good starting plot
    clu <- clu + 1
  }
}
dt_all <- rbindlist(list(dt_new, dt))
gg_tis <- windowCoveragePlot(dt_all, scoring = "zscore", title = "TIS: RCP-seq metaplot"); gg_tis
ggsave(paste0(output_folder, "TIS_windows_by_cluster_zscore_Valen", ".png"), width = 12, height = 9, 
       plot = gg_tis)
gg_tis <- windowCoveragePlot(dt_all, scoring = "mean", title = "TIS: RCP-seq metaplot"); gg_tis
ggsave(paste0(output_folder, "TIS_windows_by_cluster_mean_Valen", ".png"), width = 12, height = 9, 
       plot = gg_tis)


# whole mrna
loadRegions(df.zf.ssu, names.keep = txNames)
names(leaders) <- geneNames; names(cds) <- geneNames;names(trailers) <- geneNames;
dt <- data.table()
for (stage in df.zf.ssu$stage) {
  remove.experiments(df.zf.ssu[df.zf.ssu$stage == stage,])
  outputLibs(df.zf.ssu[df.zf.ssu$stage == stage,], type = "ofst") # output wanted cell stage
  clu <- 1
  for (i in cluster) {
    leaders_1 <- leaders[i]; cds_1 <- cds[i]; trailers_1 <- trailers[i]
    temp <- ORFik:::splitIn3Tx(leaders_1, cds_1, trailers_1,
                               SSU, fraction = clu,
                               windowSize = 30)
    dt <- rbindlist(list(dt, temp))
    # zscore gives shape, a good starting plot
    clu <- clu + 1
  }
}
gg_tis <- windowCoveragePlot(dt, scoring = "zscore", title = "TIS: RCP-seq metaplot"); gg_tis
ggsave(paste0(output_folder, "Transcript_by_cluster_zscore_Valen", ".png"), width = 12, height = 9, 
       plot = gg_tis)
gg_tis <- windowCoveragePlot(dt, scoring = "mean", title = "TIS: RCP-seq metaplot"); gg_tis
ggsave(paste0(output_folder, "Transcript_by_cluster_mean_Valen", ".png"), width = 12, height = 9, 
       plot = gg_tis)