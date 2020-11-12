#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Analysis of clusters from Andi using Chew data, September 2020

# Chew data (add dome stage)
# TE per stage per cluster plot
# 1st: fpkm rna and rfp per stage per cluster
# uORF analysis, seperate by cluster, make statistics (eivind)
# top motif per cluster (eivind)
# mean coverage plot bazzini
# start region relative, why 0 ?

# New analysis 7th of october:
# 1. Gene length per fraction
# 2. New CAI algorithm than MILC
# 3. Ribo collision, LSU, around 60 in readlength
# 4. ?

library(ORFikPipeline); library(ggplot2)

output_folder <- "/export/valenfs/projects/Hakon/Andi_genes_zebrafish/Feature_tables/"

#oe.dir <- "/export/valenfs/data/processed_data/experiment_tables_for_R/"
## Chew 13
stages = c("2to4Cell", "256Cell", "1KCell", "Dome"); # Pick you stage here, then just run the rest. "2to4Cell", "256Cell"
df.zf.rfp <- read.experiment(paste0("zf_Chew13")); df.zf.rfp
df.zf.rfp <- df.zf.rfp[df.zf.rfp$stage %in% stages,]; df.zf.rfp
df.zf.rna <- read.experiment(paste0("zf_Chew_RNA")); df.zf.rna
df.zf.rna <- df.zf.rna[df.zf.rna$stage %in% stages,]; df.zf.rna
## Subtelny
df.rfp.su <- read.experiment("zf_subtel_RFP"); df.rfp.su
df.rfp.su <- df.rfp.su[df.rfp.su$rep == 3,]; df.rfp.su
df.rna.su <- read.experiment(paste0("zf_subtel_RNA")); df.rna.su
df.rna.su <- df.rna.su[df.rna.su$rep == 3,]; df.rna.su

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load annotation
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# path to danrer10 gtf
gtf <- df.zf.rfp@txdb # Get gtf from experiment
txdb <- loadTxdb(gtf) 
fa <- ORFik:::findFa(df.zf.rfp)
# Get longest isoform per gene, CDS > 30 bp
txNames <- filterTranscripts(txdb, longestPerGene = TRUE) 
loadRegions(txdb, c("mrna", "cds"), names.keep = txNames)
cds_1 <- cds; names(cds_1) <- paste0(names(cds_1), "_1"); cds_1 <- ORFik:::removeMetaCols(cds_1)
geneNames <- ORFik:::txNamesToGeneNames(names(cds), txdb)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load gene list and create clusters
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

# Load list
gl <- "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/aligned_Zv11/aligned/QC_STATS/clusterMZgenes.txt"
gl <- data.table::fread(gl)
gl <- gl[gl$ensembl_gene_id != "", ] # only valid Ensembl ids
gl <- gl[gl$ensembl_gene_id %in% geneNames, ]
gl <- gl[!duplicated(ensembl_gene_id)]; print(nrow(gl))
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Compute features per stage
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#



#stages <- "Dome"
computeFeaturesByExperiment <- function(grl, df.rfp, df.rna, output_folder, gl) {
  ORFik:::validateExperiments(df.rfp)
  if (!dir.exists(output_folder)) 
    dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
  
  stages <- df.rfp$stage
  libs <-bplapply(stages,
    function(stage, grl, gl, df.rfp, df.rna, output_folder) {
      message(stage)
      
      print(paste("Running stage:", stage))
      # output libs
      savename <- paste0(output_folder, "Features_table_", df.rfp@experiment, stage, ".csv")
      if (!file.exists(savename)) {
        rfp.file <- filepath(df.rfp[df.rfp$stage == stage,], "pshifted")
        if (length(rfp.file) > 1) stop("Wrong subsetting of rfp experiment!")
        
        RFP <- fimport(rfp.file, grl)
        if (!is.null(df.rna)) {
          ORFik:::validateExperiments(df.rna)
          RNA <- fimport(filepath(df.rna[df.rna$stage == stage,], "bedo"), grl)
        } else RNA <- NULL
        
        dt <- computeFeatures(grl = grl, RFP = RFP, RNA = RNA, Gtf = txdb, 
                              sequenceFeatures = FALSE, grl.is.sorted = TRUE,
                              weight.RFP = "score", weight.RNA = "score")  
        dt <- cbind(geneNames, dt)
        
        # Now subset cds & mrna by gene list
        dt_sub <- dt[geneNames %in% gl$ensembl_gene_id,]; print(nrow(dt_sub))
        dt_merge <- data.table::merge.data.table(dt_sub, gl, by.x = "geneNames", by.y = "ensembl_gene_id")
        fwrite(dt_merge, file = savename)
      } else dt_merge <- fread(savename)
      #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
      # Cluster analysis
      #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
      value <- melt(dt_merge, id.vars = c("geneNames", "external_gene_name", "cluster"))
      value$value <- log2(value$value)
      gg_ngs <- ggplot(value, aes(x = variable, y = value, group = cluster)) +
        geom_boxplot(aes(color = cluster)) +
        facet_wrap( ~ variable, scales = "free") + 
        ylab("log2(value)"); plot(gg_ngs)
      ggsave(paste0(output_folder, "Feature_plot_", df.rfp@experiment,  stage, ".png"), 
             plot = gg_ngs, width = 12, height = 9)
      for (i in colnames(dt_merge)) {
        if (!is(unlist(dt_merge[,i, with = F]), "character")) {
          print(i)
          temp <- dt_merge[, .(i = mean(get(i), na.rm = TRUE)), by = cluster]
          print(temp[order(cluster),])
        }
      }
      return(invisible(NULL))
}, grl = grl, gl = gl, df.rfp = df.rfp, df.rna = df.rna, output_folder = output_folder)
}
output_folder_ngs <- "/export/valenfs/projects/Hakon/Andi_genes_zebrafish/Feature_tables/NGS_features/"
# Chew
computeFeaturesByExperiment(cds_1, df.zf.rfp, df.zf.rna, output_folder_ngs, gl)
# Subtelny
computeFeaturesByExperiment(cds_1, df.rfp.su, df.rna.su, output_folder_ngs, gl)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# TIS windows +/- 30 (per cluster per stage) 
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
### CHEW
cluster <- gl[, ensembl_gene_id, by = cluster]
cluster <- split(cluster$ensembl_gene_id, cluster$cluster)
TIS_region <- cds
TIS_region <- startRegion(TIS_region, mrna, TRUE, 29, 30); unique(widthPerGroup(TIS_region))
names(TIS_region) <- geneNames

dt <- data.table()

for (stage in stages) {
  remove.experiments(df.zf.rfp[df.zf.rfp$stage == stage,])
  outputLibs(df.zf.rfp[df.zf.rfp$stage == stage,], type = "pshifted") # output wanted cell stage
  clu <- 1
  for (i in cluster) {
    grl <- TIS_region[i]
    temp <- metaWindow(RFP, grl, scoring = NULL, zeroPosition = 30,
                       feature = as.character(clu), fraction = stage)
    dt <- rbindlist(list(dt, temp))
    # zscore gives shape, a good starting plot
    clu <- clu + 1
  }
}

gg_tis <- windowCoveragePlot(dt, scoring = "zscore", title = "TIS: Ribo-seq metaplot"); gg_tis
ggsave(paste0(output_folder, "TIS_windows_by_cluster_zscore", ".png"), width = 12, height = 9, 
       plot = gg_tis)
gg_tis <- windowCoveragePlot(dt, scoring = "mean", title = "TIS: Ribo-seq metaplot"); gg_tis
ggsave(paste0(output_folder, "TIS_windows_by_cluster_mean", ".png"), width = 12, height = 9, 
       plot = gg_tis)

### BAZZINI 
df.rfp.bz <- read.experimentl(paste0("zf_baz14_RFP")); df.rfp.bz
if (gtf != df.rfp.bz@txdb) stop("Not equal txdb!")

dt <- data.table()

for (stage in df.rfp.bz$stage) {
  remove.experiments(df.rfp.bz[df.rfp.bz$stage == stage,])
  outputLibs(df.rfp.bz[df.rfp.bz$stage == stage,], type = "pshifted") # output wanted cell stage
  clu <- 1
  for (i in cluster) {
    grl <- TIS_region[i]
    temp <- metaWindow(RFP, grl, scoring = NULL, zeroPosition = 30,
                       feature = as.character(clu), fraction = stage)
    dt <- rbindlist(list(dt, temp))
    # zscore gives shape, a good starting plot
    clu <- clu + 1
  }
}
gg_tis <- windowCoveragePlot(dt, scoring = "zscore", title = "TIS: Ribo-seq metaplot"); gg_tis
ggsave(paste0(output_folder, "TIS_windows_by_cluster_zscore_bazzini", ".png"), width = 12, height = 9, 
       plot = gg_tis)
gg_tis <- windowCoveragePlot(dt, scoring = "mean", title = "TIS: Ribo-seq metaplot"); gg_tis
ggsave(paste0(output_folder, "TIS_windows_by_cluster_mean_bazzini", ".png"), width = 12, height = 9, 
       plot = gg_tis)
### Subtelny
# Pshift
df.rfp.su <- read.experimentl(paste0("zf_subtel_RFP")); df.rfp.su
if (gtf != df.rfp.su@txdb) stop("Not equal txdb!")
# convertLibs(df.rfp.su); remove.experiments(df.rfp.su)
# shiftFootprintsByExperiment(df.rfp.su); shifts.load(df.rfp.su)
# shiftPlots(df.rfp.su, output = paste0(dirname(df.rfp.su$filepath[1]), "/QC_STATS/heatmaps/pshifting.png"))
df.rfp.su <- df.rfp.su[df.rfp.su$rep == 3,]; df.rfp.su

# TIS plot
dt <- data.table()

for (stage in df.rfp.su$stage) {
  remove.experiments(df.rfp.su[df.rfp.su$stage == stage,])
  outputLibs(df.rfp.su[df.rfp.su$stage == stage,], type = "pshifted") # output wanted cell stage
  clu <- 1
  for (i in cluster) {
    grl <- TIS_region[i]
    temp <- metaWindow(RFP, grl, scoring = NULL, zeroPosition = 30,
                       feature = as.character(clu), fraction = stage)
    dt <- rbindlist(list(dt, temp))
    # zscore gives shape, a good starting plot
    clu <- clu + 1
  }
}
gg_tis <- windowCoveragePlot(dt, scoring = "zscore", title = "TIS: Ribo-seq metaplot"); gg_tis
ggsave(paste0(output_folder, "TIS_windows_by_cluster_zscore_subtelny", ".png"), width = 12, height = 9, 
       plot = gg_tis)
gg_tis <- windowCoveragePlot(dt, scoring = "mean", title = "TIS: Ribo-seq metaplot"); gg_tis
ggsave(paste0(output_folder, "TIS_windows_by_cluster_mean_subtelny", ".png"), width = 12, height = 9, 
       plot = gg_tis)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# TTS windows +/- 30 (per cluster per stage) 
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
cluster <- gl[, ensembl_gene_id, by = cluster]
cluster <- split(cluster$ensembl_gene_id, cluster$cluster)
TTS_region <- cds
TTS_region <- stopCodons(TTS_region, TRUE)
TTS_region <- startRegion(TTS_region, mrna, TRUE, 29, 30); unique(widthPerGroup(TIS_region))
names(TTS_region) <- geneNames

# TTS plot
dt <- data.table()
for (stage in df.rfp.su$stage) {
  remove.experiments(df.rfp.su[df.rfp.su$stage == stage,])
  outputLibs(df.rfp.su[df.rfp.su$stage == stage,], type = "pshifted") # output wanted cell stage
  clu <- 1
  for (i in cluster) {
    grl <- TTS_region[i]
    temp <- metaWindow(RFP, grl, scoring = NULL, zeroPosition = 30,
                       feature = as.character(clu), fraction = stage)
    dt <- rbindlist(list(dt, temp))
    # zscore gives shape, a good starting plot
    clu <- clu + 1
  }
}
gg_tis <- windowCoveragePlot(dt, scoring = "zscore", title = "TTS: Ribo-seq metaplot"); gg_tis
ggsave(paste0(output_folder, "TTS_windows_by_cluster_zscore_subtelny", ".png"), width = 12, height = 9, 
       plot = gg_tis)
gg_tis <- windowCoveragePlot(dt, scoring = "mean", title = "TTS: Ribo-seq metaplot"); gg_tis
ggsave(paste0(output_folder, "TTS_windows_by_cluster_mean_subtelny", ".png"), width = 12, height = 9, 
       plot = gg_tis)


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# sequence features & codon optimality
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
savename_seqs <- paste0(output_folder, "sequence_features.csv")
if (!file.exists(savename_seqs)) {
  library(coRdon)
  cds_gl <- cds[geneNames %in% gl$ensembl_gene_id]; length(cds_gl)
  cds_gl <- ORFik:::removeMetaCols(cds_gl)
  # Full length
  
  seqs <- startRegionString(groupGRangesBy(stopSites(cds_gl, TRUE, TRUE, TRUE)),
                            tx = cds, upstream = 59, downstream = 0, faFile = fa)
  seqs <- DNAStringSet(seqs, use.names = TRUE)
  scores <- data.table(codon_bias = coRdon::MILC(coRdon::codonTable(seqs))); colnames(scores) <- "codon_bias_20AA" 
  seqs <- startRegionString(groupGRangesBy(stopSites(cds_gl, TRUE, TRUE, TRUE)),
                            tx = cds, upstream = 299, downstream = 0, faFile = fa)
  seqs <- DNAStringSet(seqs, use.names = TRUE)
  scores[, codon_bias_100AA := coRdon::MILC(coRdon::codonTable(seqs))]
  # Full length
  seqs <- ORFik:::txSeqsFromFa(cds_gl, faFile = fa, is.sorted = TRUE)
  scores[, codon_bias_full := coRdon::MILC(coRdon::codonTable(seqs))]
  
  scores[, kozak := kozakSequenceScore(cds_gl, mrna, fa)]
  scores[, gc := gcContent(cds_gl, fa)]
  # Start and stop codons
  starts <- startCodons(cds_gl, is.sorted = TRUE)
  stops <- stopCodons(cds_gl, is.sorted = TRUE)
  scores[, StartCodons := txSeqsFromFa(starts, fa, TRUE, FALSE)]
  scores[, StopCodons := txSeqsFromFa(stops, fa, TRUE, FALSE)]
  scores[, fractionLengths := fractionLength(cds_gl, widthPerGroup(mrna, TRUE))]
  scores[, length := widthPerGroup(cds_gl, FALSE)]
  scores <- cbind(geneNames = geneNames[geneNames %in% gl$ensembl_gene_id], scores)
  scores_merged <- data.table::merge.data.table(scores, gl, by.x = "geneNames", by.y = "ensembl_gene_id")
  fwrite(scores_merged, file = savename_seqs)
} else scores_merged <- fread(savename_seqs)


for (i in colnames(scores_merged)) {
  if (!is(unlist(scores_merged[,i, with = F]), "character")) {
    print(i)
    temp <- scores_merged[, .(i = median(get(i), na.rm = TRUE)), by = cluster]
    print(temp[order(cluster),])
  }
}


value <- melt(scores_merged, id.vars = c("geneNames", "StartCodons", "StopCodons", "external_gene_name", 
                                         "cluster"))
gg_seq <- ggplot(value, aes(x = variable, y = value, group = cluster)) +
    geom_boxplot(aes(color = cluster)) +
    facet_wrap( ~ variable, scales = "free") +
    scale_y_log10(); gg_seq
ggsave(paste0(output_folder, "Sequence_Feature_plot_Z10.png"), width = 12, height = 9, 
       plot = gg_seq)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Cluster size tests
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
library("biomaRt")
gll <- "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/aligned_Zv11/aligned/QC_STATS/clusterMZgenes.txt"
ensemblsIDS <- data.table::fread(gll)
ensemblsIDS[, .N, cluster] # 323
# Z11
tt <- read.experimentl("Chew_13_RFP_Z11")
n <- ORFik:::txNamesToGeneNames(names(loadRegion(tt, "cds")), tt)
glll <- ensemblsIDS[ensemblsIDS$ensembl_gene_id %in% n, ]
glll[, .N, cluster] # 303
# Z11 (filter out short trailers and leaders)
txNames <- filterTranscripts(tt, 1, 1, 1, longestPerGene = TRUE) 
ncds <- loadRegion(tt, "cds", txNames)
nmrna <- loadRegion(tt, "mrna", txNames)
n <- ORFik:::txNamesToGeneNames(names(loadRegion(tt, "cds", txNames)), tt)
names(ncds) <- n; names(nmrna) <- n
glll <- ensemblsIDS[ensemblsIDS$ensembl_gene_id %in% n, ]
glll[, .N, cluster] # 246
gll_some <- glll
# Z11 (filter out short trailers and leaders)
txNames <- filterTranscripts(tt,0, 1, 0, longestPerGene = TRUE) 
n <- ORFik:::txNamesToGeneNames(names(loadRegion(tt, "cds", txNames)), tt)
glll <- ensemblsIDS[ensemblsIDS$ensembl_gene_id %in% n, ]
glll[, .N, cluster] # 246
gll_none <- glll
gll_res <- gll_none[!(gll_none$ensembl_gene_id %in% gll_some$ensembl_gene_id),]
fwrite(gll_res, file = "no leaders_and_or_trailers.csv")
# Z10 (all)
n <- ORFik:::txNamesToGeneNames(names(loadRegion(df.zf.rfp, "cds")), df.zf.rfp)
glll <- ensemblsIDS[ensemblsIDS$ensembl_gene_id %in% n, ]
glll[, .N, cluster] # 282
glll[ensembl_gene_id != "", .N, cluster]
# Z10 (filter out short trailers and leaders)
txNames <- filterTranscripts(df.zf.rfp, longestPerGene = TRUE) 
n <- ORFik:::txNamesToGeneNames(names(loadRegion(df.zf.rfp, "cds", txNames)), df.zf.rfp)
glll <- ensemblsIDS[ensemblsIDS$ensembl_gene_id %in% n, ]
glll[, .N, cluster] # 245
# Z10 (filter out short trailers and leaders)
txNames <- filterTranscripts(df.zf.rfp, 1, 1, 1, longestPerGene = TRUE) 
n <- ORFik:::txNamesToGeneNames(names(loadRegion(df.zf.rfp, "cds", txNames)), df.zf.rfp)
glll <- ensemblsIDS[ensemblsIDS$ensembl_gene_id %in% n, ]
glll[, .N, cluster] # 260
loadRegions(df.zf.rfp, names.keep = txNames)
dt <- data.table(leaders = widthPerGroup(leaders, F), cds = widthPerGroup(cds, F),
                 trailers = widthPerGroup(trailers, F), mrna = widthPerGroup(mrna, F), txNames = names(mrna))
dt$gene_id <- ORFik:::txNamesToGeneNames(dt$txNames, df.zf.rfp)

dtt <- data.table::merge.data.table(dt, glll, by.x = "gene_id", by.y = "ensembl_gene_id")
value <- melt(dtt, id.vars = c("gene_id", "external_gene_name", "txNames", 
                                         "cluster"))
gg_seq <- ggplot(value, aes(x = variable, y = value, group = cluster)) +
  geom_boxplot(aes(color = cluster)) +
  facet_wrap( ~ variable, scales = "free") +
  scale_y_log10(); gg_seq


#
ensemblsIDS <- ensemblsIDS$external_gene_name
ensembl = useMart("refseq",dataset="drerio_gene_ensembl")
res <- getBM(attributes= c('refseq', 'hgnc_symbol'), 
      filters = 'hgnc_symbol', 
      values = ensemblsIDS, 
      mart = ensembl)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# TE (Translational efficiency) Multiple experiments
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# plot 5' 3'
# tcp (3 stages)
# giraldes (RFP and RNA, 2, 4, 6 stages)
# bazzini 14 TE (total) test

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Bazzini 14
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
df <- read.experimentl("zf_baz14_RNA")
dt_RNA <- countTable(df, type = "fpkm")
df <- read.experimentl("zf_baz14_RFP")
dt <- countTable(df, type = "fpkm")
dt <- (dt + 0.00001) / (dt_RNA + 0.00001)
dt$txNames <- names(loadRegion(df, "mrna"))
dt$gene_id <- ORFik:::txNamesToGeneNames(dt$txNames, df)

gll <- "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/aligned_Zv11/aligned/QC_STATS/clusterMZgenes.txt"
ensemblsIDS <- data.table::fread(gll)
glll <- ensemblsIDS[ensemblsIDS$ensembl_gene_id %in% dt$gene_id, ]
dt_merge <- data.table::merge.data.table(dt, glll, by.x = "gene_id", by.y = "ensembl_gene_id")
value <- melt(dt_merge, id.vars = c("gene_id", "external_gene_name", "cluster", "txNames"))
gg_seq <- ggplot(value, aes(x = variable, y = value, group = cluster)) +
  geom_boxplot(aes(color = cluster)) +
  facet_wrap( ~ variable, scales = "free") +
  scale_y_log10(); gg_seq
ggsave("bazzini_14_TE.png", plot = gg_seq, 
       width = 12, height = 9)

# bazzini 14 RFP (total) test
df <- read.experimentl("zf_baz14_RFP")
dt <- countTable(df, type = "fpkm")
dt$txNames <- names(loadRegion(df, "mrna"))
dt$gene_id <- ORFik:::txNamesToGeneNames(dt$txNames, df)

gll <- "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/aligned_Zv11/aligned/QC_STATS/clusterMZgenes.txt"
ensemblsIDS <- data.table::fread(gll)
glll <- ensemblsIDS[ensemblsIDS$ensembl_gene_id %in% dt$gene_id, ]
dt_merge <- data.table::merge.data.table(dt, glll, by.x = "gene_id", by.y = "ensembl_gene_id")
value <- melt(dt_merge, id.vars = c("gene_id", "external_gene_name", "cluster", "txNames"))
gg_seq <- ggplot(value, aes(x = variable, y = value, group = cluster)) +
  geom_boxplot(aes(color = cluster)) +
  facet_wrap( ~ variable, scales = "free") +
  scale_y_log10(); gg_seq
ggsave("bazzini_14_RFP.png", plot = gg_seq, 
       width = 12, height = 9)
# bazzini 14 RNA (total) test
df <- read.experimentl("zf_baz14_RNA")
dt <- countTable(df, type = "fpkm")
dt$txNames <- names(loadRegion(df, "mrna"))
dt$gene_id <- ORFik:::txNamesToGeneNames(dt$txNames, df)

gll <- "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/aligned_Zv11/aligned/QC_STATS/clusterMZgenes.txt"
ensemblsIDS <- data.table::fread(gll)
glll <- ensemblsIDS[ensemblsIDS$ensembl_gene_id %in% dt$gene_id, ]
dt_merge <- data.table::merge.data.table(dt, glll, by.x = "gene_id", by.y = "ensembl_gene_id")
value <- melt(dt_merge, id.vars = c("gene_id", "external_gene_name", "cluster", "txNames"))
gg_seq <- ggplot(value, aes(x = variable, y = value, group = cluster)) +
  geom_boxplot(aes(color = cluster)) +
  facet_wrap( ~ variable, scales = "free") +
  scale_y_log10(); gg_seq
ggsave("bazzini_14_RNA.png", plot = gg_seq, 
       width = 12, height = 9)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Chew 13 TE
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# chew 13 (total) test
df <- read.experimentl("zf_Chew_RNA")
dt <- countTable(df, type = "fpkm")
dt$txNames <- names(loadRegion(df, "mrna"))
dt$gene_id <- ORFik:::txNamesToGeneNames(dt$txNames, df)

gll <- "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/aligned_Zv11/aligned/QC_STATS/clusterMZgenes.txt"
ensemblsIDS <- data.table::fread(gll)
glll <- ensemblsIDS[ensemblsIDS$ensembl_gene_id %in% dt$gene_id, ]
dt_merge <- data.table::merge.data.table(dt, glll, by.x = "gene_id", by.y = "ensembl_gene_id")
value <- melt(dt_merge, id.vars = c("gene_id", "external_gene_name", "cluster", "txNames"))
gg_seq <- ggplot(value, aes(x = variable, y = value, group = cluster)) +
  geom_boxplot(aes(color = cluster)) +
  facet_wrap( ~ variable, scales = "free") +
  scale_y_log10(); gg_seq
ggsave("pauli_13.png", plot = gg_seq, 
       width = 12, height = 9)

df <- read.experimentl("zf_Chew_RNA")
dt <- countTable(df, type = "fpkm")
dt$txNames <- names(loadRegion(df, "mrna"))
dt$gene_id <- ORFik:::txNamesToGeneNames(dt$txNames, df)

# TE test
df <- read.experimentl("zf_Chew_RNA")
dt_rna <- countTable(df, type = "fpkm")
df.zf.rfp <- read.experimentl(paste0("zf_Chew13"))
dt_rfp <- countTable(df.zf.rfp, type = "fpkm")

dt <- (dt_rfp + 0.00001) / (dt_rna + 0.00001)

dt$txNames <- names(loadRegion(df, "mrna"))
dt$gene_id <- ORFik:::txNamesToGeneNames(dt$txNames, df)

gg_seq <- ggplot(value, aes(x = variable, y = value, group = cluster)) +
  geom_boxplot(aes(color = cluster)) +
  facet_wrap( ~ variable, scales = "free") +
  scale_y_log10(); gg_seq
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Lee 13 TE
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# lee 13 (total) test
df <- read.experimentl("Lee13")
dt <- countTable(df, type = "fpkm")
dt$txNames <- names(loadRegion(df, "mrna"))
dt$gene_id <- ORFik:::txNamesToGeneNames(dt$txNames, df)

gll <- "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/aligned_Zv11/aligned/QC_STATS/clusterMZgenes.txt"
ensemblsIDS <- data.table::fread(gll)
glll <- ensemblsIDS[ensemblsIDS$ensembl_gene_id %in% dt$gene_id, ]
dt_merge <- data.table::merge.data.table(dt, glll, by.x = "gene_id", by.y = "ensembl_gene_id")
value <- melt(dt_merge, id.vars = c("gene_id", "external_gene_name", "cluster", "txNames"))
gg_seq <- ggplot(value, aes(x = variable, y = value, group = cluster)) +
  geom_boxplot(aes(color = cluster)) +
  facet_wrap( ~ variable, scales = "free") +
  scale_y_log10(); gg_seq
ggsave("lee_13.png", plot = gg_seq, 
       width = 12, height = 9)

# lee 13 (polyA) test
#create.experimentl("/export/valenfs/data/processed_data/RNA-seq/lee_2013_zebrafish/polyA_RNA/aligned_GRCz10/", exper = "Lee13_polyA", viewTemplate = FALSE)
df <- read.experimentl("Lee13_polyA")
makeSummarizedExperimentFromBam(df, saveName = "/export/valenfs/data/processed_data/RNA-seq/lee_2013_zebrafish/polyA_RNA/aligned_GRCz10/QC_STATS/countTable_mrna")
dt <- countTable(df, type = "fpkm")
dt$txNames <- names(loadRegion(df, "mrna"))
dt$gene_id <- ORFik:::txNamesToGeneNames(dt$txNames, df)

gll <- "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/aligned_Zv11/aligned/QC_STATS/clusterMZgenes.txt"
ensemblsIDS <- data.table::fread(gll)
glll <- ensemblsIDS[ensemblsIDS$ensembl_gene_id %in% dt$gene_id, ]
dt_merge <- data.table::merge.data.table(dt, glll, by.x = "gene_id", by.y = "ensembl_gene_id")
value <- melt(dt_merge, id.vars = c("gene_id", "external_gene_name", "cluster", "txNames"))
gg_seq <- ggplot(value, aes(x = variable, y = value, group = cluster)) +
  geom_boxplot(aes(color = cluster)) +
  facet_wrap( ~ variable, scales = "free") +
  scale_y_log10(); gg_seq
ggsave("lee_13_polyA.png", plot = gg_seq, 
       width = 12, height = 9)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# TE bazzini 12
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
star.index <- "/export/valenfs/data/references/Zv10_zebrafish_allsteps"
align <- STAR.align.folder(input.dir = "/export/valenfs/data/raw_data/Ribo-Seq/bazzini_2012_zebrafish/merged",
                           output.dir = "/export/valenfs/data/processed_data/Ribo-seq/bazzini_2012_zebrafish/aligned_Zv10/",
                           index.dir = star.index,
                           adapter.sequence = "ATCTCGTATGCCGTCTTCTGCTTG", steps = "tr-ge",
                           min.length = 20, star.path = STAR.install(version = "2.7.0c"))
create.experimentl(paste0(align, "/aligned/"), exper = "zf_baz12_RFP", viewTemplate = FALSE)
df <- read.experimentl("zf_baz12_RFP")
ORFikQC(df)
shiftFootprintsByExperiment(df)
shifts.load(df)
shiftPlots(df, output = paste0(dirname(df$filepath[1]), "/QC_STATS/heatmaps/pshifts.png"),
           scoring = "transcriptNormalized")

create.experimentl("/export/valenfs/data/processed_data/RNA-seq/bazzini_2012_zebrafish/aligned_GRCz10/merge", exper = "zf_baz12_RNA", viewTemplate = FALSE)

df <- read.experimentl("zf_baz12_RNA")
makeSummarizedExperimentFromBam(df, paste0(dirname(df$filepath[1]), "/QC_STATS/countTable_mrna"))
dt_RNA <- countTable(df, type = "fpkm")
df <- read.experimentl("zf_baz12_RFP")
shiftFootprintsByExperiment(df)
dt <- countTable(df, type = "fpkm")
dt <- (dt + 0.00001) / (dt_RNA + 0.00001)
dt$txNames <- names(loadRegion(df, "mrna"))
dt$gene_id <- ORFik:::txNamesToGeneNames(dt$txNames, df)

gll <- "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/aligned_Zv11/aligned/QC_STATS/clusterMZgenes.txt"
ensemblsIDS <- data.table::fread(gll)
glll <- ensemblsIDS[ensemblsIDS$ensembl_gene_id %in% dt$gene_id, ]
dt_merge <- data.table::merge.data.table(dt, glll, by.x = "gene_id", by.y = "ensembl_gene_id")
value <- melt(dt_merge, id.vars = c("gene_id", "external_gene_name", "cluster", "txNames"))
gg_seq <- ggplot(value, aes(x = variable, y = value, group = cluster)) +
  geom_boxplot(aes(color = cluster)) +
  facet_wrap( ~ variable, scales = "free") +
  scale_y_log10(); gg_seq
ggsave("bazzini_12_TE.png", plot = gg_seq, 
       width = 12, height = 9)

all <- c("zebrafish_bazzini2012_2h_2", "zebrafish_bazzini2012_2h", "zebrafish_bazzini2012_4h_2",
"zebrafish_bazzini2012_4h", "zebrafish_bazzini2012_6h_2", "zebrafish_bazzini2012_6h")

