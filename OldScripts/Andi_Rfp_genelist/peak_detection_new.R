#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO 
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Peak-Analysis data of Ribo-seq made by Chew et al 2013
# Data analysis for Andi
# Analysis done: August 2020
# Output 3 .csv files.
# 1. Peaks per gene, all positions, zscore. Filter out bottom 50% genes by reads (must be minimum 10 reads on gene).
# 2. Max zscore peak per gene, from list 1. 
# 3. Peaks per gene, filtered by zscore > median of list 2.

library(ORFikPipeline)
library(ggplot2)
output_folder <- "/export/valenfs/projects/Hakon/Andi_genes_zebrafish/peaks/"
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load and prepare cds
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
## First load NGS experiments (ORFik experiment class)
# Chew 13
stages_c = c("2to4Cell", "256Cell", "1KCell", "Dome"); # Pick you stage here, then just run the rest. "2to4Cell", "256Cell"
df.zf.rfp <- read.experimentl(paste0("zf_Chew13")); df.zf.rfp
df.zf.rfp <- df.zf.rfp[df.zf.rfp$stage %in% stages_c,]; df.zf.rfp
df.zf.rna <- read.experimentl(paste0("zf_Chew_RNA")); df.zf.rna
df.zf.rna <- df.zf.rna[df.zf.rna$stage %in% stages_c,]; df.zf.rna
## Subtelny
df.rfp.su <- read.experimentl(paste0("zf_subtel_RFP")); df.rfp.su
df.rfp.su <- df.rfp.su[df.rfp.su$rep == 3,]; df.rfp.su
df.rna.su <- read.experimentl(paste0("zf_subtel_RNA")); df.rna.su
df.rna.su <- df.rna.su[df.rna.su$rep == 3,]; df.rna.su

# path to danrer10 gtf
gtf <- df.zf.rfp@txdb # Get gtf from experiment
txdb <- loadTxdb(gtf) 
# Get longest isoform per gene, CDS > 30 bp
txNames <- filterTranscripts(txdb, 0, 30, 0, longestPerGene = TRUE) 
loadRegions(txdb, c("mrna", "cds"), names.keep = txNames)
# Remove 10 bases from 3' end
cds <- startRegion(cds, upstream = 0, downstream = widthPerGroup(cds, FALSE) - 10 - 1)
# convert names to GENE_IDs
names(cds) <- ORFik:::txNamesToGeneNames(names(cds), txdb)
names(mrna) <- ORFik:::txNamesToGeneNames(names(mrna), txdb)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load gene list
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
gl <- "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/aligned_Zv11/aligned/QC_STATS/clusterMZgenes.txt"
gl <- data.table::fread(gl)
gl <- gl[gl$ensembl_gene_id != "", ] # only valid Ensembl ids
gl <- gl[!duplicated(gl$ensembl_gene_id), ]
print(gl[, .N, by = cluster])
# Now subset cds & mrna by gene list
cds <- cds[names(cds) %in% gl$ensembl_gene_id]; print(length(cds))
mrna <- mrna[names(mrna) %in% gl$ensembl_gene_id]; print(length(mrna))
identical(names(cds), names(mrna)) # If false something is wrong

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Stage loop
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# for each stage 
peak_analysis <- function(cds, mrna, df.rfp, df.rna, output_folder, gl) {
  if (!dir.exists(output_folder)) 
    dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
  
  dt_all <- data.table(); dt_max <- data.table(); dt_counts <- data.table(); dt_counts_all <- data.table()
  # output libs
  outputLibs(df.rfp, type = "pshifted") # output wanted cell stage
  outputLibs(df.rna, type = "bedo") # output wanted cell stage
  stages <- df.rfp$stage
  i <- 1 # lib count
  for (stage in stages) {
    print(paste("Running stage:", stage))
    if (!file.exists(paste0(output_folder,"peaks_per_gene_zscore>median_", stage,".csv"))) {
      ## Get libs
      RFP <- get(bamVarName(df.rfp)[i], envir = .GlobalEnv)
      RNA <- get(bamVarName(df.rna)[i], envir = .GlobalEnv)
      #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
      # All peaks and all positions, filtered by max(quantile(sum_per_gene, 0.50), 10) (list 1)
      #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
      # Filter out bottom 50% genes by reads (must be minimum 10 reads on gene). 
      coverage <- findPeaksPerGene(cds, RFP, type = "all")
      #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
      # Max zscore peak per gene of list 1 (list 2 & 5)
      #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
      coverage_max <- findPeaksPerGene(cds, RFP, type = "max")
      
      RNA <- convertToOneBasedRanges(RNA, addScoreColumn = TRUE)
      # List 4 here ->
      dt <- translationalEff(cds, RNA, RFP, mrna, with.fpkm = TRUE,
                             weight.RFP = "score", weight.RNA = "score")
      dt$gene_id <- names(cds)
      fwrite(dt, paste0(output_folder, "FPKM_tables_all_", stage, ".csv"))
      
      ff <- data.table::merge.data.table(coverage_max, dt, by = "gene_id"); dim(dt); dim(ff)
      print(cor.test(ff$fpkmRFP, ff$fpkmRNA))
      #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
      # Peaks per gene, removed all filtered by zscore over median (list 3)
      #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
      # 
      cov <- copy(coverage)
      cov <- cov[zscore > median(ff$zscore), ]
      ff_group <- cov[, .(peaks = .N), by = genes]
      ff_group$gene_id <- names(cds)[ff_group$genes];
      if (anyNA(ff_group$gene_id)) stop("NA found as gene_id, check naming!")
      ff_group$genes <- NULL
      #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
      # Peaks per gene, filtered by zscore over median, all kept (list 4)
      #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
      # 
      cov <- copy(coverage)
      cov[, is_peak := zscore > median(ff$zscore)]
      ff_per <- cov[, .(peaks = sum(is_peak)), by = genes]
      ff_per$gene_id <- names(cds)[ff_per$genes];
      if (anyNA(ff_per$gene_id)) stop("NA found as gene_id, check naming!")
      ff_per$genes <- NULL
      
      ## Save
      fwrite(coverage, paste0(output_folder, "Peaks_per_gene_zscore_", stage, ".csv"))
      fwrite(ff, paste0(output_folder, "Max_peak_per_gene_zscore_", stage, ".csv"))
      fwrite(ff_group, paste0(output_folder,"peaks_per_gene_zscore>median_", stage,".csv"))
      fwrite(ff_per, paste0(output_folder,"peaks_per_gene_All_kept_zscore>median_", stage,".csv"))
    } else { # Load pre-existing
      coverage <- fread(paste0(output_folder, "Peaks_per_gene_zscore_", stage, ".csv"))
      ff <- fread(paste0(output_folder, "Max_peak_per_gene_zscore_", stage, ".csv"))
      ff_group <- fread(paste0(output_folder,"peaks_per_gene_zscore>median_", stage,".csv"))
      ff_per <- fread(paste0(output_folder,"peaks_per_gene_All_kept_zscore>median_", stage,".csv"))
    }
    coverage$stage <- stage; ff$stage <- stage; ff_per$stage <- stage; ff_group$stage <- stage
    dt_all <- rbindlist(list(dt_all, coverage))
    dt_max <- rbindlist(list(dt_max, ff))
    dt_counts <- rbindlist(list(dt_counts, ff_group))
    dt_counts_all <- rbindlist(list(dt_counts_all, ff_per))
    i <- 1 + 1
  }
  
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # Plots
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  
  {
    # list 1 plot
    dt_merge <- data.table::merge.data.table(dt_all, gl, by.x = "gene_id", by.y = "ensembl_gene_id")
    value <- dt_merge
    gg_peaks <- ggplot(value, aes(x = cluster, y = zscore, group = cluster)) +
      geom_violin(aes(color = cluster)) +
      geom_boxplot(aes(color = cluster), alpha = 0.5) +
      facet_wrap( ~ stage, scales = "free") +
      ylab("zscore of all peaks"); gg_peaks
    ggsave(paste0(output_folder, "All_peaks_list1_plot_zscore",".png"), width = 12, height = 9, 
           plot = gg_peaks)
    # list 2 plot
    dt_merge <- data.table::merge.data.table(dt_max, gl, by.x = "gene_id", by.y = "ensembl_gene_id")
    value <- dt_merge
    gg_peaks <- ggplot(value, aes(x = cluster, y = zscore, group = cluster)) +
      geom_violin(aes(color = cluster)) +
      geom_boxplot(aes(color = cluster), alpha = 0.5) +
      facet_wrap( ~ stage, scales = "free") +
      ylab("zscore of max peaks"); gg_peaks
    ggsave(paste0(output_folder, "Max_plot_zscore", stage, ".png"), width = 12, height = 9, 
           plot = gg_peaks)
    # list 3 plot
    dt_merge <- data.table::merge.data.table(dt_counts, gl, by.x = "gene_id", by.y = "ensembl_gene_id")
    value <- dt_merge
    value$cluster <- as.factor(value$cluster)
    value$peak_group <- as.factor(value$peaks)
    value[, sum_peak_cluster := sum(peaks), by = .(cluster, stage)]
    value[, genes_cluster := .N , by = .(cluster, stage)]
    value[, genes_per_peak := .N, by = .(cluster, peak_group, stage)]
    value <- value[, .(sum_peaks_per_cluster = sum(peaks)), by = .(cluster, peak_group, genes_cluster, genes_per_peak, stage)]
    value <- value[order(cluster, stage),]
    value <- value[, norm := 100*(genes_per_peak / genes_cluster)]
    gg_peaks <- ggplot(value, aes(x=peak_group, y = norm, group = cluster, fill = cluster)) +
      geom_bar(aes(color = cluster), stat="identity", position=position_dodge())+
      scale_fill_brewer(palette="Paired")+
      ylab("% of genes in peak group") + 
      xlab("Peak Group") +
      facet_grid(  ~ stage) +
      theme_minimal(); gg_peaks
    ggsave(paste0(output_folder, "peaks_per_gene_plot_zscore_all_stages.png"), width = 12, height = 9, 
           plot = gg_peaks)
    # list 4 plot
    dt_merge <- data.table::merge.data.table(dt_counts_all, gl, by.x = "gene_id", by.y = "ensembl_gene_id")
    value <- dt_merge
    value$cluster <- as.factor(value$cluster)
    value$peak_group <- as.factor(value$peaks)
    value[, sum_peak_cluster := sum(peaks), by = cluster]
    value[, genes_cluster := .N , by = cluster]
    value[, genes_per_peak := .N, by = .(cluster, peak_group)]
    value <- value[, .(sum_peaks_per_cluster = sum(peaks)), by = .(cluster, peak_group, genes_cluster, genes_per_peak)]
    value <- value[order(cluster),]
    value <- value[, norm := 100*(genes_per_peak / genes_cluster)]
    gg_peaks <- ggplot(value, aes(x=peak_group, y = norm, group = cluster, fill = cluster)) +
      geom_bar(aes(color = cluster), stat="identity", position=position_dodge())+
      scale_fill_brewer(palette="Paired")+
      ylab("% of genes in peak group") + 
      xlab("Peak Group") +
      theme_minimal()
    ggsave(paste0(output_folder, "peaks_per_all_genes_plot_zscore_merged_stages.png"), width = 12, height = 9, 
           plot = gg_peaks)
    # list 5 plot: TE plots
    files <- dir(output_folder, pattern = "FPKM_tables_all", full.names = TRUE)
    dt_merge <- rbindlist(lapply(files, function(x) cbind(fread(x), stage = x)))
    dt_merge$stage <- gsub(paste0(output_folder, "/FPKM_tables_all_"), "", dt_merge$stage)
    dt_merge$stage <- gsub("\\.csv", "", dt_merge$stage)
    dt_merge <- data.table::merge.data.table(dt_merge, gl, by.x = "gene_id", by.y = "ensembl_gene_id")
    value <- melt(dt_merge, id.vars = c("gene_id", "external_gene_name", "cluster", "stage"))
    gg_te <- ggplot(value, aes(x = variable, y = value, group = cluster)) +
      geom_boxplot(aes(color = cluster)) +
      facet_wrap(  ~ stage + variable, scales = "free", ncol = 3) +
      scale_y_log10() + 
      ylab("log10(value)"); gg_te
    ggsave(paste0(output_folder, "TE&FPKM_plot_all_stages.png"), width = 12, height = 9, 
           plot = gg_te)
    # Split in 3 
    for (i in unique(value$variable)) {
      v <- value[variable == i,]
      gg_seq <- ggplot(v, aes(x = stage, y = value, group = cluster)) +
        geom_boxplot(aes(color = cluster)) +
        facet_grid(  ~ stage, scales = "free") +
        scale_y_log10() + 
        ylab(paste0("log10(", i, ")")); gg_seq
      ggsave(paste0(output_folder, i,"_plot_all_stages", ".png"), width = 12, height = 9, 
             plot = gg_seq)
    }
  }
}
# Chew 13
output_folder <- "/export/valenfs/projects/Hakon/Andi_genes_zebrafish/peaks/chew/"
peak_analysis(cds, mrna, df.zf.rfp, df.zf.rna, output_folder, gl)
# Subtelny
output_folder <- "/export/valenfs/projects/Hakon/Andi_genes_zebrafish/peaks/subtelny/"
peak_analysis(cds, mrna, df.rfp.su, df.rna.su, output_folder, gl)
