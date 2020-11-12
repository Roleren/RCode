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
gl <- gl[,2:3]
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

ints <- loadRegion(txdb, "introns", by = "gene", names.keep = geneNames)
gen <- transcriptsBy(txdb, by = "gene")
gen <- gen[names(ints)]

dt <- data.table(gene_id = names(ints), gene_length = widthPerGroup(gen), introns_length_total = widthPerGroup(ints), 
                 num_introns = numExonsPerGroup(ints, FALSE))
dt_merge <- data.table::merge.data.table(dt, gl, by.x = "gene_id", by.y = "ensembl_gene_id")


value <- melt(dt_merge, id.vars = c("gene_id", "cluster"))
gg_seq <- ggplot(value, aes(x = variable, y = value, group = cluster)) +
  geom_boxplot(aes(color = cluster)) +
  facet_wrap( ~ variable, scales = "free") +
  scale_y_discrete(trans = scales::pseudo_log_trans()); gg_seq

ggsave(paste0(output_folder, "gene_lengths_clusters.png"), width = 12, height = 9, 
       plot = gg_seq)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# NEW AIC algorithms
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

library(coRdon)

# 30 AA (90 nt) 3' end
cds_gl <- cds[geneNames %in% gl$ensembl_gene_id]; length(cds_gl)
cds_gl <- ORFik:::removeMetaCols(cds_gl)

seqs <- startRegionString(groupGRangesBy(stopSites(cds_gl, TRUE, TRUE, TRUE)),
                          tx = cds, upstream = 59, downstream = 0, faFile = fa)
seqs <- DNAStringSet(seqs, use.names = TRUE)
scores <- data.table(codon_bias_20AA_MILC = coRdon::MILC(coRdon::codonTable(seqs)),
                     )

# 100 AA (300 nt) 5' end
seqs <- startRegionString(cds_gl,
                          tx = cds, upstream = 299, downstream = 0, faFile = fa)
seqs <- DNAStringSet(seqs, use.names = TRUE)
scores[, codon_bias_100AA_5p_MILC := coRdon::MILC(coRdon::codonTable(seqs))]

# 100 AA (300 nt) 3' end
seqs <- startRegionString(groupGRangesBy(stopSites(cds_gl, TRUE, TRUE, TRUE)),
                          tx = cds, upstream = 299, downstream = 0, faFile = fa)
seqs <- DNAStringSet(seqs, use.names = TRUE)
scores[, codon_bias_100AA := coRdon::MILC(coRdon::codonTable(seqs))]

# Full length
seqs <- ORFik:::txSeqsFromFa(cds_gl, faFile = fa, is.sorted = TRUE)
scores[, codon_bias_full := coRdon::MILC(coRdon::codonTable(seqs))]