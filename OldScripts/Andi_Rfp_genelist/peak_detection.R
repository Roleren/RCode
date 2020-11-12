library(ORFik)
library(data.table)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load and prepare cds
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# path to danrer11 gtf
gtf <- "/export/valenfs/data/references/Zv11_zebrafish/primary_assembly/Danio_rerio.GRCz11.101_ensembl.gtf.db"
txdb <- loadTxdb(gtf) 
# Get longest isoform per gene, CDS > 30 bp
txNames <- filterTranscripts(txdb, 0, 30, 0, longestPerGene = TRUE) 
loadRegions(txdb, "cds", names.keep = txNames)
# Remove 10 bases from 3' end
cds <- startRegion(cds, upstream = 0, downstream = widthPerGroup(cds, FALSE) - 10 - 1)
# convert names to GENE_IDs
names(cds) <- ORFik:::txNamesToGeneNames(names(cds), txdb)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load gene list
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
gl <- "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/aligned_Zv11/aligned/QC_STATS/clusterMZgenes.txt"
gl <- data.table::fread(gl)
gl <- gl[gl$ensembl_gene_id != "", ] # only valid Ensembl ids

# Now subset cds by gene list
cds <- cds[names(cds) %in% gl$ensembl_gene_id]; print(length(cds))

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Load reads (wig files, forward and reverse)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
reads <- fimport(c("/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/aligned_Zv11/aligned/QC_STATS/Chew13_RFP_1K_pshifted_forward.wig",
                   "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/aligned_Zv11/aligned/QC_STATS/Chew13_RFP_1K_pshifted_reverse.wig"))


# The peak detection, this is in ORFik now, but lets do it all here:
coverage <- coveragePerTiling(cds, reads, TRUE, as.data.table = TRUE)
coverage[, sum_per_gene := sum(count), by = genes]
coverage <- coverage[sum_per_gene >= max(quantile(sum_per_gene, 0.50), 10),]
coverage[, mean_per_gene := mean(count), by = genes]
coverage[, sd_per_gene := sd(count), by = genes]
coverage[, zscore := (count - mean_per_gene) / sd_per_gene, by = genes]
summary(coverage$zscore)
coverage[, gene_id := names(cds)[genes]]
coverage <- coverage[coverage[, .I[zscore == max(zscore)], by = genes]$V1]
length(unique(coverage$genes))

fwrite(coverage, "~/peak_coverage_Chew13_1KCell.csv")

# optional, filter out low max peaks
coverage_filtered <- coverage[count > 20,] # some number here
length(unique(coverage_filt$genes))