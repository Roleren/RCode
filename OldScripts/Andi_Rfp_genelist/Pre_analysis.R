#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO 
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Pre-Analysis data of Ribo-seq made by Chew et al 2013
# Analysis done: August 2020
library(ORFikPipeline)
library(ORFik)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Download annotation (Drerio 11 (93))
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Custome contaminant genomes: phix, ncRNA, (rRNA and tRNA download manually)
# Do depletion on phix, ncRNA and rRNA
annotations <- getGenomeAndAnnotation(
  organism = "Danio rerio", 
  GTF = TRUE, # Set to FALSE if you already have annotation
  genome = TRUE, # Set to FALSE if you already have genome
  output.dir = "/export/valenfs/data/references/Zv11_zebrafish/primary_assembly", 
  phix = TRUE, 
  ncRNA = "auto",
  tRNA = "", # Skip depletion of tRNA, waste of time
  rRNA = "/export/valenfs/data/references/rrna/SILVA_119_bothSURef.fasta"
)

# Index genome
index <- STAR.index(annotations, wait = FALSE)

### RNA-seq alignment
# Libraries downloaded from uninett filesender (fix your files here and put them somewhere):
input.dir.rna <- "/export/valenfs/data/raw_data/RNA-Seq/pauli_2012_zebrafish/"
output.dir.rna <- "/export/valenfs/data/processed_data/RNA-seq/pauli_2012_zebrafish/aligned_Zv11/"
STAR.align.folder(input.dir.rna, output.dir.rna, index,
                  paired.end = "yes", include.subfolders = "y",
                  steps = "ge", # Check steps argument for what this means (trim not needed, good quality)
                  adapter.sequence = "TGGAATTCTCGGGTGCCAAGG",
                  max.cpus = 90, trim.front = 3, min.length = 20)

### Ribo-seq alignment
input.dir.rfp <- "/export/valenfs/data/raw_data/Ribo-Seq/chew_2013_zebrafish"
output.dir.rfp <- "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/aligned_Zv11/"
STAR.align.folder(input.dir.rfp, output.dir.rfp, index,
                  steps = "rR-ge", # Check steps argument for what this means
                  adapter.sequence = "TGGAATTCTCGGGTGCCAAGG",
                  max.cpus = 90, trim.front = 0, min.length = 20) 
system(paste("/Home/ii/hakontj/R/x86_64-redhat-linux-gnu-library/R-3.6.0/library/ORFik/STAR_Aligner/cleanup_folders.sh", output.dir.rfp))

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Create ORFik experiments of libraries
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
txdb_file <- paste0(annotations["gtf"], ".db") # Get txdb file, not raw gtf
fa <- annotations["genome"]
create.experimentl(exper = "Chew_13_RFP_Z11",
                  dir = paste0(output.dir.rfp, "/aligned/"),
                  txdb = txdb_file, fa = fa, organism = "Danio rerio", viewTemplate = FALSE)


create.experimentl(exper = "Pauli_12_RNA_Z11",
                   dir = "/export/valenfs/data/processed_data/RNA-seq/Preeti_Jain_2020_Rattus_norvegicus/aligned",
                   txdb = txdb_file, fa = fa, organism = "Danio rerio", pairedEndBam = TRUE,
                   viewTemplate = FALSE)
# Fix columns manually in Libre office / Microsoft Office if needed

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# QC report (RNA-seq & Ribo-seq)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
df.rfp <- read.experimentl(paste0("Chew_13_RFP_Z11"))
convertLibs(df.rfp)
outputLibs(df.rfp, type = "ofst")
ORFik::QCreport(df.rfp)

df.rna <- read.experimentl("Pauli_12_RNA_Z11")
convertLibs(df.rna, type = "bedo", method = "5prime")
ORFik::QCreport(df.rna)
ORFik:::transcriptWindow1(df = df.rna, idName = "fullReads", outdir = paste0(dirname(df.rna$filepath[1]), "/QC_STATS/"))
makeSummarizedExperimentFromBam(df.rna, longestPerGene = FALSE,
                                saveName = "/export/valenfs/data/processed_data/RNA-seq/pauli_2012_zebrafish/aligned_Zv11/aligned/QC_STATS/countTable_mrna",
                                region = "mrna")

counts <- countTable(df.rna, region = "mrna", collapse = T)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# P-shifting
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
shiftFootprintsByExperiment(df.rfp, output_format = c("ofst", "wig"), accepted.lengths = 28:32)

outputLibs(df.rfp[1,], type = "ofst")
df <- detectRibosomeShifts(RFP, df.rfp, heatmap = T)
df <- data.frame(fraction = c(28, 29, 30, 31), offsets_start = c(-11, -11, -12, -12))
a <- shiftFootprints(RFP, df)
export.wiggle(a, "~/Chew13_RFP_1K.wig")
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Heatmap validation of P-shifting
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
outputLibs(df.rfp, type = "pshifted")
heatMapRegion(df.rfp, region = "TIS", acceptedLengths = NULL, scores = "transcriptNormalized")


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Peak analysis test
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
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

# Filter out bottom 50% genes by reads (must be minimum 10 reads on gene). 
#Find max zscore position per gene
outputLibs(df.rfp[1,], type = "ofst")
reads <- RFP
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



# fpkms deseq
counts_rna <- countTable(df.rna, region = "mrna", type = "fpkm", collapse = T)
counts_rfp <- countTable(df.rfp, region = "mrna", type = "fpkm", collapse = T)

dt <- data.table(rna_fpkm = counts_rna$RNA_1KCell, rfp_fpkm = counts_rfp$RFP_1KCell, 
                 a = rownames(counts_rna), b = rownames(counts_rfp))
identical(dt$a, dt$b)
dt[, te := ((rfp_fpkm + 0.000001) / (rna_fpkm + 0.000001))]

# Single isoform per gene
dt <- dt[a %in% txNames, ]
dt$a <- ORFik:::txNamesToGeneNames(dt$a, txdb)
final <- dt[a %in% coverage$gene_id,]
final[, gene_id := a, ]
final$a <- NULL; final$b <- NULL

ff <- data.table::merge.data.table(coverage, final, by = "gene_id")
fwrite(ff, "~/peak_coverage_Chew13_1KCell_with_fpkm.csv")
coverage <- coverage[count > 20,]
length(unique(coverage$genes))

# fpkm danrer10
cds2 <- loadRegion(df.zf.rna, "cds")
names(cds2) <- ORFik:::txNamesToGeneNames(names(cds2), df.zf.rna)
countOverlaps(cds2["ENSDARG00000000729"], RNA)

df.zf.rna <- read.experimentl("zf_Chew_RNA")
counts_zf_rna <- countTable(df.zf.rna, region = "mrna", type = "fpkm", collapse = T)

dt2 <- data.table(rna_fpkm2 = counts_zf_rna$RNA_1KCell,  a = rownames(counts_zf_rna))
dt2 <- dt2[a %in% txNames, ]
dt2$a <- ORFik:::txNamesToGeneNames(dt2$a, txdb)
dt2 <- dt2[a %in% coverage$gene_id,]
dt2[, gene_id := a, ]; dt2$a <- NULL


dt <- cbind(final[gene_id %in% dt2$gene_id,], dt2)
dt[,ncol(dt)] <- NULL
dt[, rna_fpkm := rna_fpkm2]; dt[, rna_fpkm2 := NULL]
dt[, te := ((rfp_fpkm + 0.000001) / (rna_fpkm + 0.000001))]
identical(dt$gene_id, dt$gene_id); summary(dt)
coverage2 <- coverage[gene_id %in% dt$gene_id,]
coverage2 <- coverage2[!duplicated(gene_id),]; any(duplicated(coverage2$gene_id))

dim(coverage2); dim(dt)
ff_new <- data.table::merge.data.table(coverage2, dt, by = "gene_id"); dim(ff_new)
fwrite(ff_new, "~/peak_coverage_Chew13_1KCell_with_fpkm_RNA10.csv")
# New mix

counts_rna <- countTable(df.zf.rna, region = "mrna", type = "fpkm", collapse = T)
counts_rfp <- countTable(df.rfp, region = "mrna", type = "fpkm", collapse = T)

dt <- data.table(rna_fpkm = counts_rna$RNA_1KCell, rfp_fpkm = counts_rfp$RFP_1KCell, 
                 a = rownames(counts_rna), b = rownames(counts_rfp))
identical(dt$a, dt$b)
dt[, rna_fpkm := rna_fpkm2]; dt[, rna_fpkm2 := NULL]

dt[, te := rfp_fpkm + 0.000001 / (rna_fpkm + 0.000001)]

# Single isoform per gene
dt <- dt[a %in% txNames, ]
dt$a <- ORFik:::txNamesToGeneNames(dt$a, txdb)
final <- dt[a %in% coverage$gene_id,]
final[, gene_id := a, ]
final$a <- NULL; final$b <- NULL

ff <- data.table::merge.data.table(coverage, final, by = "gene_id")
fwrite(ff, "~/peak_coverage_Chew13_1KCell_with_fpkm.csv")
coverage <- coverage[count > 20,]
length(unique(coverage$genes))

# Peaks zscore over median

cov<- coveragePerTiling(cds, reads, TRUE, as.data.table = TRUE)
cov[, sum_per_gene := sum(count), by = genes]
cov <- cov[sum_per_gene >= max(quantile(sum_per_gene, 0.50), 10),]
cov[, mean_per_gene := mean(count), by = genes]
cov[, sd_per_gene := sd(count), by = genes]
cov[, zscore := (count - mean_per_gene) / sd_per_gene, by = genes]
summary(cov$zscore)

cov[, gene_id := names(cds)[genes]]
cov <- cov[zscore > median(ff_new$zscore), ]
ff_group <- cov[, .(peaks = .N), by = genes]
ff_group$gene_id <- names(cds)[ff_group$genes];anyNA(ff_group$gene_id)
ff_group$genes <- NULL

max(ff_group$peaks); nrow(ff_group)

fwrite(ff_group, "~/peaks_per_gene_zscore>median.csv")

# correlation test
df.zf.rfp <- read.experimentl("zf_Chew13")
counts_rna_c <- countTable(df.zf.rna, region = "cds", type = "fpkm", collapse = T)
counts_rfp_c <- countTable(df.zf.rfp, region = "mrna", type = "fpkm", collapse = T)
cds23 <- loadRegion(df.zf.rna, "cds")
outputLibs(df.zf.rfp[1,], type = "bedo")
RFP2 <- collapseDuplicatedReads(RFP)
seqlevels(RFP2) <- seqlevels(RFP)
ov <- countOverlapsW(cds23, RFP)

cor.test(counts_rna_c$RNA_1KCell, ov)
