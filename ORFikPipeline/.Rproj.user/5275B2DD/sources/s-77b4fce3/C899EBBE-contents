#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Analysis of test data of C11 made by Yamilla, MARCH 2020
#devtools::document("/export/valenfs/projects/uORFome/RCode1/ORFikPipeline/")
############################## CREATE DATA ###################################################
library(ORFikPipeline)
plotFolder <- "/export/valenfs/projects/Hakon/C11/"

df <- read.experimentl("Val20HumMer")
outputLibs(df, type = "bedo")
dt <- readRDS("/export/valenfs/projects/Hakon/C11/matrix.rds")
# Cage
txdb <- loadTxdb("/export/valenfs/projects/Hakon/C11/GRCH38.CAGE.Mammary.gtf.db")
txNames100 <- filterTranscripts(txdb, 100, 100, 100)

# heatmaps
heatMapRegion(df, c("TIS","TSS", "TTS"))
# With cage
heatMapRegion(df, c("TIS","TSS"), paste0(dirname(df$filepath[1]), "/QC_STATS/heatmaps/TSS_cage/"),
              cage = "/export/valenfs/projects/uORFome/DATA/CAGE/human/Mammary%20Epithelial%20Cell%2c%20donor1.CNhs11077.11273-116H4.hg38.nobarcode.ctss.bed.gz")


# Merged Fraction 11-13 SSU
library(ORFikPipeline)
df <- read.experimentl("Val20HumMer")
outputLibs(df, type = "bedo")
heatMapRegion(df, c("TSS"), paste0(dirname(df$filepath[1]), "/QC_STATS/heatmaps/TSS_cage/"),
              cage = "/export/valenfs/projects/uORFome/DATA/CAGE/human/Mammary%20Epithelial%20Cell%2c%20donor1.CNhs11077.11273-116H4.hg38.nobarcode.ctss.bed.gz")


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Metacoverage, IR & TOP motif
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Analysis march 2020

# 1 Metacoverage
library(ORFikPipeline)
df <- read.experimentl("Val20HumMer", expInVarName = TRUE)
outputLibs(df, type = "bedo")

rnaFilt <- dt$txNames[(dt$RNA_MRNA_FPKM > 1 & dt$SSU_LEADERS_FPKM > 55)]
mrna <- loadRegion(txdb, "mrna", names.keep = rnaFilt)
filter <- names(mrna) %in% txNames100
mrna <- mrna[filter]

t <- removeBadTxByRegion(mrna, get(bamVarName(df)[4]))
loadRegions(txdb, names.keep = names(t))
transcriptWindow(leaders, get("cds", mode = "S4"),
                 trailers, df = df, outdir = paste0(dirname(df$filepath[1]), "/QC_STATS/metacoverage/"),
                 allTogether = TRUE, colors = c(rep("orange", 2),rep("skyblue4", 2)),
                 scores = c("sum", "zscore"), idName = "5prime_SSU55")
slackrUpload(filename = paste0(dirname(df$filepath[1]), "/QC_STATS/metacoverage/","val20humMer_cp_all_zscore_5prime.png"), channels = "visualizations")

# Use BAM files
df <- read.experimentl("Val20HumMer")
outputLibs(df, chrStyle = t)
transcriptWindow(leaders, get("cds", mode = "S4"),
                 trailers, df = df, outdir = paste0(dirname(df$filepath[1]), "/QC_STATS/metacoverage/"),
                 allTogether = TRUE, colors = c(rep("orange", 2),rep("skyblue4", 2)),
                 scores = c("sum", "zscore"), idName = "bams")
slackrUpload(filename = paste0(dirname(df$filepath[1]), "/QC_STATS/metacoverage/","val20humMer_cp_all_zscore_.png"), channels = "visualizations")

transcriptWindow(leaders, get("cds", mode = "S4"),
                 trailers, df = df[df$condition == "WT",], outdir = paste0(dirname(df$filepath[1]), "/QC_STATS/metacoverage/"),
                 allTogether = TRUE, colors = c(rep("orange", 1),rep("skyblue4", 1)),
                 scores = c("sum", "zscore"), idName = "bams_WT", )

gg_sum <- transcriptWindow(leaders, get("cds", mode = "S4"),
                 trailers, df = df[df$condition == "rp",], outdir = NULL,
                 allTogether = TRUE, colors = c(rep("orange", 1),rep("skyblue4", 1)),
                 scores = c("sum"), idName = "bams_rp")
gg_zscore <- transcriptWindow(leaders, get("cds", mode = "S4"),
                              trailers, df = df[df$condition == "rp",], outdir = NULL,
                              allTogether = TRUE, colors = c(rep("orange", 1),rep("skyblue4", 1)),
                              scores = c("zscore"), idName = "bams_rp")
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# IR
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

seqs <- startRegionString(cds[dt$txNames], tx = mrna, faFile = df, upstream = 6, downstream = 5)
hits <- nchar(seqs) == max(nchar(seqs)); summary(hits)
seqs <- seqs[hits]
kozakHeat <- kozakHeatmap(seqs = seqs, rate = dt$IR[hits],
                          start = 1, stop = max(nchar(seqs)), center = 7, type = "IR");kozakHeat
ggslackR()
ggsave(filename = p(plotFolder, "kozakHeatmap.png"), kozakHeat)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# TOP motif (Effect on SE FPKM)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
cage <- "/export/valenfs/projects/uORFome/DATA/CAGE/human/Mammary%20Epithelial%20Cell%2c%20donor1.CNhs11077.11273-116H4.hg38.nobarcode.ctss.bed.gz"

leadersCage <- reassignTSSbyCage(leaders[dt$txNames], cage)

seqs <- startRegionString(leadersCage, NULL, df, 0, 4)
rate <- dt$SE

comb <- TOP.Motif.ecdf(seqs, rate)
ggslackR(plot = comb)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# SSU processivity (loss over leader)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

library(ORFikPipeline)


df <- read.experimentl("Val20HumMer")
#filt <- dt$txNames[(dt$RNA_MRNA_FPKM > 1 & dt$SSU_LEADERS_FPKM > 55)]
filt <- dt$txNames[(dt$RNA_MRNA_FPKM > 1)]
loadRegions(txdb, c("mrna", "leaders", "cds"), names.keep = filt)

readsSSUAll <- readRDS(file = "/export/valenfs/projects/Hakon/C11/merged_SSU_bam_overlapping.rds")
leaders200 <- leaders[widthPerGroup(leaders) >= 202]
windowsOne <- startRegion(leaders200, extendLeaders(leaders200, 100), upstream = 0, downstream = 201, is.sorted = TRUE)
windowsOnetest <- removeBadTxByRegion(windowsOne, readsSSUAll, upstream = -150, downstream = 190, median_multiplier = 20, min_cutoff = 20000)
hitMapTSS <- metaWindow(readsSSUAll, windowsOnetest, scoring = "sumPos", withFrames = FALSE,
                        fraction = "SSU", feature = "TSS", zeroPosition = 0)

windowCoveragePlot(hitMapTSS, scoring = "sum", type = "relation to transcript start (nt)",
                   output = paste0(plotFolder, "leaderLoss_TSS.pdf"))

cds200 <- cds[names(leaders200)]
cds200 <- cds200[widthPerGroup(cds200, F) > 201]
windowsOne <- startRegion(cds200, mrna, upstream = 200, downstream = 201, is.sorted = TRUE)
windowsOnetest <- removeBadTxByRegion(windowsOne, readsSSUAll, median_multiplier = 20, min_cutoff = 10000)
hitMapTIS <- metaWindow(readsSSUAll, windowsOnetest, scoring = "sumPos", withFrames = FALSE,
                        fraction = "SSU", feature = "TIS", zeroPosition = 200)
windowCoveragePlot(hitMapTIS, scoring = "sum", type = "relation to start codon (nt)",
                   output = paste0(plotFolder, "leaderLoss_TIS.pdf"))

gg <- windowCoveragePlot(rbindlist(list(hitMapTSS, hitMapTIS)), scoring = "sum", type = "relation to start codon (nt)", setMinToZero = TRUE)
#gg <- gg + xlab("") + theme_classic() + theme(panel.spacing = unit(2, "lines")) + theme(strip.background = element_blank(), strip.text.x = element_blank(), title = element_blank(), legend.position = "none", strip.text = element_blank()) + ylab("Sum") +  facet_grid(fraction ~ feature, scales = "free_x",space = "free")
ggsave(gg, filename = paste0(plotFolder, "leaderLoss_TISandTSS.pdf"), dpi = 300, width = 300, height = 75, units = "mm")
ggslackR(plot = gg, width = 300, height = 75)

gg <- windowCoveragePlot(rbindlist(list(hitMapTSS, hitMapTIS)), scoring = "transcriptNormalized", type = "relation to start codon (nt)", setMinToZero = TRUE)
#gg <- gg + xlab("") + theme_classic() + theme(panel.spacing = unit(2, "lines")) + theme(strip.background = element_blank(), strip.text.x = element_blank(), title = element_blank(), legend.position = "none", strip.text = element_blank()) + ylab("Sum") +  facet_grid(fraction ~ feature, scales = "free_x",space = "free")
ggsave(gg, filename = paste0(plotFolder, "leaderLoss_TISandTSS_zscore.pdf"), dpi = 300, width = 300, height = 75, units = "mm")
ggslackR(plot = gg, width = 300, height = 75)


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Correlation plot (SE vs TE)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
dtt <- dt[TE > 0 & SE > 0, .(TE, SE)]
gg <- ggplot(data = dtt, aes(x = log2(SE), y = log2(TE), alpha = 0.5)) +
  geom_point(color = "blue2")
ggslackR(gg)
gg <- GGally::ggpairs(log2(dtt))
ggslackR(gg)
gg <- GGally::ggpairs(dtt)
ggslackR(gg)

