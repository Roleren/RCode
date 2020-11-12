#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Prepare data for C11 analysis

library(ORFikPipeline)
df <- read.experimentl("Val20HumMer", expInVarName = TRUE)
outputLibs(df, type = "bedo")
dfr <- read.experimentl("Val20Hum")
tables <- countTable(dfr, region = "mrna", type = "fpkm")
tables <- tables[, "RNA_"]; summary(tables)
mrnaOri <- loadRegion(df, "mrna")
rnaFilt <- (tables$RNA_ > 1)
mrna <- mrnaOri[rnaFilt]

loadRegions(df, parts = c("leaders", "cds"))

if (!all(names(leaders) %in% names(cds))) stop("Not all leaders present!")


dt <- data.table()
tablesLeaders <- countTable(df, "leaders", type = "fpkm")
tablesCDS <- countTable(df, "cds", type = "fpkm")
tablesCDS <- tablesCDS[(names(cds) %in% names(leaders)) & (names(cds) %in% names(mrna)),]
tablesLeaders <- tablesLeaders[names(leaders) %in% names(mrna),]
if (nrow(tablesLeaders) != nrow(tablesCDS)) stop("not matching rows!")

dt$txNames <- names(leaders)[names(leaders) %in% names(mrna)]
dt$LSU_CDS_FPKM <- tablesCDS$LSU_rp + tablesCDS$LSU_WT
dt$SSU_LEADERS_FPKM <- tablesLeaders$LSU_rp + tablesLeaders$LSU_WT
dt$RNA_MRNA_FPKM <- tables$RNA_[names(mrnaOri) %in% dt$txNames]
dtOri <- copy(dt)
dt <- dt[SSU_LEADERS_FPKM > 0,]
dt$IR <- dt$LSU_CDS_FPKM / dt$SSU_LEADERS_FPKM
dt$SE <- dt$SSU_LEADERS_FPKM / dt$RNA_MRNA_FPKM
dt$TE <- dt$LSU_CDS_FPKM / dt$RNA_MRNA_FPKM

saveRDS(dt, "/export/valenfs/projects/Hakon/C11/matrix.rds")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# CAGE
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

cage <- "/export/valenfs/projects/uORFome/DATA/CAGE/human/Mammary%20Epithelial%20Cell%2c%20donor1.CNhs11077.11273-116H4.hg38.nobarcode.ctss.bed.gz"
txdb <- loadTxdb(df)
txdb <- reassignTxDbByCage(txdb, cage)
saveDb(txdb, "/export/valenfs/projects/Hakon/C11/GRCH38.CAGE.Mammary.gtf.db")


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Merged SSU BAM
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
remove.experiments(df)
df <- read.experimentl("Val20HumMer")
outputLibs(df[df$libtype == "SSU",], chrStyle = txdb)
loadRegions(txdb, c("mrna"), names.keep = dt$txNames)
SSU_rp <- SSU_rp[unique(from(findOverlaps(SSU_rp, mrna, type = "within")))]
SSU_WT <- SSU_WT[unique(from(findOverlaps(SSU_WT, mrna, type = "within")))]
readsSSUAll <- c(SSU_rp, SSU_WT)
saveRDS(readsSSUAll, file = "/export/valenfs/projects/Hakon/C11/merged_SSU_bam_overlapping.rds")
