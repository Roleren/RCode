#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# INFO
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Pre-Analysis data of Ribo-seq made by Preeti, May 2020
library(ORFikPipeline)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Download annotation (Rnor 6.0)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Custome contaminant genomes: phix, ncRNA, (rRNA and tRNA download manually)
annotations <- getGenomeAndFasta(organism = "Rattus norvegicus",
                           output.dir = "/export/valenfs/data/references/Rnor_6.0_rat",
                           phix = TRUE, ncRNA = "rat",
                           tRNA = "/export/valenfs/data/references/Rnor_6.0_rat/rn6-mature-tRNAs.fa",
                           rRNA = "/export/valenfs/data/references/rrna/SILVA_119_bothSURef.fasta")

index <- STAR.index(annotations, wait = FALSE)

### RNA-seq alignment
# Libraries downloaded from uninett filesender (fix your files here and put them somewhere):
input.dir.rna <- "/export/valenfs/data/raw_data/RNA-Seq/Preeti_Jain_2020_Rattus_norvegicus"
output.dir.rna <- "/export/valenfs/data/processed_data/RNA-seq/Preeti_Jain_2020_Rattus_norvegicus/"
STAR.align.folder(input.dir.rna, output.dir.rna, index,
                  paired.end = "yes",
                  steps = "tr-ph-rR-nc-ge",
                  adapter.sequence = "TGGAATTCTCGGGTGCCAAGG",
                  max.cpus = 90, trim.front = 3, min.length = 20)

### Ribo-seq alignment
input.dir.rfp <- "/export/valenfs/data/raw_data/Ribo-Seq/Preeti_Jain_2020_Rattus_norvegicus/"
output.dir.rfp <- "/export/valenfs/data/processed_data/Ribo-seq/Preeti_Jain_2020_Rattus_norvegicus/"
STAR.align.folder(input.dir.rfp, output.dir.rfp, index,
                  steps = "tr-ph-rR-nc-ge",
                  adapter.sequence = "TGGAATTCTCGGGTGCCAAGG",
                  max.cpus = 90, trim.front = 0, min.length = 20)


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Create ORFik experiments of libraries
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
txdb_file <- "/export/valenfs/data/references/Rnor_6.0_rat/Rattus_norvegicus.Rnor_6.0.100_ensembl.gtf.db"
fa <- "/export/valenfs/data/references/Rnor_6.0_rat/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa"
create.experimentl(exper = "Preeti_RFP",
                   dir = "/export/valenfs/data/processed_data/Ribo-seq/Preeti_Jain_2020_Rattus_norvegicus/aligned",
                   txdb = txdb_file, fa = fa, organism = "Rattus norvegicus")


create.experimentl(exper = "Preeti_RNA",
                   dir = "/export/valenfs/data/processed_data/RNA-seq/Preeti_Jain_2020_Rattus_norvegicus/aligned",
                   txdb = txdb_file, fa = fa, organism = "Rattus norvegicus", pairedEndBam = TRUE)
# Fix columns manually in Libre office if needed

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# QC report (RNA-seq & Ribo-seq)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
df.rfp <- read.experimentl("Preeti_RFP")
ORFik::QCreport(df.rfp)

df.rna <- read.experimentl("Preeti_RNA")
ORFik::QCreport(df.rna)
ORFik:::transcriptWindow1(df = df.rna, idName = "fullReads", outdir = paste0(dirname(df.rna$filepath[1]), "/QC_STATS/"))

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# P-shifting
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
shiftFootprintsByExperiment(df.rfp)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Heatmap validation of P-shifting
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
outputLibs(df.rfp, type = "pshifted")
heatMapRegion(df.rfp, region = "TIS", acceptedLengths = NULL, scores = "transcriptNormalized")



