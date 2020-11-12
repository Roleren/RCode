#' Read standard experiment on server
#' @param relPath name of experiment, run list.experiments() for candidates
#' @inheritParams list.experiments
#' @param expInVarName keep experiment name in variable names? defualt (FALSE)
#' @importFrom ORFik read.experiment
#' @export
read.experimentl <- function(relPath, dir = "/export/valenfs/data/processed_data/experiment_tables_for_R/",
                             expInVarName = FALSE) {
  if (tools::file_ext(relPath) == "") {
    df <- read.experiment(p(dir, relPath, ".csv"))
  } else if (tools::file_ext(relPath) == "csv"){
    df <- read.experiment(p(dir, relPath))
  } else stop("invalid file type of file")
  df@expInVarName <- expInVarName
  return(df)
}

#' Create standard zebrafish experiment on server
#' Using Z10
#' @inheritParams ORFik::create.experiment
#' @param txdb ("/export/valenfs/projects/uORFome/Annotations/Zebrafish/zebrafish_GRCh10_81.gtf.db")
#' @param fa a fasta index path ("/export/valenfs/projects/uORFome/Annotations/Zebrafish/zebrafish_GRCh10_81.gtf.db")
#' @param saveDir directory for ORFik experiments: default ("/export/valenfs/data/processed_data/experiment_tables_for_R/")
#' @param organism Danio rerio
#' @importFrom ORFik create.experiment
#' @export
create.experimentz <- function(dir, exper,
                               saveDir = "/export/valenfs/data/processed_data/experiment_tables_for_R/", types = c("bam", "bed", "wig"),
                               txdb = "/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.81_chr.gtf.db",
                               fa = "/export/valenfs/data/references/Zv10_zebrafish/Danio_rerio.GRCz10.fa",
                               organism = "Danio rerio",
                               viewTemplate = FALSE, pairedEndBam = FALSE) {
  create.experiment(dir, exper, saveDir, txdb, fa, organism, pairedEndBam, viewTemplate, types)
}




#' Create standard human experiment on server
#' Using hg38
#' @inheritParams ORFik::create.experiment
#' @param txdb ("/export/valenfs/data/references/Homo_sapiens_GRCh38_110/Homo_sapiens.GRCh38.101_ensembl.gtf.db")
#' @param fa a fasta index path ("/export/valenfs/data/references/Homo_sapiens_GRCh38_110/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
#' @param saveDir directory for ORFik experiments: default ("/export/valenfs/data/processed_data/experiment_tables_for_R/")
#' @param organism "Homo sapiens"
#' @importFrom ORFik create.experiment
#' @export
create.experimenth <- function(dir, exper,
                               saveDir = "/export/valenfs/data/processed_data/experiment_tables_for_R/", types = c("bam", "bed", "wig"),
                               txdb = "/export/valenfs/data/references/Homo_sapiens_GRCh38_110/Homo_sapiens.GRCh38.101_ensembl.gtf.db",
                               fa = "/export/valenfs/data/references/Homo_sapiens_GRCh38_110/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
                               organism = "Homo sapiens",
                               viewTemplate = FALSE, pairedEndBam = FALSE) {
  create.experiment(dir, exper, saveDir, txdb, fa, organism, pairedEndBam, viewTemplate, types)
}


# For legacy reasons

#' @inherit create.experimentz
#' @export
create.experimentl <- create.experimentz

#' Save template experiment to disc
#' @param temp a temperary ORFik experiment template
#' @param saveDir directory for ORFik experiments: default ("/export/valenfs/data/processed_data/experiment_tables_for_R/")
#' @importFrom ORFik save.experiment
#' @export
save.experimentl <- function(temp, saveDir = "/export/valenfs/data/processed_data/experiment_tables_for_R/") {
  save.experiment(temp, file = p(saveDir, temp[1,2]))
}

#' List current experiment available
#'
#' Will only search .csv extension, also exclude any experiment with the word template.
#' @param dir directory for ORFik experiments: default ("/export/valenfs/data/processed_data/experiment_tables_for_R/")
#' @param pattern allowed patterns in experiment file name: default ("*", all experiments)
#' @param libtypeExclusive search for experiments with exclusivly this libtype, default (NULL, all)
#' @importFrom ORFik list.experiments
#' @export
list.experiments <- function(dir = "/export/valenfs/data/processed_data/experiment_tables_for_R/",
                             pattern = "*", libtypeExclusive = NULL) {
  return(ORFik::list.experiments(dir, pattern, libtypeExclusive))
}
