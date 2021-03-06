
#' Get the Genomic transcript format, currently using GRch38 data
getGTF = function(assignIt = T){
  library(AnnotationDbi)
  if(exists("Gtf") == F){
    print("loading human GTF GRch38")
    if(is.null(gtfdb)){
      if(file.exists(gtfName)){
        Gtf = makeTxDbFromGFF(gtfName)
      } else {
        stop("could not find gtf file, check path")
      }
    }
    Gtf = loadDb(gtfdb)
    if(assignIt)
      assign("Gtf",Gtf,envir = .GlobalEnv)
  } else if(!dbIsValid(Gtf$conn)) {
    Gtf = ORFik:::loadTxdb(Gtf$conn@dbname)
    assign("Gtf",Gtf,envir = .GlobalEnv)
  }
}

#' Get transcripts from gtf
getTx <- function(assignIt = F){
  if (exists("tx",mode = "S4") == F) {
    getGTF()
    tx <- exonsBy(Gtf, by = "tx", use.names = TRUE)
    if (assignIt) {
      assign("tx",tx, envir = .GlobalEnv)
      return(tx)
    } else {
      return(tx)
    }
  }
}

#Get the coding sequences from the gtf file
getCDS = function(assignIt = T){
  if (!assignIt) {
    if (!file.exists(p(dataFolder, "/cds.rdata"))) {
      getGTF()
      cds = cdsBy(Gtf,"tx", use.names = TRUE)
      save(cds, file = p(dataFolder, "/cds.rdata"))
      print("saving cds rdata file for quicker reload")
    } else {
      load(p(dataFolder, "/cds.rdata"), envir = .GlobalEnv)
    }
    return(cds)
  }
  if (exists("cds",mode = "S4") == F) {
    if (!file.exists(p(dataFolder, "/cds.rdata"))) {
      getGTF()
      cds = cdsBy(Gtf,"tx", use.names = TRUE)
      save(cds, file = p(dataFolder, "/cds.rdata"))
      print("saving cds rdata file for quicker reload")
    } else {
      load(p(dataFolder, "/cds.rdata"), envir = .GlobalEnv)
    }
    if (assignIt) {
      assign("cds", cds, envir = .GlobalEnv)
      return(cds)
    } else {
      return(cds)
    }
  }
}

#Get the 3' sequences from the gtf file
getThreeUTRs = function(){
  if (!exists("threeUTRs")) {
    if (!file.exists(p(dataFolder, "/threeUTRs.rdata"))) {
      getGTF()
      threeUTRs = threeUTRsByTranscript(Gtf, use.names = TRUE)
      assign("threeUTRs", threeUTRs, envir = .GlobalEnv)
      save(threeUTRs, file = p(dataFolder, "/threeUTRs.rdata"))
      print("saving 3'UTRs rdata file for quicker reload")
      return(NULL)
    } else {
      load(p(dataFolder, "/threeUTRs.rdata"), envir = .GlobalEnv)
    }
  }
}

#' Get the 5' leaders, either from gtf, cage data to reassign
#' the transcription start site(TSS), or load from existing data
#' either as .rdata or .bed (bed6)
getLeaders = function(cageName = NULL, assignLeader = T, exportUorfRegions = T) {
  if(exists("fiveUTRs") == F){
    cat("creating leader from scratch\n")

    if (file.exists(p(dataFolder,"/leader.rdata"))) {
      load(p(dataFolder,"/leader.rdata"))
    } else {
      getGTF()
      cat("loading Leader from gtf\n")
      fiveUTRs = fiveUTRsByTranscript(Gtf,use.names = T)
      save(fiveUTRs, file = p(dataFolder,"/leader.rdata"))
    }

    if (!is.null(cageName)) {
      print("Using cage.. ")

      fiveUTRs = ORFik::reassignTSSbyCage(fiveUTRs, cageName, filterValue = 3)
      if (exportUorfRegions) {
        getCDS()
        uORFSeachRegion <- ORFik:::addCdsOnLeaderEnds(fiveUTRs, cds, onlyFirstExon = F)
        uORFSeachRegion <- sortPerGroup(uORFSeachRegion)
        print("exporting new uorf regions")
        exportNamerdata = paste0(regionUORFsFolder,
                                 basename(p(cageName, ".regionUORF.rdata")))
        save(uORFSeachRegion, file = exportNamerdata)
      }

      print("exporting new leaders")
      exportNamerdataLeader = paste0(leadersFolder,
                               basename(p(cageName, ".leader.rdata")))
      save(fiveUTRs,file = exportNamerdataLeader)
      print("finished new 5' UTRs")
    }
  }
  else{
    print("fiveUTRs already exists! cancel if this is wrong!")
  }

  if(assignLeader)
    assign("fiveUTRs",fiveUTRs, envir = .GlobalEnv)

  print("finished loading leaders")
}

leaderCage <- function(with.cds = TRUE){
  if(with.cds) {
    load(p(dataBaseFolder,"/CageFiveUTRsWithCDS.rdata"))
    return(CageFiveWithCDS)
  }
  load(p(dataBaseFolder,"/CageFiveUTRs.rdata"))

  return(CageFiveUTRs)
}

#' Convenience wrapper from findMapORFs
#' Get the upstream open reading frames from the 5' leader sequences, given as GRangesList
getUnfilteredUORFs = function(uORFSeachRegion, assignRanges = TRUE, isSorted = FALSE,
                              startCodons = "ATG", groupByTx = FALSE, minimumLength = 2){
  print("making leader sequences")
  getSequencesFromFasta(uORFSeachRegion, isSorted)

  print("find uORFs")
  rangesOfuORFs = findMapORFs(grl = uORFSeachRegion, seqs = seqs, longestORF = FALSE,
                                     minimumLength = minimumLength, startCodon = startCodons,
                                     groupByTx = groupByTx)

  if(assignRanges)
    assign("rangesOfuORFs",rangesOfuORFs,envir = .GlobalEnv)
  print("finished unfiltered UORFs")
  return(rangesOfuORFs)
}

#' Get the fasta indexed file
#'
#' if assignIt is TRUE, the object is not return to local scope
#' Only assigned to globalenvir
getFasta = function(filePath = NULL, assignIt = T){

  if(exists("fa") == F){ #index files
    if (is.null(filePath)){
      fa = FaFile(faiName)
    } else {
      fa = FaFile(filePath)
    }
    if (assignIt){
      assign("fa",fa,envir = .GlobalEnv)
    } else {
      return(fa)
    }
  }
}

#' get sequences from a GRangeslist
getSequencesFromFasta = function(grl, isSorted = T){
  getFasta() #get .fai
  if(!isSorted) grl <- ORFik:::sortPerGroup(grl)
  seqs = extractTranscriptSeqs(fa, transcripts = grl)
  assign("seqs",seqs,envir = .GlobalEnv)
}

getAll <- function(include.cage = T, cdsOnFiveEnd = F){
  getFasta()

  getCDS()
  getThreeUTRs()
  getLeaders()

  #or with extension
  if (include.cage) {
    cageFiveUTRs <- leaderCage(cdsOnFiveEnd)
    assign("cageFiveUTRs", cageFiveUTRs,  envir = .GlobalEnv)
    getCageTx()
  } else {
    getTx(T)
  }
  return(NULL)
}

#' Make directoy structure for orf finding
#'
#' The main Path is ./.. relative to RCode1/ location
orfikDirs <- function(mainPath, makeDatabase = FALSE){
  setwd(mainPath)
  print(paste("main path for project will be: ", mainPath))
  resultsLoc <- resultsFolder
  if (!dir.exists(resultsLoc)) dir.create(resultsLoc)

  dir.create(p(resultsLoc,"/New_Cage_Leaders"))
  dir.create(p(resultsLoc,"/regionUORFs"))
  dir.create(p(resultsLoc,"/rangesOfUORFs"))
  dir.create(p(resultsLoc,"/fasta"))
  dir.create(p(resultsLoc,"/uorfIDs"))

  if (makeDatabase) {
    dir.create("dataBase")
    dir.create("dataBase/forests/")
    dir.create("dataBase/forests/predicateTables")
  }

  print("directories created successfully")
}

#Check if uorfRanges exist already, or if must be created.
###########Should make this more failsafe!!!!!!!!!! add possibility to give ranges!!!!!!!!
UorfRangesNotExists <- function(assignUorf = F, givenCage = NULL){
  if(!exists("rangesOfuORFs")){
    if(file.exists(getUORFRDataName(givenCage))){#!!!Will not work for single run now!!!
      if(assignUorf){
        cat("loading rangesOfuorf from folder\n",uorfFolder,"\n")
        load(getUORFRDataName(givenCage),envir = .GlobalEnv)
      }
      return(F)
    }else{return(T)}
  }
  return(F)
}

#Get the number of overlaps between the upstream open reading frames and the coding sequences
getUOrfOverlaps = function(){
  overlap1 = findOverlaps(cds,rangesOfuORFs)
  overlapCount = countOverlaps(cds,rangesOfuORFs)
  numberOfOverlaps = sum(overlapCount >= 1)
  overlapHitsIndex = overlapCount[overlapCount == 1]
  return(numberOfOverlaps)
}

getCageTx <- function() {
  if (file.exists(p(dataBaseFolder, "/cageTx.rdata"))) {
    load(p(dataBaseFolder, "/cageTx.rdata"), envir = .GlobalEnv)
  } else {
    tx <- getTx()
    cageFiveUTRs <- leaderCage()
    tx[names(cageFiveUTRs)] <- ORFik:::extendLeaders(tx, cageFiveUTRs)
    assign("tx", tx,  envir = .GlobalEnv)
    save(tx, file = p(dataBaseFolder, "/cageTx.rdata"))
  }
  return(NULL)
}
