# Gets some required meta data for datasets such as Pubmed ID, Contributors, last update date, etc.
# Also uses GEOmetadb to look up the appropriate Bioconductor re-annotation packages, i.e. at
# https://www.bioconductor.org/packages/3.3/data/annotation/

# library(GEOmetadb)
# library(data.table)

datasetsAnnotation <- function(con, GSELIST) {
  rs <- lapply(GSELIST, function(x) 
    { dbGetQuery(con, paste0("select gse,type,title,contributor,pubmed_id,submission_date,last_update_date from gse where gse='", x, "'")) } )
  rs <- rbindlist(rs)
  gpl <- lapply(GSELIST, function(x) { dbGetQuery(con, paste0("select gse,gpl from gse_gpl where gse='", x, "'")) } )
  gpl <- rbindlist(gpl)
  species <- lapply(gpl$gpl, function(x) { dbGetQuery(con, paste0("select gpl,organism from gpl where gpl='", x, "'")) } )
  species <- rbindlist(species)
  gpl[, species := species$organism]
  rs <- merge(rs, gpl, by = "gse")
  names(rs) <- c("GSE", "DatasetType", "Title", "DatasetContrib", "PMID", "SubmissionDate", "LastUpdateDate", "GPL", "Species") 
  # Query for BioC packages linked to given GPLs
  packages <- sapply(gpl$gpl, 
                     function(gpl) {
                       package <- dbGetQuery(con, paste0("select gpl,bioc_package from gpl where gpl='", gpl, "'"))$bioc_package
                       if(!is.na(package)) package <- paste0(package, ".db") else package <- NA
                     })
  packages <- unlist(packages)
  rs$Reannotation <- packages
  return(rs)
}

writeDatasetName <- function(contrast, info, controlterms = c("", "(WT)", "(Veh)", NA)) {
  Analysis <- ifelse(grepl("TS", contrast$xpType[1]), "Time course analysis ", "Analysis ")
  Node <- contrast$Node[!contrast$Node %in% controlterms]
  if(length(Node)) {
    Node <- paste(Node, collapse = ", ")
    Node <- paste0(Node, "-dependent ")
  }
  BSM <- contrast$BSM[!contrast$BSM %in% controlterms]
  if(length(BSM)) {
    BSM <- paste(BSM, collapse = ", ")  
    BSM <- paste0(BSM, "-regulated ")
  }
  and <- ifelse(Node != "" & BSM != "", "and ", "")
  species <- ifelse(info$Species[1] == "Homo sapiens", "human ", "mouse ")
  biosample <- paste(unique(contrast$BioSampName), collapse = " and ")
  Name <- paste0(Analysis, "of the ", Node, and, BSM, "transcriptome in ", species, biosample, ".")
  Name
}

writeDatasetDesc <- function(contrast, info) {
  biosample <- paste(unique(contrast$BioSampName), collapse = " and ")
  hasDoseInfo <- contrast$BSMDCD2 != ""
  BSM <- paste(contrast$BSMDCD2[hasDoseInfo], contrast$BSM[hasDoseInfo])
  BSM <- paste(unique(BSM), collapse = ", ")
  time <- paste(unique(contrast$BSMDCD[hasDoseInfo], collapse = ", "))
  Desc <- paste0(biosample, " were treated with ", BSM, " for ", time, ".")
  Desc
}

# Installs packages required and loads them
getBioCPackages <- function(packages) {
  for(p in packages) {
    if(!require(p)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite(p)
    }
  library(p, character.only = TRUE)
  }
}
  
geneAnnotate <- function(fit, package, gpl) {
    data <- read.table(paste0(results, fit, "_fit.txt"), sep = "\t", header = T, check.names = F)
    probes <- row.names(data)
    if(!is.na(package)) {
      data <- geneAnnotateBioC(probes, package)
    } else {
      data <- geneAnnotateGPL(fit, gpl)
    }
    return(data)
}

geneAnnotateBioC <- function(probes, package) {
  genes <- select(get(package), keys = probes, columns = "SYMBOL", keytype = "PROBEID")
  genes <- split(genes$SYMBOL, factor(genes$PROBEID))
  symbol <- sapply(genes, function(x) x[1])
  data$genes <- symbol
  return(data)
}

# Annotation with GPL requires more supervision
geneAnnotateGPL <- function(probes, gpl) {
  gpl <- getGEO(filename = paste0("GEOtemp/", gpl, ".soft"))
  cols <- colnames(Table(gpl))
  geneCol <- grep("(GENE)?_?SYMBOL", cols, ignore.case = T, val = T)
  if(!length(geneCol)) geneCol <- grep("^GENE$", ignore.case = T, cols, val = T)
  gpl <- Table(gpl)[, c("ID", geneCol)]
  genes <- gpl[[geneCol]][match(probes, gpl$ID)]
  data$genes <- genes
  return(data)
}