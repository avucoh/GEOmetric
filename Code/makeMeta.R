# Simply takes a list of GSE ids and pulls the metadata. 

library(data.table)
library(limma)
library(GEOquery)
library(magrittr)

file <- "TLR/Start/tlr7.txt"

tlr <- fread(file, sep = "\t")
gses <- tlr[select == 1, gse]

getSamples <- function(gse) {
  eset <- getGEO(GEO = gse, destdir = "GEOtemp")[[1]]
  pData(eset)
}

results <- lapply(gses, function(x) try(getSamples(x)))
names(results) <- gses

is.error <- function(x) inherits(x, "try-error")
ok <- !vapply(results, is.error, logical(1))
table(ok)

# select only a subset of p data
subP <- function(data) {
  cols <- grep("title|geo_accession|source_name|description|platform|characteristics", names(data))
  data <- data[, cols]
  data <- lapply(data, function(x) gsub("\n", " ", x))
  as.data.table(data)
} 

GSExGSM <- unlist(mapply(rep, gses, sapply(results, nrow))) 

export <- lapply(results, subP)
export <- rbindlist(export, fill = T)


export[, c("Group", "Node", "NodeFunction", "BSM", "BSMDCD", "BSMDCD2", "BSMMapped", "BioSampName", "BioSampID", "BioSampExt", "Comment") := ""]
export[, c("Ignore", "Batch", "isControl") := 0]
export[, GSE := GSExGSM]

# Uses some simple rules to assign values; these need to be confirmed manually after exporting


inferBSMDCD <- function(s) {
  m <- regexpr("[0-9]+ ?(d|h|D|H)|(D|d)(ay) ?[0-9]+", s)
  s[m == -1] <- ""
  s[m != -1] <- regmatches(s, m)
  s <- tolower(gsub(" ", "", s))
}


allbio <- fread("allbiosamples.txt")
setnames(allbio, c("ID", "System", "Organ", "TissueCell1", "TissueCell2", "TissueCell3", "TissueCell4", "TissueCell5"))
for(col in names(allbio)) set(allbio, j = col, value = trimws(allbio[[col]]))
biosamps <- unique(allbio[, c(TissueCell1, TissueCell2, TissueCell3)]) 

inferBiosamp <- function(s) {
  m <- lapply(s, function(x) which(sapply(biosamps, function(y) grep(y, x, fixed = T)) == 1))
  m <- sapply(m, function(x) trimws(paste(biosamps[x], collapse = " ")))
  m
}

export[, isControl := as.numeric(1:.N %in% grep("wt|wildtype|control|ctrl|healthy|media|medium|veh|unstim|untreated|mock", title, ignore.case = T))]
export[, BioSampName := inferBiosamp(paste(title, source_name_ch1, source_name_ch2))]
export[, BSMDCD := inferBSMDCD(title)]
neworder <- c("GSE", "geo_accession", "Comment", "Group", "Ignore", "Batch", "isControl", 
              "Node", "NodeFunction", "BSM", "BSMDCD", "BSMDCD2", "BSMMapped", "BioSampName",          
              "BioSampID", "BioSampExt", "title", grep("^source", names(export), val = T),
              grep("^char", names(export), val = T), "platform_id", grep("^desc", names(export), val = T))
setcolorder(export, neworder)

write.table(export, "tlr7_metadata.txt", sep = "\t", row.names = F, quote = F)




