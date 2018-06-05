# Gets some required meta data for datasets such as Pubmed ID, Contributors, last update date, etc.
# Also uses GEOmetadb to look up the appropriate Bioconductor re-annotation packages, i.e. at
# https://www.bioconductor.org/packages/3.3/data/annotation/

library(GEOmetadb)
library(data.table)

setwd("~/NURSA/GEOmetric")
results <- "TLR/Results1/"
outfile <- ""
  
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")

gselist <- unique(gsub("_.*", "", list.files(results)))

rs <- lapply(gselist, function(x) 
  { dbGetQuery(con, paste0("select gse,type,title,contributor,pubmed_id,submission_date,last_update_date from gse where gse='", x, "'")) } )
rs <- rbindlist(rs)
gpl <- lapply(gselist, function(x) { dbGetQuery(con, paste0("select gse,gpl from gse_gpl where gse='", x, "'")) } )
gpl <- rbindlist(gpl)
species <- lapply(gpl$gpl, function(x) { dbGetQuery(con, paste0("select gpl,organism from gpl where gpl='", x, "'")) } )
species <- rbindlist(species)
gpl[, species := species$organism]
rs <- merge(rs, gpl, by = "gse")

packages <- sapply(gpl$gpl, 
            function(gpl) {
                package <- dbGetQuery(con, paste0("select gpl,bioc_package from gpl where gpl='", gpl, "'"))$bioc_package
                if(!is.na(package)) package <- paste0(package, ".db") else package <- NA
            })
packages <- unlist(packages)
rs$Reannotation <- packages
# View the packages required and download them.
for (p in packages) {
  source("https://bioconductor.org/biocLite.R")
  biocLite(package)
  library(package, character.only = TRUE)
}

write.table(rs, file = "Datasets_Experiments.txt", sep = "\t", row.names = F)
  
Annotate <- function(fit, package, gpl) {
    data <- read.table(paste0(results, fit, "_fit.txt"), sep = "\t", header = T, check.names = F)
    probes <- row.names(data)
    if(!is.na(package)) {
      genes <- select(get(package), keys = probes, columns = "SYMBOL", keytype = "PROBEID")
      genes <- split(genes$SYMBOL, factor(genes$PROBEID))
      symbol <- sapply(genes, function(x) x[1])
      data$genes <- symbol
    } else {
      gpl <- getGEO(filename = paste0("GEOtemp/", gpl, ".soft"))
      cols <- colnames(Table(gpl))
      geneCol <- grep("(GENE)?_?SYMBOL", cols, ignore.case = T, val = T)
      if(!length(geneCol)) geneCol <- grep("^GENE$", ignore.case = T, cols, val = T)
      gpl <- Table(gpl)[, c("ID", geneCol)]
      genes <- gpl[[geneCol]][match(probes, gpl$ID)]
      data$genes <- genes
    }
    write.table(data, file = paste0(results, fit, "_annotated.txt"))
} 
         
for (i in nrow(rs)) Annotate(rs$gse[i], package = rs$Reannotation[i], gpl = rs$gpl[i])

