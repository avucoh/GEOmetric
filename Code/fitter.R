library(data.table)
library(GEOquery)
library(limma)
library(magrittr)

datadir <- "GEOtemp/"
resultsdir <- "Results/"
DT <- fread("reviewed_tlr7.txt")
DT <- DT[GSE != "", 1:17]
setkey(DT, GSE)
contrasts <- fread("contrasts_tlr7.txt")
contrasts <- split(contrasts, by = c("GSE", "BioSampName"))

### Module Functions
fitDesign <- function(contrast) {
  gse <- contrast$GSE[1]
  dt <- DT[gse]
  bio <- contrast$BioSampName[1]
  xpType <- contrast$xpType[1]
  path <- paste0(datadir, gse, "_series_matrix.txt.gz")
  if(!file.exists(path)) getGEO(gse, parseCharacteristics = F)
  gse <- getGEO(filename = path, parseCharacteristics = F)
  # make sure to select samples of same biosample type because gses will include multiple cell types
  select <- (dt$BioSampName == bio) & (dt$Ignore == 0) 
  eset <- exprs(gse)
  eset <- eset[, select]
  dt <- dt[select]
  batch <- factor(dt$Batch) 
  group <- factor(dt$Group)
  if(grepl("1C", xpType)) {
    if(nlevels(batch) > 1) { # Account for batch/block design by including term in formula
      design <- model.matrix(~0 + group + batch)
      colnames(design) <- gsub("group|batch", "", colnames(design))
    } else {
      design <- model.matrix(~0 + group)
      colnames(design) <- levels(group)
    }
    fit <- lmFit(eset, design)
  } else if(grepl("2CDirect|2CCommonRef", xpType)) {
    targets <- dt[, .(geo_accession, RefDye, Group)]
    targets[, Cy3 := ifelse(RefDye == "Cy3", "Ref", Group)]
    targets[, Cy5 := ifelse(RefDye == "Cy5", "Ref", Group)]
    targets <- targets[, .(geo_accession, Cy3, Cy5)]
    design <- modelMatrix(targets, ref = "Ref")
    fit <- lmFit(eset, design)
  } else if(grepl("Htseq", xpType)) {
    ## TO DO: call voom -- need to review how to process HT data first
    return(NULL)
  } else {
    return(NULL)
  }
  return(fit)
}

fitContrast <- function(fit, contrasts) {
  cont.matrix <- makeContrasts(contrasts = contrasts$Formula, levels = fit$design)
  colnames(cont.matrix) <- contrasts$Contrast
  rownames(cont.matrix) <- colnames(design)
  fit2 <- contrasts.fit(fit, cont.matrix) %>% eBayes()
  return(fit2)
}

# Wrapper for above calls
try2getFit <- function(contrast) {
  try({
    fit <- fitDesign(contrast)
    if(is.null(fit)) return(NULL)
    fit2 <- try(fitContrast(fit, contrast))
    return(fit2)
  })
}

results <- lapply(contrasts, try2getFit)
names(results) <- names(contrasts)
results <- results[!sapply(results, is.null)] # The Htseq and other xp types that will be handled later
error <- which(sapply(results, class) == "try-error") 
# No errors
# results <- results[-error] 

for(gse in names(results)) {
write.fit(results[[gse]], results = NULL, file = paste0(resultsdir, gse, "_fit.txt"), digits = 10, 
          adjust = "fdr", method = "separate", F.adjust = "none", sep="\t")
}
