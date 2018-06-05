### Assign a group (factor)

file <- "tlr7_reviewed2.txt"

tlr <- fread(file)
tlr[Ignore == 1, Group := ""]
tlr[Ignore == 0 & isCTRL == 0, Group := LETTERS[factor(paste(BSMDCD, BSM, Node, NodeFunction))], by = .(GSE, BioSampName, XP)]
tlr[Ignore == 0 & isCTRL == 1, Group := letters[factor(paste(BSMDCD, BSM, Node, NodeFunction))], by = .(GSE, BioSampName, XP)]

write.table(tlr, "tlr7_reviewed2.txt", sep = "\t", row.names = F, quote = F)

tlr <- tlr[Ignore != 1]

key <- unique(tlr[, .(GSE, Group, xpType, XP, Node, NodeFunction, BSM, BSMDCD, BSMDCD2, BioSampName)])
write.table(key, "tlr7_key.txt", sep = "\t", row.names = F, quote = F, na = "")

fpaste <- function(A, a) {
  A <- A[names(A) != "BSMDCD2"]
  a <- a[names(a) != "BSMDCD2"]
  constant <- which(A == a)
  var <- which(A != a)
  s <- paste(paste(A[var], collapse = " "), "vs", paste(a[var], collapse = " "), "|", paste(A[constant], collapse = " "))  
  s <- gsub("\\s+", " ", trimws(s))
}

possibleContrasts <- function(dt) {
  if(nrow(dt) == 1) {
    return(c(Contrast = "", Formula = "", dt[, .(Node, NodeFunction, BSM, BSMDCD, BSMDCD2)]))
  } else if(nrow(dt) == 2) {
    group <- sort(dt$Group, decreasing = T)
    formula <- paste(group, collapse = "-")
    A <- dt[Group == group[1], .(Node, NodeFunction, BSM, BSMDCD, BSMDCD2)]
    a <- dt[Group == group[2], .(Node, NodeFunction, BSM, BSMDCD, BSMDCD2)]
    contrast <- fpaste(unlist(A), unlist(a))
    return(c(Contrast = contrast, Formula = formula, as.list(A)))
  } else {
    Node <- unique(unlist(strsplit(dt$Node, ";")))
    Node <- Node[Node != "(WT)"]
    BSM <- unique(unlist(strsplit(dt$BSM, ";")))
    BSM <- BSM[BSM != "(Veh)"]
    DCD <- unique(unlist(strsplit(dt$BSMDCD, ";")))
    allF <- c(Node, BSM, DCD)
    cX <- lapply(c(Node, BSM), function(x) allF == x) # there will always be either a Node or BSM, so not necessary to check cX 
    cT <- lapply(DCD, function(x) allF == x) # sometimes there isn't cT -> list()
    x0 <- x1 <- x2 <- list()
    if(length(Node)) x0 <- setNames(lapply(Node, function(f) grepl(f, dt$Node)), Node)
    if(length(BSM)) x1 <- setNames(lapply(BSM, function(f) grepl(f, dt$BSM)), BSM)
    if(length(DCD)) {
      x2 <- setNames(lapply(DCD, function(f) grepl(f, dt$BSMDCD)), DCD)
      cT <- suppressWarnings(mapply(function(v1, v2) setNames(v1 | v2, allF[v1]), cX, cT, SIMPLIFY = F)) # Expect arguments to be recycled
    }
    d <- as.data.frame(c(x0, x1, x2), optional = T)
    # Calculate the results of different contrast combinations
    m <- as.matrix(d)
    pairs <- combn(nrow(m), 2)
    pairs <- cbind(pairs, pairs[2:1, ])
    res <- Map(function(i, j) m[i, ] - m[j, ], pairs[1, ], pairs[2, ])
    # Find which contrast combinations isolate the desired variable 
    valid <- lapply(c(cX, cT), function(x) as.matrix(pairs[, sapply(res, function(y) all(y == x))]))
    npairs <- sapply(valid, ncol)
    effectVar <- unlist(mapply(rep, c(Node, BSM, names(cT)), npairs))
    effectNode <- effectBSM <- effectVar
    effectNode[!effectVar %in% Node] <- "" 
    effectBSM[!effectVar %in% BSM] <- ""
    valid <- do.call(cbind, valid)
    A <- dt[valid[1, ], .(Node, NodeFunction, BSM, BSMDCD, BSMDCD2)]
    a <- dt[valid[2, ], .(Node, NodeFunction, BSM, BSMDCD, BSMDCD2)]
    formula <- mapply(function(g1, g2) paste(dt[g1, "Group"], dt[g2, "Group"], sep = "-"), valid[1, ], valid[2, ])
    contrast <- mapply(function(g1, g2) fpaste(unlist(g1), unlist(g2)), split(A, 1:nrow(A)), split(a, 1:nrow(a)))
    return(list(Contrast = contrast, Formula = formula, Node = effectNode, 
                NodeFunction = A$NodeFunction, BSM = effectBSM, BSMDCD = A$BSMDCD, BSMDCD2 = A$BSMDCD2))
  }
}

contrasts <- key[, possibleContrasts(.SD), by = .(GSE, BioSampName, xpType, XP)]
write.table(contrasts, "contrasts_2.txt", sep = "\t", row.names = F, quote = F)
# manual: 5, 15*, 19, 22

# List method
dtlist <- split(key, by = c("GSE", "BioSampName", "xpType", "XP"))
results <- lapply(dtlist, function(x) try(possibleContrasts(x)))
bug <- which(sapply(results, class) == "try-error")
