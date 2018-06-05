
allOK <- function(results) {
  is.error <- function(x) inherits(x, "try-error")
  ok <- !sapply(results, is.error)
  ok
}