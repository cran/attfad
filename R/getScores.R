getScores <-
function(x, y, from, to, bins, enrichment) {
  x <- c(x)
  y <- c(y)
  if(length(na.omit(x)) < 5 | length(na.omit(y)) < 5)
    return(rep(NA, 4))
  osc1 <- calcOverlapScore(x, y, from, to, bins, enrichment, shift=FALSE)
  osc2 <- calcOverlapScore(x, y, from, to, bins, enrichment, shift=TRUE)
  csc1 <- calcCorScore(x, y, from, to, bins, enrichment, shift=FALSE)
  csc2 <- calcCorScore(x, y, from, to, bins, enrichment, shift=TRUE)  
  c(overlap=osc1, overlap.shifted=osc2, correlation=csc1, correlation.shifted=csc2)
}
