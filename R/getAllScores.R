getAllScores <-
function(results.1, results.2) {
  defaults <- getDefaults()
  scores <- matrix(NA, nrow=0, ncol=4)
  for(i in names(results.1)) {
    if(i %in% names(results.2) & i != 'expression.data') {
      scores <- rbind(scores, i=getScores(results.1[[i]], results.2[[i]], defaults[i,"from"], defaults[i,"to"], defaults[i,"bins"], defaults[i,"enrichment"]))
      rownames(scores)[nrow(scores)] <- i
    }
  }
  scores
}
