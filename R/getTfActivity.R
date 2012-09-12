getTfActivity <-
function(nmat, targets, cutoff) {
  p.vals <- c()
  for(i in 1:ncol(nmat)) {
    if(sum(!is.na(nmat[targets,i])) > 1) {
      p.vals <- append(p.vals, wilcox.test(nmat[targets,i], nmat[! rownames(nmat) %in% as.character(targets),i])$p.value)
    }
  }

  ### multiple testing correction
  p.vals <- p.adjust(p.vals, method='BH')
  if(length(p.vals)>0)
    return(mean(p.vals < cutoff))
  NA
}
