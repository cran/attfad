plotEnrichment <-
function(x, y, x.name, y.name, xlab, bins, title, from, to, skip=2) {

  plot(seq(from, to, len=ceiling(length(x)/skip)), x[seq(1,length(x),skip)], t='l', col='darkgreen',lwd=2, main=title, xlab=xlab, ylab='Enrichment', ylim=c(min(c(x,y)), max(c(x,y))))
  lines(seq(from, to, len=ceiling(length(y)/skip)), y[seq(1,length(y),skip)], col='purple',lwd=2)
  abline(0,0,lty=2)

  hsc <- calcOverlapScore(x,y,enrichment=TRUE)
  csc1 <- calcCorScore(x, y, shift=FALSE, enrichment=TRUE)
  csc2 <- calcCorScore(x, y, shift=TRUE, enrichment=TRUE)

  par(mar=c(0, 0, 0, 0))
  plot.new()
  legend('center','groups', c(x.name, y.name), col=c('darkgreen', 'purple'), lty=c(1,1), ncol=2, bty='n', title=paste('Scalar score: ',round(csc2, 3),' / ', round(csc1, 3), sep='') )
  par(mar=c(5, 4, 4, 2) + 0.1)
}
