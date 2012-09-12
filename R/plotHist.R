plotHist <-
function(x, y, x.name, y.name, xlab, bins, log, title, from, to) {

  if(log) {
    min.bin <- min(c(x[which(x!=0)], y[which(y!=0)]))
    x[which(x==0)] <- min.bin
    y[which(y==0)] <- min.bin
    x <- log10(x)
    y <- log10(y)
  }
  from <- ifelse(is.na(from),min(c(x,y),na.rm=T),from)
  to <- ifelse(is.na(to),max(c(x,y),na.rm=T),to)

  #legend('top', '(x,y)', c(paste(c('deming 25%-75%', 'diagonal'),c(round(scDm,3), round(scDi,3)),sep=': ')), col=c('red', 'blue'),lty=c(1,2))
  hx <- hist(x,seq(from,to,len=bins),plot=F)
  hy <- hist(y,seq(from,to,len=bins),plot=F)

  first.x <- max(1,which(hx$density>0)[1]-1)
  last.x <- min(length(hx$density), rev(which(hx$density>0))[1]+1)
  first.y <- max(1,which(hy$density>0)[1]-1)
  last.y <- min(length(hy$density), rev(which(hy$density>0))[1]+1)
  hx$mids <- hx$mids[first.x:last.x]
  hx$density <- hx$density[first.x:last.x]*(to-from)
  hy$mids <- hy$mids[first.y:last.y]
  hy$density <- hy$density[first.y:last.y]*(to-from)
  x.rel.freq <- hx$density
  y.rel.freq <- hy$density
  if(all(is.nan(x.rel.freq))) x.rel.freq <- rep(1e6, length(x.rel.freq))
  if(all(is.nan(y.rel.freq))) y.rel.freq <- rep(1e6, length(y.rel.freq))
  
  plot(hx$mids,x.rel.freq,t='l',col='darkgreen',axes=!log, main=title,lwd=2,xlab=xlab,ylab='density',xlim=c(from, to),ylim=c(0,max(c(x.rel.freq,y.rel.freq))))
  lines(hy$mids,y.rel.freq,col='purple',lwd=2)
  

  if(log) {
    axis(2)
    from <- min(c(hx$mids,hy$mids))
    to <- max(c(hx$mids,hy$mids))
    at <- seq(from, to, len=5)
    pos <- log10(c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1))
    at <- pos[pos>(from-0.5) & pos<to]
    axis(1, at=at, labels=10^at)
    box()
  }
  
  hsc1 <- calcOverlapScore(x, y, from, to, bins, shift=FALSE, enrichment=FALSE)
  hsc2 <- calcOverlapScore(x, y, from, to, bins, shift=TRUE, enrichment=FALSE)
  csc1 <- calcCorScore(x, y, from, to, bins, shift=FALSE, enrichment=FALSE)
  csc2 <- calcCorScore(x, y, from, to, bins, shift=TRUE, enrichment=FALSE)
  
  par(mar=c(0, 0, 0, 0))
  plot.new()
  legend('center','groups', c(x.name, y.name), col=c('darkgreen', 'purple'), lty=c(1,1), ncol=2, bty='n', title=paste('Overlap score: ',round(hsc1, 3),' / ', round(hsc2, 3),'\nScalar score: ',round(csc1, 3),' / ', round(csc2, 3), sep='') )
  par(mar=c(5, 4, 4, 2) + 0.1)
}
