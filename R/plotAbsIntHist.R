plotAbsIntHist <-
function(x, y, x.name, y.name, xlab, bins, log, title, from, to, mat1, mat2) {

  from <- ifelse(is.na(from),min(c(x,y),na.rm=T),from)
  to <- ifelse(is.na(to),max(c(x,y),na.rm=T),to)

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

  sd1 <- apply(mat1, 2, sd)
  sd2 <- apply(mat2, 2, sd)
  hsd1l <- hist(mat1[,which(abs(sd1-quantile(sd1,0.01))==min(abs(sd1-quantile(sd1,0.01))))],seq(min(from,min(mat1,na.rm=TRUE)),max(to,max(mat1,na.rm=T)),len=floor(bins/2)),plot=F)
  hsd2l <- hist(mat2[,which(abs(sd2-quantile(sd2,0.01))==min(abs(sd2-quantile(sd2,0.01))))],seq(min(from,min(mat2,na.rm=TRUE)),max(to,max(mat2,na.rm=T)),len=floor(bins/2)),plot=F)
  hsd1h <- hist(mat1[,which(abs(sd1-quantile(sd1,0.99))==min(abs(sd1-quantile(sd1,0.99))))],seq(min(from,min(mat1,na.rm=TRUE)),max(to,max(mat1,na.rm=T)),len=floor(bins/2)),plot=F)
  hsd2h <- hist(mat2[,which(abs(sd2-quantile(sd2,0.99))==min(abs(sd2-quantile(sd2,0.99))))],seq(min(from,min(mat2,na.rm=TRUE)),max(to,max(mat2,na.rm=T)),len=floor(bins/2)),plot=F)
  hsd1ldens <- hsd1l$density[ceiling(first.x/2):floor(last.x/2)]*(to-from)
  hsd2ldens <- hsd2l$density[ceiling(first.y/2):floor(last.y/2)]*(to-from)
  hsd1hdens <- hsd1h$density[ceiling(first.x/2):floor(last.x/2)]*(to-from)
  hsd2hdens <- hsd2h$density[ceiling(first.y/2):floor(last.y/2)]*(to-from)

  plot(hx$mids,x.rel.freq,t='l',col='darkgreen',axes=!log, main=title,lwd=2,xlab=xlab,ylab='density',xlim=c(from, to),ylim=c(0,max(c(x.rel.freq,y.rel.freq,hsd1ldens,hsd2ldens,hsd1hdens,hsd2hdens))))
  lines(hy$mids,y.rel.freq,col='purple',lwd=2)
  
  lines(hsd1l$mids[ceiling(first.x/2):floor(last.x/2)], hsd1ldens,lty=2,lwd=1,col='gray10')
  lines(hsd2l$mids[ceiling(first.y/2):floor(last.y/2)], hsd2ldens,lty=2,lwd=1,col='gray10')
  lines(hsd1h$mids[ceiling(first.x/2):floor(last.x/2)], hsd1hdens,lty=4,lwd=1,col='gray40')
  lines(hsd2h$mids[ceiling(first.y/2):floor(last.y/2)], hsd2hdens,lty=4,lwd=1,col='gray40')

  hsc1 <- calcOverlapScore(x, y, from, to, bins, shift=FALSE, enrichment=FALSE)
  hsc2 <- calcOverlapScore(x, y, from, to, bins, shift=TRUE, enrichment=FALSE)
  csc1 <- calcCorScore(x, y, from, to, bins, shift=FALSE, enrichment=FALSE)
  csc2 <- calcCorScore(x, y, from, to, bins, shift=TRUE, enrichment=FALSE)

  par(mar=c(0, 0, 0, 0))
  plot.new()
  legend('center','groups', c(x.name, y.name), col=c('darkgreen', 'purple'), lty=c(1,1), ncol=2, bty='n', title=paste('Overlap score: ',round(hsc1, 3),' / ', round(hsc2, 3),'\nScalar score: ',round(csc1, 3),' / ', round(csc2, 3), sep='') )
  par(mar=c(5, 4, 4, 2) + 0.1)
}
