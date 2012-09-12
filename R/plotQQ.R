plotQQ <-
function(x, y, xlab, ylab, title, maxpoints=1e6) {
  
  #if(!is.installed('MethComp')) install.packages('MethComp')
  #require('MethComp')

  qq <- qqplot(sample(x,min(length(x),maxpoints)), sample(y,min(length(y),maxpoints)), plot.it=T, cex=0.5, main=title, xlab=xlab, ylab=ylab, axes=T)
  abline(0,1,col='blue',lty=2)
  points(quantile(qq$x,c(0.25,0.5,0.75)), quantile(qq$y,c(0.25,0.5,0.75)), pch=c('(','x',')'), col='green', cex=2)

  dm2 <- Deming(qq$x[as.integer(0.25*length(qq$x)):as.integer(0.75*length(qq$x))],qq$y[as.integer(0.25*length(qq$y)):as.integer(0.75*length(qq$y))])
  if(is.finite(dm2[1]) & is.finite(dm2[2])) {
    abline(dm2[1:2],col='red',lty=2)
  }
}
