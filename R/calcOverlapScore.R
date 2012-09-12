calcOverlapScore <-
function(x, y, from, to, bins, shift=FALSE, enrichment=FALSE) {
  if(enrichment) {
    m <- cbind(x,y)
    all <- sum(apply(m,1,function(x) abs(max(0,x)) + abs(min(0,x)) ))
    same.sign <- abs(m[sign(m[,1])==sign(m[,2]),])
    return(sum(apply(same.sign,1,min))/all)
  }
  if(shift) y <- y - mean(y, na.rm=TRUE) + mean(x, na.rm=TRUE)
  hx <- hist(x, seq(min(c(x,y), na.rm=TRUE), max(c(x,y), na.rm=TRUE), len=bins), plot=FALSE)
  hy <- hist(y, seq(min(c(x,y), na.rm=TRUE), max(c(x,y), na.rm=TRUE), len=bins), plot=FALSE)
  #if(length(hx$density) != length(hy$density)) browser()
  sum(apply(cbind(hx$density, hy$density),1,min))/sum(apply(cbind(hx$density, hy$density),1,max))
}
