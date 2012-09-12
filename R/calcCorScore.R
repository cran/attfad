calcCorScore <-
function(x, y, from=NA, to=NA, bins=NA, shift=FALSE, enrichment=FALSE) {
  if(enrichment) {
      if(shift) return(cor(x,y))
      else return((x %*% y) / (sqrt(sum(x^2)) * sqrt(sum(y^2))))
  }
  else {
    if(shift) y <- y - mean(y, na.rm=TRUE) + mean(x, na.rm=TRUE)
    hx <- hist(x, seq(min(c(x,y), na.rm=TRUE), max(c(x,y), na.rm=TRUE), len=bins), plot=FALSE)$density
    hy <- hist(y, seq(min(c(x,y), na.rm=TRUE), max(c(x,y), na.rm=TRUE), len=bins), plot=FALSE)$density
    return((hx %*% hy) / (sqrt(sum(hx^2)) * sqrt(sum(hy^2))))
  }
}
