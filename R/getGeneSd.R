getGeneSd <-
function(gene.exp, exp.group) {
  #df <- data.frame(ex=gene.exp, group=exp.group)
  #mean(ddply(df, 'group', function(x) sd(x$ex)*(nrow(x)-1))[,2], na.rm=TRUE)
  sds <- c()
  tab <- table(exp.group)
  for(ex in unique(exp.group)) {
    if(tab[ex]>1)
      sds <- append(sds, rep(var(gene.exp[exp.group==ex]),tab[ex]-1) )
  }
  sqrt(mean(sds))
}
