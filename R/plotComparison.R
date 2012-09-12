plotComparison <-
function(results.1, results.2, name.1, name.2, stat, type='both', extendedAbsoluteIntensity=FALSE, bins=NA) {
  
  defaults <- getDefaults()
  stat <- rownames(defaults)[pmatch(stat, rownames(defaults))]
  if(is.na(bins))
    bins = defaults[stat,"bins"]

  x <- results.1[[stat]]
  y <- results.2[[stat]]
  if(sum(!is.na(x) & is.finite(x) & is.finite(1/x)) < 5 | sum(!is.na(y) & is.finite(y) & is.finite(1/y)) < 5) {
    cat(paste('Not enough finite values to plot ',stat,'.\n',sep=''))
    plot.new()
    text(0.5,0.5, paste('Not enough finite values to plot\n',stat,' for\n',name.1,' and ',name.2,'!',sep=''))
    return()
  }
  
  if(type=='both')
    layout(matrix(1:4,2,2,byrow=FALSE), heights=c(7,1,7,1))
  else
    layout(matrix(1:2,2,1,byrow=FALSE), heights=c(7,1))

  if ((type == 'hist' | type=='both') & stat == "Absolute.intensity.distribution" & extendedAbsoluteIntensity == TRUE)
    plotAbsIntHist(x=x, y=y, x.name=name.1, y.name=name.2, title=gsub("."," ",stat,fixed=TRUE), xlab=defaults[stat,"xlab"], bins=bins, from=defaults[stat,"from"], to=defaults[stat,"to"], log=defaults[stat,"log"], mat1=results.1[['expression.data']], mat2=results.2[['expression.data']])

  else if ((type == 'hist' | type=='both') & defaults[stat,"enrichment"] == FALSE)
    plotHist(x=x, y=y, x.name=name.1, y.name=name.2, title=gsub("."," ",stat,fixed=TRUE), xlab=defaults[stat,"xlab"], bins=bins, from=defaults[stat,"from"], to=defaults[stat,"to"], log=defaults[stat,"log"])
  
  else if ((type == 'hist' | type=='both') & defaults[stat,"enrichment"] == TRUE)
    plotEnrichment(x=x, y=y, x.name=name.1, y.name=name.2, title=gsub("."," ",stat,fixed=TRUE), xlab=defaults[stat,"xlab"], bins=bins, from=defaults[stat,"from"], to=defaults[stat,"to"])

  if (type == 'qq' | type == 'both')
    plotQQ(x=x, y=y, xlab=name.1, ylab=name.2, title=gsub("."," ",stat,fixed=TRUE))

}
