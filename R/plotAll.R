plotAll <-
function(results.1, results.2, name.1, name.2, type='both', outdir='./', format='') {
  for(i in names(results.1)) {
    underscore.name <- gsub('[.]','_',i)
    if(i %in% names(results.2) & i != 'expression.data') {
      if(format == '')
	dev.new()
      if(format == 'png') {
	dir.create(outdir, showWarnings=FALSE)
	png(paste(outdir,'/', underscore.name,'.png', sep=''), ifelse(type=='both',1200,700),700)
      }
      if(format == 'pdf') {
	dir.create(outdir, showWarnings=FALSE)
	pdf(file=paste(outdir,'/', underscore.name,'.pdf', sep=''), width=ifelse(type=='both',12,7), height=7)
      }
      plotComparison(results.1, results.2, name.1, name.2, i, type=type)
      if(format == 'png' | format == 'pdf')
	dev.off()
    }
  }
}
