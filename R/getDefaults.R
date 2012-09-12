getDefaults <-
function() {
  defaults <- data.frame(bins=numeric(0), from=numeric(0), to=numeric(0), log=logical(0), enrichment=logical(0), stringsAsFactors=FALSE)
  defaults <- rbind(defaults, Replicate.noise.distribution       = data.frame(bins=30,  from=NA, to=NA, log=FALSE, enrichment=FALSE, tabID="1 ",  xlab="correlation"))
  defaults <- rbind(defaults, Absolute.intensity.distribution    = data.frame(bins=100, from=NA, to=NA, log=FALSE, enrichment=FALSE, tabID="2 ",  xlab="concentration"))
  defaults <- rbind(defaults, Range.of.gene.expression           = data.frame(bins=50,  from=NA, to=NA, log=FALSE, enrichment=FALSE, tabID="3 ",  xlab="90% range / median of all 90% ranges"))
  defaults <- rbind(defaults, Silhouette.coefficient.genes       = data.frame(bins=20,  from=0,  to=1,  log=FALSE, enrichment=FALSE, tabID="4 ",  xlab="silhouette scores" ))
  defaults <- rbind(defaults, Silhouette.coefficient.microarrays = data.frame(bins=30,  from=NA, to=NA, log=FALSE, enrichment=FALSE, tabID="4 ",  xlab="silhouette scores" ))
  defaults <- rbind(defaults, TF.TG.correlation		         = data.frame(bins=30,  from=-1, to=1,  log=FALSE, enrichment=FALSE, tabID="5a ", xlab="correlation"  ))
  defaults <- rbind(defaults, TF.TG.correlation.enrichment       = data.frame(bins=NA,  from=-1, to=1,  log=FALSE, enrichment=TRUE,  tabID="5a ", xlab="correlation"   ))
  defaults <- rbind(defaults, TG.TG.correlation                  = data.frame(bins=50,  from=-1, to=1,  log=FALSE, enrichment=FALSE, tabID="5b ", xlab="correlation"  ))
  defaults <- rbind(defaults, TG.TG.correlation.enrichment       = data.frame(bins=NA,  from=-1, to=1,  log=FALSE, enrichment=TRUE,  tabID="5b ", xlab="correlation"   ))
  defaults <- rbind(defaults, TF.activity.distribution           = data.frame(bins=20,  from=NA, to=NA, log=TRUE,  enrichment=FALSE, tabID="6 ",   xlab="log fraction of chips with target vs non-target at 0.05\nfor tfs with >= 5 targets"))
  defaults <- rbind(defaults, Gene.correlation                   = data.frame(bins=50,  from=-1, to=1,  log=FALSE, enrichment=FALSE, tabID="",    xlab="correlation"  ))
  defaults <- rbind(defaults, Chip.correlation                   = data.frame(bins=100, from=-1, to=1,  log=FALSE, enrichment=FALSE, tabID="",    xlab="correlation"  ))
  defaults
}
