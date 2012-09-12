networkStatistics <-
function(expression.data, network, gene.cor, bins=50, alpha=0.05, min.tars=5) {

  network.results <- list()

  tfs <- list()
  for(i in unique(network[,2])) { tfs[[i]] <- sort(network[network[,2]==i,1]) }
  
  targets <- list()
  for(i in unique(network[,1])) { targets[[i]] <- sort(network[network[,1]==i,2]) }


  tg.tg <- matrix(NA, 0, 2)    
  for(i in names(which(sapply(tfs,length)!=0))) {
    for(j in names(which(sapply(tfs,length)!=0))) {
      if(j > i & length(tfs[[i]]) == length(tfs[[j]])) {
	if(all(tfs[[i]] == tfs[[j]]))
	  tg.tg <- rbind(tg.tg, c(i, j))
      }
    }
  }

  
  gene.cor.hist <- hist(gene.cor, seq(-1,1,len=bins), plot=F)
  tf.tg.cor.hist <- hist(c(gene.cor[network]), seq(-1,1,len=bins), plot=F)
  tg.tg.cor.hist <- hist(c(gene.cor[tg.tg]), seq(-1,1,len=bins), plot=F)

  network.results[['TF.TG.correlation']] <- c(gene.cor[network])
  network.results[['TG.TG.correlation']] <- c(gene.cor[tg.tg])

  network.results[['TF.TG.correlation.enrichment']] <- tf.tg.cor.hist$density-gene.cor.hist$density
  network.results[['TG.TG.correlation.enrichment']] <- tg.tg.cor.hist$density-gene.cor.hist$density

  #cat('begin tf.act\n')
  many.targets <- targets[sapply(targets,length) >= min.tars]
  tf.act <- sapply(many.targets,function(x) getTfActivity(t(scale(t(expression.data))),x,alpha))
  #cat('end tf.act\n')

  network.results[['TF.activity.distribution']] <- na.omit(tf.act)
  network.results
}
