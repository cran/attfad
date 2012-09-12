clusterStatistics <-
function(expression.data, gene.cor, chip.cor) {
  cluster.results <- list()
  #gene.cluster <- hclust(as.dist(1-gene.cor), method='average')
  #chip.cluster <- hclust(as.dist(1-chip.cor), method='average')

  #cluster.results[['Gene cluster heights']] <- 1-gene.cluster$height
  #cluster.results[['Chip cluster heights']] <- 1-chip.cluster$height
  #cluster.results[['Gene - inner node cluster ranks']] <- seq(0,1,len=length(gene.cluster$height))[gene.cluster$merge[,1]*gene.cluster$merge[,2]<0]

  ### silhouettes
#   write.table(expression.data, file='/tmp/tmp32', sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
#   out <- system(paste('/home/m/maierr/workspace/expressionStats/cluster-1.35/src/cluster -m a -g 2 -f /tmp/tmp32'), intern=TRUE)
#   gene.silhouettes <- as.numeric(sapply(strsplit(out,'\t'),'[[',2))
  gene.silhouettes <- .Call("palcluster", 1-gene.cor)


#   write.table(t(expression.data), file='/tmp/tmp32', sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
#   out <- system(paste('/home/m/maierr/workspace/expressionStats/cluster-1.35/src/cluster -m a -g 2 -f /tmp/tmp32'), intern=TRUE)
#   chip.silhouettes <- as.numeric(sapply(strsplit(out,'\t'),'[[',2))
#   file.remove('/tmp/tmp32')
  chip.silhouettes <- .Call("palcluster", 1-chip.cor)

  cluster.results[['Silhouette.coefficient.genes']] <- gene.silhouettes
  cluster.results[['Silhouette.coefficient.microarrays']] <- chip.silhouettes
  cluster.results
}
