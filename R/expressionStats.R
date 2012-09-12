expressionStats <-
function(expression.data, network, exp.group=NA, permute.genes=FALSE) {

  if(length(exp.group) <= 1)
    exp.group <- seq(nrow(expression.data))

  expression.data <- t(expression.data)
  if(permute.genes == TRUE) {
    n <- rownames(expression.data)
    expression.data <- apply(expression.data,2,sample)
    rownames(expression.data) <- n
  }
  network <- network[network[,1] %in% rownames(expression.data) & network[,2] %in% rownames(expression.data),1:2]
  rl <- list()
  
  ### pairwise correlations
  method <- ifelse(all(!is.na(expression.data)), 'na.or.complete', 'pairwise.complete.obs')

  gene.cor <- cor(t(expression.data), use=method)
  diag(gene.cor) <- NA
  rownames(gene.cor) <- rownames(expression.data)
  colnames(gene.cor) <- rownames(expression.data)

  chip.cor <- cor(expression.data, use=method)
  diag(chip.cor) <- NA

  rl[['Range.of.gene.expression']] <- levelsPerGene(expression.data)

  if(length(unique(exp.group)) != length(exp.group)) {
      ### replicates exist
      #rl[['Differential expression']] <- differentialExpression(expression.data, exp.group)
      rl[['Replicate.noise.distribution']] <- replicateGroupCorrelation(chip.cor, exp.group)
  }
  
  ### clustering
  #cat('clustering\n')
  rl <- c(rl, clusterStatistics(expression.data, gene.cor, chip.cor))
  
  ### tf-tg interactions
  #cat('tftg\n')
  rl <- c(rl, networkStatistics(expression.data, network, gene.cor))
  
  rl[['Absolute.intensity.distribution']] <- expression.data
  #rl[['Scaled absolute values']] <- scale(expression.data)
  rl[['Gene.correlation']] <- sample(gene.cor[upper.tri(gene.cor)],min(1e6,length(gene.cor[upper.tri(gene.cor)])))
  rl[['Chip.correlation']] <- chip.cor[upper.tri(chip.cor)]
  #rl[['Gene-sd/overall sd']] <- apply(expression.data,1,sd)/sd(c(expression.data))
  #cat('end es\n')
  rl

}
