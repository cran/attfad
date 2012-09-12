levelsPerGene <-
function(expression.data) {
  gene.ranges <- apply(expression.data,1, function(x) quantile(x,0.95)-quantile(x,0.05))
  gene.ranges/median(gene.ranges)
}
