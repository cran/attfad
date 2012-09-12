replicateGroupCorrelation <-
function(chip.cor, exp.group) {
  grid <- expand.grid(exp.group, exp.group) 
  chip.cor[which(grid$Var1 != grid$Var2, arr.ind=T)] <- NA
  na.omit(chip.cor[upper.tri(chip.cor)])
}
