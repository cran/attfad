differentialExpression <-
function(expression.data, exp.group, cutoff=2.5) {
  gene.sds <- apply(expression.data,1,function(x) getGeneSd(x,exp.group))
  rep.z.score <- matrix(NA, nrow=nrow(expression.data), ncol=length(unique(exp.group)))
  for(i in 1:length(unique(exp.group))) {
    #cat('group ', i, '\r')
    rep.z.score[,i] <- apply(expression.data,1, function(x) sqrt(sum(exp.group==unique(exp.group)[i]))*mean(x[exp.group==unique(exp.group)[i]]-mean(x)))/gene.sds 

    ### multiple testing correction
    rep.p.val <- pnorm(abs(rep.z.score), lower.tail=FALSE)*2
    rep.p.val <- apply(rep.p.val, 2, function(x) p.adjust(x, method='BH'))
    rep.z.score <- qnorm(1-rep.p.val)

    apply(rep.z.score,1, function(x) mean(abs(x)>cutoff))
  }
}
