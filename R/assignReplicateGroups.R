assignReplicateGroups <-
function(reps) {
  group <- 0
  last <- 2
  groups <- c()
  for (i in reps) {
    if(i <= last) {
      group <- group + 1
    }
    groups <- append(groups, group)
    last <- i
  }
  groups
}
