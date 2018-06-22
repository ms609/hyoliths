FormatTree <- function(tree) {
  edge <- tree$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  tipLabels <- tree$tip.label
  tree.ntip <- length(tipLabels)
  descendants <- Descendants(tree)
  nDescendants <- vapply(descendants, length, integer(1))
  MinKid <- function (tips) min(tipLabels[tips])
  swaps <- vapply(tree.ntip + 1:Nnode(tree), function(node) {
    kids <- child[parent == node]
    descs <- nDescendants[kids]
    if (all(descs == 1L)) {
      order(tipLabels[kids])[1] == 1
    } else if (descs[1] == descs[2]) {
      order(vapply(descendants[kids], MinKid, character(1)))[1] == 1
    } else {
      descs[1] < descs[2]
    }
  }, logical(1))
  for (node in tree.ntip + rev(which(swaps))) {
    childEdges <- parent==node
    kids <- child[childEdges]
    child[childEdges][2:1] <- kids
  }
  #tree2 <- tree
  #tree2$edge <- cbind(parent, child)
  #dev.new(); par(mfrow=c(2, 1), cex=0.75, mar=rep(0.2, 4))
  #plot(tree); nodelabels(tree.ntip + rev(which(swaps)), node=tree.ntip + rev(which(swaps)), bg='green')
  #plot(tree2); nodelabels(tree.ntip + rev(which(swaps)), node=tree.ntip + rev(which(swaps)), bg='green')
  #tree2
  tree$edge[, 1] <- parent
  tree$edge[, 2] <- child
  tree
}

ew.best <- list.files('TreeSearch', pattern='hy_ew_\\d*\\.nex', full.names=TRUE)
ew.trees <- read.nexus(file=ew.best[which.max(file.mtime(ew.best))])
ew.trees <- if (class(ew.trees) == 'multiPhylo') unique(ew.trees) else ew.trees
tipIndex <- sort(ew.trees[[1]]$tip.label)
ew.trees <- as.multiPhylo(ew.trees)
ew.trees <- lapply(ew.trees, RenumberTips, tipIndex)
ew.trees <- lapply(ew.trees, FormatTree)
class(ew.trees) <- 'multiPhylo'
iw.trees <- lapply(kValues, function (k) {
  iw.best <- list.files('TreeSearch',
                        pattern=paste0('hy_iw_k',
                                       gsub('\\.', '\\\\.', k),
                                       '_\\d+\\.?\\d*\\.all\\.nex'),
                        full.names=TRUE)
  # Return:
  if (length(iw.best) == 0) {
    list()
  } else {
    loadedTrees <- read.nexus(iw.best[which.max(order(vapply(iw.best, ApeTime, character(1))))]) # TODO: SHould use dates from file write times
    ret <- lapply(if (class(loadedTrees) == 'multiPhylo') unique(loadedTrees) else as.multiPhylo(loadedTrees),
           RenumberTips, tipIndex)
    ret <- lapply(ret, FormatTree)
    class(ret) <- 'multiPhylo'
    ret
  }
})
iw.treesLoaded <- vapply(iw.trees, length, 0) > 0
iw.exist <- iw.trees[iw.treesLoaded]

allTrees <- ew.trees
for (i in seq_along(iw.exist)) {
  allTrees <- c(allTrees, iw.exist[[i]])
}
