ew.best <- list.files('TreeSearch', pattern='hy_ew_\\d*\\.nex', full.names=TRUE)
ew.trees <- read.nexus(file=ew.best[which.max(file.mtime(ew.best))])
ew.trees <- if (class(ew.trees) == 'multiPhylo') unique(ew.trees) else ew.trees
ew.trees <- as.multiPhylo(ew.trees)
tipIndex <- sort(ew.trees[[1]]$tip.label)
ew.trees <- lapply(ew.trees, RenumberTips, tipIndex)
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
    lapply(if (class(loadedTrees) == 'multiPhylo') unique(loadedTrees) else as.multiPhylo(loadedTrees),
           RenumberTips, tipIndex)

  }
})
iw.treesLoaded <- vapply(iw.trees, length, 0) > 0
iw.exist <- iw.trees[iw.treesLoaded]

allTrees <- ew.trees
for (i in seq_along(iw.exist)) {
  allTrees <- c(allTrees, iw.exist[[i]])
}
