# Sorts each node into a consistent order, so similar trees look similar.
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
    # if (all(kids <= tree.ntip)) child[childEdges][2:1] <- kids
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
  attr(tree, 'order') <- NULL
  Cladewise(Renumber(tree))
}

SplitNumber <- function (tips, tr, tipIndex) {
  included <- tipIndex %in% tr$tip.label[tips]
  as.character(c(sum(powersOf2[included]), sum(powersOf2[!included])))
}

SplitSupport <- function(tr, splitOccurrences, tipIndex) {
  nTip <- length(tipIndex)
  splits <- Descendants(tr, nTip + seq_len(Nnode(tr)), type='tips')
  splitNumbers <- vapply(splits, SplitNumber, character(2), tr, tipIndex)[1, ]
  splitOccurrences[splitNumbers]
}

GetSplits <- function (trees, tipIndex) {
  nTip <- length(tipIndex)
  if (class(trees) == 'phylo') trees <- list(trees)
  table(vapply(trees, function (tr) {
    vapply(Descendants(tr, nTip + seq_len(nTip - 1L), type='tips'),
           SplitNumber, character(2), tr, tipIndex)
  }, character(2 * (nTip - 1L))))
}

as.multiPhylo <- phytools::as.multiPhylo

# This function will be included in a future release of TreeSearch.
TNTText2Tree <- function (treeText) {
  treeText <- gsub("(\\d+)", "\\1,", treeText, perl = TRUE)
  treeText <- gsub(")(", "),(", treeText, fixed = TRUE)
  treeText <- gsub("*", ";", treeText, fixed = TRUE)
  # Return:
  read.tree(text = gsub(", )", ")", treeText, fixed = TRUE))
}

GetJacks <- function (jackFile) {
  jackLines <- readLines(jackFile)
  jackTree <- TNTText2Tree(jackLines[3])
  jackTipOrder <- order(as.integer(jackTree$tip.label) + 1L)
  jackNodeOrder <- unique(unlist(Ancestors(jackTree, jackTipOrder)))[-1]
  nTntNode <- jackTree$Nnode

  treeFile <- gsub("\\.sym$", ".tre", jackFile)
  tntTrees <- ReadTntTree(treeFile, relativePath='.')
  tipLabel <- if (class(tntTrees) == 'phylo') tntTrees$tip.label else tntTrees[[1]]$tip.label
  jackTree$tip.label <- tipLabel

  jackScores <- trimws(gsub("ttag \\+\\d+ (.*); *", "\\1",
                            jackLines[length(jackLines) - (nTntNode - 2L):0]))[order(jackNodeOrder)]
  return (list(freq=gsub("^(\\d+)/.*", "\\1", jackScores),
              gc=gsub("^\\d+/(\\[?\\d+\\]?)$", "\\1", jackScores),
              tree=jackTree))
}
