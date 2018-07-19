# Included in TreeSearch v0.1.3.
SortTree <- function(tree) {
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
  tree$edge[, 1] <- parent
  tree$edge[, 2] <- child
  attr(tree, 'order') <- NULL
  Cladewise(Renumber(tree))
}

# This function should be replaced by the equivalent in TreeSearch 0.1.3.
SplitNumber <- function (tips, tr, tipIndex) {
  included <- tipIndex %in% tr$tip.label[tips]
  as.character(c(sum(powersOf2[included]), sum(powersOf2[!included])))
}

# This function should be replaced by the equivalent in TreeSearch 0.1.3.
SplitSupport <- function(tr, splitOccurrences, tipIndex) {
  nTip <- length(tipIndex)
  splits <- Descendants(tr, nTip + seq_len(Nnode(tr)), type='tips')
  splitNumbers <- vapply(splits, SplitNumber, character(2), tr, tipIndex)[1, ]
  splitOccurrences[splitNumbers]
}

# This function should be replaced by the equivalent in TreeSearch 0.1.3.
GetSplits <- function (trees, tipIndex) {
  nTip <- length(tipIndex)
  if (class(trees) == 'phylo') trees <- list(trees)
  table(vapply(trees, function (tr) {
    vapply(Descendants(tr, nTip + seq_len(nTip - 1L), type='tips'),
           SplitNumber, character(2), tr, tipIndex)
  }, character(2 * (nTip - 1L))))
}

as.multiPhylo <- phytools::as.multiPhylo

# This function will be included in TreeSearch 0.1.3.
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
