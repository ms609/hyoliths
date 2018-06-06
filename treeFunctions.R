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
  table(vapply(trees, function (tr) {
    vapply(Descendants(tr, nTip + seq_len(nTip - 1L), type='tips'),
           SplitNumber, character(2), tr, tipIndex)
  }, character(2 * (nTip - 1L))))
}
