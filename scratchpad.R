```{R tnt-validate-results, echo=FALSE, eval=FALSE}
tntLog <- readLines('tnt.log')
scoreStarts <- which(grepl("Tree \\d+, total adjusted homoplasy", tntLog))
TntScores <- function (start) {
  scoreLines <- start + (3 * seq_len(ceiling(n_char / 5)))
  characterLines <- scoreLines - 1L
  firstBlank <- tntLog[characterLines] == " "
  if (any(firstBlank)) {
    beforeEnd <- seq_len(which(firstBlank) - 1L)
    scoreLines <- scoreLines[beforeEnd]
    characterLines <- characterLines[beforeEnd]
  }
  scores <- unlist(sapply(tntLog[scoreLines], function (line) {
    as.double(unlist(strsplit(line, "\\s+"))[-1])
  }))
  
  names(scores) <- unlist(sapply(tntLog[characterLines], function (line) {
    as.integer(unlist(strsplit(line, "\\s+"))[-1])
  }))
  # Return: 
  scores
}
  
tntScores <- lapply(scoreStarts, TntScores)
fullScores <- as.double(unlist(lapply(tntScores, sum)))
reportedScores <- as.double(gsub(".*?(\\d+(\\.\\d+)?).*", "\\1", tntLog[scoreStarts - 1L]))
k <- 2L
tntSteps <- round(k * tntScores[[1]] / (1 - tntScores[[1]]), 2)
if (all(tntSteps != as.integer(tntSteps))) warning("Step counting error?")
sum(tntSteps)

minSteps <- attr(fitch_data, 'min.steps')
weights <- attr(fitch_data, 'weight')
tr <- tnTrees[[1]][[1]]
treeSteps <- rep(FitchSteps(tr, fitch_data) - minSteps, weights)
treeSteps - tntSteps
sum(rep(FitchSteps(tr, fitch_data) - minSteps, weights))
rep(FitchSteps(tr, fitch_data) - minSteps, weights)
IWScore
lapply(seq_along(tnTrees), function (i) unique(unlist(lapply(tnTrees[[i]], IWScore, fitch_data, concavity=i))))
tntNames <- rep(kValues, lapply(tnTrees, length))
scores <- vapply(kValues, function (k) {
  ret <- vapply(unlist(tnTrees, recursive = FALSE), IWScore, double(1), fitch_data, concavity=k)
  names(ret) <- tntNames
  return(ret)
}, double(length(tntNames)))
colnames(scores) <- kValues
bestScore <- apply(scores, 2, min)
apply(scores, 2, function (col) rownames(scores)[col == min(col)])

```
