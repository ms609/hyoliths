figHeight <- 525

if (!exists('allTrees')) source('loadTrees.R')
tipIndex <- sort(allTrees[[1]]$tip.label)
if (!exists('allSplits')) allSplits <- GetSplits(allTrees, tipIndex)
if (!exists('tntSplits')) {
  tnTrees <- lapply(tntFiles, ReadTntTree, relativePath='.', keepEnd=1)
  tnTrees <- lapply(tnTrees, as.multiPhylo)
  tnConsensi <- lapply(tnTrees, consensus)
  allTnt <- list()
  for (i in seq_along(tnTrees)) {
    allTnt <- c(allTnt, tnTrees[[i]])
  }
  tntSplits <- GetSplits(allTnt, tipIndex)
}

majCon <- RootTree(consensus(allTrees, p=0.5), rootingTips)
supporters <- SplitSupport(majCon, allSplits, tipIndex)
tnSupporters <- SplitSupport(majCon, tntSplits, tipIndex) / length(allTnt)
conLabel1 <- signif(100 * supporters / length(allTrees), 2)
conLabel2 <- ifelse(is.na(tnSupporters), '-', signif(100 * tnSupporters, 2))
conLabelColour1 <- NodeColour(round(conLabel1 / 100, 2), FALSE)
conLabelColour2 <- NodeColour(round(tnSupporters, 2), FALSE)

nConTip <- length(conTips)

conInternal <- nConTip + 1:Nnode(majCon)
conTerminal <- seq_len(nConTip)
conLabel <- character(nConTip)
conEdge <- dim(majCon$edge)[1]
eachConEdge <- seq_len(conEdge)
conParent <- majCon$edge[, 1]
conChild  <- majCon$edge[, 2]

conYStep <- round((figHeight - 7L) / nConTip, 1)

ancestors <- Ancestors(majCon, type='all')
outgroupNode <- vapply(ancestors, function (x) 56 %in% x, TRUE)
nOutgroupTips <- sum(outgroupNode) / 2L # nodes
lifters <- conTips %in% c("Amathia", 'Flustra')
nLifters <- sum(lifters)
nAncestors <- vapply(ancestors, length, 1)
conXStep <- round(seq(1L, to=svgWidth - 115L, len=max(nAncestors) + 1L), 1)
conX <- conXStep[nAncestors + 1L]
conY <- double(length(nAncestors))
lineHeight <- ((figHeight - 7L) - 10L) / (nConTip - 1 - nOutgroupTips - nLifters)
yLines <- seq(7 - (nOutgroupTips * lineHeight),
              figHeight - 7 + (nLifters * lineHeight),
              by=lineHeight)

conY[conChild[conChild <= nConTip]] <- yLines
conY[outgroupNode] <- conY[outgroupNode] + (25 * lineHeight)
conY[lifters] <- conY[lifters] - (14 * lineHeight)

for (nodeI in rev(conInternal)) {
  conY[nodeI] <- mean(conY[Children(majCon, nodeI)])
}
conY <- round(conY, 1)

edges <- paste0('<path d="', paste0('M', conX[conParent], ',', conY[conParent],
                                    'V', conY[conChild], 'H', conX[conChild],
                                    collapse=''), '" stroke="#888" fill="none"></path>')
tips <- paste0('<text x="', (conX[conTerminal] + 0L),
               '" y="', conY[conTerminal] + 4L,
               '" fill="', taxonColour[conTips[conTerminal]],
               '" class="taxonLabel', ifelse(conTips == 'Pedunculotheca_diania', ' bold', ''),
               '">',
               gsub('_', ' ', conTips[conTerminal], fixed=TRUE), '</text>',
               collapse='')
nodes <- paste0('<text x="', conX[conInternal] + 2L, '" y="', conY[conInternal] + 4L,
                '" fill="', conLabelColours ,'" class="nodeLabel">',
                conLabels, '</text>', collapse='')
tipLegend <- paste0('<g transform="translate(', svgWidth, ' ', figHeight, ')">',
                    '<text x="0" y="0" text-anchor="end" class="stepsLabel">',
                    paste0('<tspan x="0" dy="-1.2em" fill="', groupCol, '">',
                           names(groupCol), '</tspan>', collapse=''),
                    '<tspan x="0" dy="-1.2em" class="bold">Key to colours:</tspan>',
                    '</text></g>')
svgSource <- paste0('<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="',
                    svgWidth, '" height="', figHeight, '"><defs><style type="text/css">',
                    '<![CDATA[text.taxonLabel{font-style:italic}',
                    'text.nodeLabel{font-size:8pt}',
                    'svg {font-family: "Arial", sans-serif;font-size:9pt}',
                    '.bold{font-weight:bold}',
                    ']]></style></defs>',
                    tipLegend, tips, edges,  nodes, '</svg>')
write(svgSource, file="textFigure-raw.svg")
