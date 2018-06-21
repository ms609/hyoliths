figHeight <- 525

rootingTips <- c('Tonicella', 'Serpula', 'Loxosomella')
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

if (!exists('bayesSplits')) {
  mbFiles <- list.files('MrBayes', pattern='^hyo\\..*\\.t$', full.names=TRUE)
  bayesTrees <- lapply(mbFiles, read.nexus)
  nBayesTrees <- length(bayesTrees[[1]])
  postBurnin <- lapply(bayesTrees, function (x) x[ceiling(nTrees * 0.1):nBayesTrees])
  bayesThinned <- c(postBurnin[[1]], postBurnin[[2]], postBurnin[[3]], postBurnin[[4]])
  class(bayesThinned) <- 'multiPhylo'
  sampleSize <- min(length(bayesThinned), 10000L) # Not all trees, but ought to be enough for our purposes
  bayesSplits <- GetSplits(sample(bayesThinned, sampleSize), tipIndex)
}

supporters <- SplitSupport(mbCon, mbSplits, tipIndex)
nodeSupport <- c('', signif(supporters[-1] / sampleSize * 100L, 3))

majCon <- RootTree(consensus(allTrees, p=0.5), rootingTips)
supporters <- SplitSupport(majCon, allSplits, tipIndex)
tnSupporters <- SplitSupport(majCon, tntSplits, tipIndex) / length(allTnt)
bayesSupporters <- SplitSupport(majCon, bayesSplits, tipIndex) / sampleSize
conLabel1 <- signif(100 * supporters / length(allTrees), 2)
#conLabel2 <- ifelse(is.na(tnSupporters), '0', signif(100 * tnSupporters, 2))
conLabel2 <- ifelse(is.na(bayesSupporters), '0', round(100 * bayesSupporters))
conLabelColour1 <- NodeColour(round(conLabel1 / 100, 2), TRUE)
conLabelColour2 <- NodeColour(round(bayesSupporters, 2), TRUE)
conLabel1[conLabel1 == 100] <- '.'
conLabel2[conLabel2 == 100] <- '.'

conTips <- majCon$tip.label
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
tips <- paste0('<text x="',  conX[conTerminal],
               '" y="', 4L + conY[conTerminal],
               '" fill="', taxonColour[conTips[conTerminal]],
               '" class="taxonLabel', ifelse(conTips == 'Pedunculotheca_diania', ' bold', ''),
               '">',
               gsub('_', ' ', conTips[conTerminal], fixed=TRUE), '</text>',
               collapse='')
nodes <- paste0('<text x="', conX[conInternal][-1] + 2L, '" y="', conY[conInternal][-1] + 4L,
                '" class="nodeLabel"><tspan  fill="', conLabelColour1[-1] ,'">',
                conLabel1[-1], '</tspan>/<tspan  fill="', conLabelColour2[-1] ,'">',
                conLabel2[-1], '</tspan></text>', collapse='')
tipLegend <- paste0('<g transform="translate(', svgWidth, ' ', figHeight, ')">',
                    '<text x="0" y="0" text-anchor="end" class="stepsLabel">',
                    paste0('<tspan x="0" dy="-1.2em" fill="', groupCol, '">',
                           names(groupCol), '</tspan>', collapse=''),
                    '<tspan x="0" dy="-1.2em" class="bold">Key to colours:</tspan>',
                    '</text></g>')
notes <- c(list(c(conXStep[11] + 5L, lineHeight * 16.5, '<tspan>Crown group</tspan><tspan dx="-6.2em" dy="1.2em">Brachiopoda</tspan>'),
              c(0,0,0)),
lapply(2:21, function (x) c(conXStep[x], lineHeight * 1, x)),
lapply(2:21, function (x) c(conXStep[x], lineHeight * 24, x)),
lapply(2:21, function (x) c(conXStep[x], lineHeight * 47, x)),
lapply(1:60, function (x) c(10, lineHeight * x, x)),
lapply(1:60, function (x) c(svgWidth-20, lineHeight * x, x)),
lapply(1:60, function (x) c(svgWidth/2, lineHeight * x, x))
)
annotations <- paste0(vapply(notes, function (note)
  sprintf('\n\n<text x="%s" y="%s" class="annotation">%s</text>\n\n', note[1], note[2], note[3]),
  character(1)), collapse='')

CircleX <- function (taxa) mean(conX[which(conTips %in% taxa)]) - 4L
CircleY <- function (taxa) mean(conY[which(conTips %in% taxa)]) + 4L

circles <- paste0(sprintf('<circle cx="%s" cy="%s"></circle>',
                          conXStep[c(6)] + ((conXStep[2] - conXStep[1]) / 2),
                          vapply(list('Halkieria_evangelista'), CircleY, double(1)) - 4L),
                  collapse='')

svgSource <- paste0('<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="',
                    svgWidth, '" height="', figHeight, '"><defs><style type="text/css">',
                    '<![CDATA[text.taxonLabel{font-style:italic}',
                    'text.nodeLabel{font-size:8pt}',
                    'svg {font-family: "Arial", sans-serif;font-size:9pt}',
                    '.bold,.annotation{font-weight:bold}',
                    '.annotation{fill:#444}',
                    'circle{r:5px;fill:#f2e259;stroke:black;stroke-width:2px;}',
                    ']]></style></defs>',
                    tipLegend, tips, edges, nodes, annotations, circles, '</svg>')
write(svgSource, file="textFigure-raw.svg")
