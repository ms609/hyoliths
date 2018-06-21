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
               '"><tspan class="hidden">[', conTerminal, '] </tspan>',
               gsub('_', ' ', conTips[conTerminal], fixed=TRUE), '</text>',
               collapse='')
nodes <- paste0('<text x="', conX[conInternal][-1] + 2L, '" y="', conY[conInternal][-1] + 4L,
                '" class="nodeLabel"><tspan class="hidden">[', conInternal[-1],
                '] </tspan><tspan  fill="', conLabelColour1[-1] ,'">',
                conLabel1[-1], '</tspan>/<tspan  fill="', conLabelColour2[-1] ,'">',
                conLabel2[-1], '</tspan></text>', collapse='')
tipLegend <- paste0('<g transform="translate(', svgWidth, ' ', figHeight, ')">',
                    '<text x="0" y="0" text-anchor="end" class="stepsLabel">',
                    paste0('<tspan x="0" dy="-1.2em" fill="', groupCol, '">',
                           names(groupCol), '</tspan>', collapse=''),
                    '<tspan x="0" dy="-1.2em" class="bold">Key to colours:</tspan>',
                    '</text></g>')
halfEdge <- ((conXStep[2] - conXStep[1]) / 2)
onEdge <- conXStep + halfEdge
onNode <- conY + 4L
overEdge <-  round(conY - (lineHeight / 2), 2)
underEdge <- round(conY + (lineHeight / 2), 2)
onLine <- round(lineHeight * seq_len(nConTip) + 4L, 2)
legendKey <- conXStep[2] + halfEdge - 3
notes <- c(list(c(onEdge[6 ], overEdge[66], "Pe+")
                , c(onEdge[9 ], overEdge[98], "Pe-") # Hyolith loss
                , c(onEdge[12], overEdge[80], "Pe-") # Crani. loss
                , c(onEdge[9 ], overEdge[41], "Um") # Yugano.
                , c(onEdge[9 ], overEdge[10], "Um") # PEd
                , c(onEdge[12], overEdge[84], "Um") # Salany+
                , c(onEdge[17], overEdge[89], "Co") # Top 4 rhunchs
                , c(onEdge[15], overEdge[75], "D")  # Discinids & pals
                , c(onEdge[11], overEdge[83], "St") # Pater + Rhync
                , c(onEdge[10], overEdge[70], "*") # Crown
                , c(onEdge[9 ], overEdge[69], "L") # Main Brachs
                , c(onEdge[8 ], overEdge[68], "A") # Yugano + Brachs
                #, c(onEdge[3 ], overEdge[48], "c") # Namacal
                #, c(onEdge[11], underEdge[ 9], "c") # PAramicro
                #, c(onEdge[11], underEdge[ 5], "c") # Cupi
                #, c(onEdge[10], overEdge[95], "c,p") # Micrina gp
                #, c(onEdge[16], underEdge[79], "c") # Disc
                #, c(onEdge[18], overEdge[36], "c") # Disc
                #, c(onEdge[13], underEdge[40], "p") # Lingula
                #, c(onEdge[20], overEdge[20], "p") # Terebrata
                #, c(onEdge[15], overEdge[27], "p") # Novocrania
                #, c(onEdge[18], overEdge[78], "p") # Siph/Acanth
                #, c(conXStep[11] + 42L, lineHeight * 16.5,
                # '<tspan>Crown group</tspan><tspan dx="-6.2em" dy="1.2em">Brachiopoda</tspan>')
           )
#, lapply(2:21, function (x) c(conXStep[x], lineHeight * 1, x))
#, lapply(2:21, function (x) c(conXStep[x], lineHeight * 24, x))
#, lapply(2:21, function (x) c(conXStep[x], lineHeight * 47, x))
#, lapply(1:60, function (x) c(10, lineHeight * x, x))
#, lapply(1:60, function (x) c(svgWidth-20, lineHeight * x, x))
#, lapply(1:60, function (x) c(svgWidth/2, lineHeight * x, x))
)
key <- c(' ' = "Carbonate mineralogy"
         , Pe = "Pedicle"
         , Um = "Migration of pedicle to valve umbo"
         , Co = "Coelom lost in pedicle"
         , D = "Delthyrium surrounds pedicle"
         , St = "Strophic hinge; planar cardinal area"
         , A = "Anterior-coiling lophophore"
         , L = "Low ventral interarea; attached"
         , '  ' = 'setulose larvae'
         #, c = "Canaliculate structure"
         #, p = "Punctae"
         , '*' = "Brachiopod crown group"
         )

noteLegend <- paste0('<g id="legend" transform="translate(', onEdge[1], ' ',
                     round(lineHeight * 3, 2), ')"><text x="0" y="0">',
                     paste0(
                       sprintf('<tspan x="0" dy="%s" class="annot">%s</tspan>', round(lineHeight * 1.5, 2), names(key)),
                       sprintf('<tspan x="2em">%s</tspan>', key), collapse=''),
                     '</text></g>', collapse = '')

annotations <- paste0(vapply(notes, function (note)
  sprintf('\n\n<text x="%s" y="%s" class="annot">%s</text>\n\n', note[1], note[2], note[3]),
  character(1)), collapse='')

circles <- paste0(sprintf('<circle cx="%s" cy="%s"></circle>',
                          onEdge[c(1, 6, 8, 12, 13, 17)],
                          c(lineHeight * 3L - 4L, conY[c(43, 97, 80, 85, 31)])),
                  collapse='')

svgSource <- paste0('<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="',
                    svgWidth, '" height="', figHeight, '"><defs><style type="text/css">',
                    '<![CDATA[text.taxonLabel{font-style:italic}',
                    'text.nodeLabel{font-size:8pt}',
                    'svg {font-family: "Arial", sans-serif;font-size:9pt}',
                    '.bold,.annot{font-weight:bold}',
                    '.hidden{display:none}',
                    '.annot{text-anchor:middle;opacity:0.55}',
                    'circle{r:5px;fill:#f2e259;stroke:black;stroke-width:2px;}',
                    ']]></style></defs>',
                    tipLegend, noteLegend, tips, edges, circles, nodes, annotations, '</svg>')
write(svgSource, file="textFigure-raw.svg")
