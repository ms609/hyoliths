library('TreeSearch')
library('ape')
filename <- "X24932.nex"
# First replace all spaces in taxon names with underscores, remove ''s, and replace
# polymorphisms with "?".
raw_data <- ape::read.nexus.data(filename)

library(phangorn)
my_data <- phyDat(raw_data, type='USER', levels=c(0:9, '-'))

nj.tree <- NJTree(my_data)
Fitch(nj.tree, my_data)
rooted.tree <- EnforceOutgroup(nj.tree, 'Namacalathus')
excluded.taxa <- c('Orthis', 'Gasconsia', 'Glyptoria')
rooted.tree <- drop.tip(rooted.tree, excluded.taxa)
my_data[names(my_data) %in% excluded.taxa] <- NULL
plot(rooted.tree)
better.tree <- TreeSearch(tree=rooted.tree, dataset=my_data, 
                           EdgeSwapper=RootedNNISwap, verbosity=3)
best.tree <- Ratchet(better.tree, my_data, verbosity=0, k=5,
                     swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap))

for (i in 1:10) {
  plot(best.tree <- Ratchet(best.tree, my_data, verbosity=0,
                       swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap)))
  text(0.5, 0.5, paste0("Search iteration: ", i, "; Score: ", Fitch(best.tree, my_data)), pos=4)
}


my.consensus <- RatchetConsensus(best.tree, my_data, ratchHits = 100, nSearch=100, swappers=
                                 list(RootedTBRSwap, RootedNNISwap)) 
par(mar=rep(0.25, 4), cex=0.75) # make plot easier to read
plot(ape::consensus(my.consensus))                               
