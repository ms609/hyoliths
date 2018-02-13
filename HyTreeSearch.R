library('TreeSearch'); library('ape'); library('phangorn')
library('Inapp') #  For parse.nexus
filename <- "../Hyoliths/X24932.nex"
# First replace all spaces in taxon names with underscores, remove ''s, and replace
# polymorphisms with "?".
#raw_data <- ape::read.nexus.data(filename)
#my_data <- phyDat(raw_data, type='USER', levels=c(0:9, '-'))
my_data <- Inapp::read.as.phydat(filename)
#excluded.taxa <- c('Orthis', 'Gasconsia', 'Glyptoria')
nj.tree <- NJTree(my_data)
Fitch(nj.tree, my_data)
rooted.tree <- EnforceOutgroup(nj.tree, 'Namacalathus')
plot(rooted.tree)
#rooted.tree <- drop.tip(rooted.tree, excluded.taxa)
#my_data[names(my_data) %in% excluded.taxa] <- NULL
iw_data <- PrepareDataIW(my_data)


better.tree <- TreeSearch(tree=rooted.tree, dataset=my_data, 
                           EdgeSwapper=RootedNNISwap, verbosity=3)
plot(best.tree <- Ratchet(better.tree, my_data, verbosity=0, k=5,
                     swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap)))

for (i in 1:10) {
  plot(best.tree <- Ratchet(best.tree, my_data, verbosity=0,
                       swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap)))
  text(0.5, 0.5, paste0("Score: ", Fitch(best.tree, my_data)), pos=4)
}
write.nexus(best.tree, file=paste0("hy_ew_", Fitch(best.tree, my_data), ".nex", collapse=''))

my.consensus <- RatchetConsensus(best.tree, my_data, swappers=
                                 list(RootedTBRSwap, RootedNNISwap)) 
par(mar=rep(0.25, 4), cex=0.75) # make plot easier to read
plot(ape::consensus(my.consensus))                               
pdf(file=paste0("hy_ew_", Fitch(best.tree, my_data), ".pdf", collapse=''))
par(mar=rep(0.25, 4), cex=0.75) # make plot easier to read
plot(ape::consensus(my.consensus))
dev.off()

iw.tree <- best.tree
attr(iw.tree, 'score') <- NULL
k=4

iw.tree <- IWTreeSearch(tree=iw.tree, dataset=iw_data, concavity=k,
                           EdgeSwapper=RootedNNISwap, verbosity=3)
for (i in 1:100) {
  cat("\n\n############################# I =", i)
  plot(iw.tree <- IWRatchet(iw.tree, iw_data, verbosity=2, concavity=k,
                       swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap)))
  text(0.5, 0.5, paste0("IW Score: ", IWScore(iw.tree, iw_data, concavity=k)), pos=4)
}
write.nexus(iw.tree, file=paste0("../Hyoliths/hy_iw_k4_", round(IWScore(iw.tree, iw_data, concavity=k), 3), ".tre", collapse=''))

iw.consensus <- IWRatchetConsensus(iw.tree, iw_data, 
                swappers=list(RootedTBRSwap, RootedNNISwap),
                suboptimal=0.1,
                nSearch=20) 
par(mar=rep(0.25, 4), cex=0.75) # make plot easier to read
plot(ape::consensus(iw.consensus))                               
write.nexus(iw.consensus, file=paste0("../Hyoliths/hy_iw_k4_", round(IWScore(iw.tree, iw_data, concavity=k), 3), ".con.tre", collapse=''))

pdf(file=paste0("../Hyoliths/hy_iw_", round(IWScore(iw.tree, my_data, concavity=k), 3), ".pdf", collapse=''))
par(mar=rep(0.25, 4), cex=0.75) # make plot easier to read
plot(ape::consensus(iw.consensus))
dev.off()
