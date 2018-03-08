library('TreeSearch'); library('ape'); library('phangorn')
library('Inapp') #  For parse.nexus
setwd("C:/Research/R/Hyoliths")
cd <- "../Hyoliths/"
files <- list.files(cd, pattern='mbank_.*\\.nex', full.names=TRUE)
filename <- files[which.max(file.mtime(files))]
cat("Reading data from", filename)
my_data <- Inapp::read.as.phydat(filename)
if (TRUE) {
#  cat("\n ! SUPRESSING TOMTELUVA") # !!!
#  my_data$Tomteluva_perturbata <- NULL # !!!
  cat("\n ! SUPRESSING NAMACALATHUS") # !!!
  my_data$Namacalathus <- NULL # !!!
}

iw_data <- PrepareDataIW(my_data)
nj.tree <- NJTree(my_data)
Fitch(nj.tree, my_data)
rooted.tree <- EnforceOutgroup(nj.tree, 'Dailyatia')
plot(rooted.tree, main="NJ tree")
par(mar=rep(0.25, 4), cex=0.75) # make plot easier to read

plot(better.tree <- TreeSearch(tree=rooted.tree, dataset=my_data, maxIter=3000,
                           EdgeSwapper=RootedNNISwap, verbosity=2))
text(0.5, 1.4, paste0("Score after first NNI swaps: ", Fitch(better.tree, my_data)), pos=4, cex=0.8)
plot(best.tree <- Ratchet(better.tree, my_data, verbosity=1,
                     swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap)))
text(0.5, 1.4, paste0("Score after first Ratchet: ", Fitch(best.tree, my_data)), pos=4, cex=0.8)
text(0.5, 0.5, Sys.time(), pos=4, cex=0.7)

i <- 0
bestScore <- 1e+9
if (FALSE) {
while (i < 22) {
  plot(best.tree <- Ratchet(best.tree, my_data, verbosity=1, searchHits=55,
                       swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap)))
  score <- Fitch(best.tree, my_data)
  i <- if (score + 1e-09 < bestScore) 0 else i + 1
  bestScore <- score
  text(0.5, 1.4, paste0("Score: ", score), pos=4)
  text(0.5, 0.5, Sys.time(), pos=4)
}
write.nexus(best.tree, file=paste0("hy_ew_", Fitch(best.tree, my_data), ".nex", collapse=''))

my.consensus <- RatchetConsensus(best.tree, my_data, nSearch=150, 
                                 swappers=list(RootedTBRSwap, RootedNNISwap))
par(mar=rep(0.25, 4), cex=0.75) # make plot easier to read
plot(ape::consensus(my.consensus))
no_acro <- lapply(my.consensus, drop.tip, tip='Clupeafumosus_socialis')
plot(ape::consensus(no_acro))
legend(0.5, 4, '? Clupeafumosus_socialis', lty=1, bty='n', text.font=3, cex=0.8)
text(0.5, 1.4, paste0("Score: ", Fitch(best.tree, my_data)), pos=4)                               
text(0.5, 0.5, Sys.time(), pos=4)

no_vagrants <- lapply(my.consensus, drop.tip, tip=
  c('Clupeafumosus_socialis', 'Yuganotheca_elegans', 'Longtancunella_chengjiangensis', 'Salanygolina', 'Tomteluva_perturbata')[
  c(3)])
plot(ape::consensus(no_vagrants))

# PDF
pdf(file=paste0("hy_ew_", Fitch(best.tree, my_data), ".pdf", collapse=''))
par(mar=rep(0.25, 4), cex=0.75)
plot(ape::consensus(no_vagrants))
legend(0.5, 4, c('? Longtancunella chengjiangensis'), lty=1, bty='n', text.font=3, cex=0.8)
text(0.5, 1.4, paste0("Score: ", Fitch(best.tree, my_data)), pos=4)                               
text(0.5, 0.5, Sys.time(), pos=4)

#plot(ape::consensus(no_acro))  # No ACROTRETIDS
#legend(0.5, 4, '? Clupeafumosus socialis', lty=1, bty='n', text.font=3, cex=0.8)
#text(0.5, 1.4, paste0("Score: ", Fitch(best.tree, my_data)), pos=4)                               
#text(0.5, 0.5, Sys.time(), pos=4)

plot(ape::consensus(my.consensus))
text(0.5, 1.4, paste0("Score: ", Fitch(best.tree, my_data)), pos=4)                               
text(0.5, 0.5, Sys.time(), pos=4)

dev.off()
}
iw.tree <- best.tree
kValues <- c(2, 3, 4.5, 7, 10.5, 16, 24)
latestFile <- vapply(kValues, function (k) {
  allFiles <- list.files(cd, pattern=paste0('hy_iw_k', k, '_.*\\.tre', sep=''), full.names=TRUE)
  allFiles[which.max(file.mtime(allFiles))]
  }, character(1))

for (k in kValues[order(file.mtime(latestFile))]) {

  cat("\n\n\n************************ Concavity constant k =", k, "**********************\n")
  attr(iw.tree, 'score') <- NULL
  iw.tree <- IWTreeSearch(tree=iw.tree, dataset=iw_data, concavity=k,
                             EdgeSwapper=RootedNNISwap, verbosity=3)
  i <- 0
  bestScore <- 1e+07
  while (i < 6) {
    cat("\n\n########################## I =", i, "#############")
    plot(iw.tree <- IWRatchet(iw.tree, iw_data, verbosity=2, concavity=k,
                         ratchHits = 4,
                         swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap)))
    score <- IWScore(iw.tree, iw_data, concavity=k)
    i <- if (score + 1e-09 < bestScore) 0 else i + 1
    text(0.5, 1.4, paste0("IW Score: ", score), pos=4)
    text(0.5, 0.5, Sys.time(), pos=4)
    bestScore <- score
  }
  score <- IWScore(iw.tree, iw_data, concavity=k)
  write.nexus(iw.tree, file=paste0("../Hyoliths/hy_iw_k", k, "_", round(score, 3), ".tre", collapse=''))

  suboptFraction = 0.02
  cat("\n\nEstimating consensus...")
  iw.consensus <- IWRatchetConsensus(iw.tree, iw_data, concavity=k,
                  swappers=list(RootedTBRSwap, RootedNNISwap),
                  searchHits=4,
                  suboptimal=score * suboptFraction,
                  nSearch=150, verbosity=1L) 
  plot(ape::consensus(iw.consensus))
  text(0.5, 1.4, paste0("k = ", k, "; IW Score: ", signif(score * (1 - suboptFraction), 5), "-", signif(score, 5)), pos=4)
  text(0.5, 0.5, Sys.time(), pos=4)  
  write.nexus(iw.consensus, file=paste0("../Hyoliths/hy_iw_k", k, "_", round(IWScore(iw.tree, iw_data, concavity=k), 3), ".con.tre", collapse=''))

  pdf(file=paste0("../Hyoliths/hy_iw_k", k, "_", round(IWScore(iw.tree, my_data, concavity=k), 3), ".pdf", collapse=''))
  par(mar=rep(0.25, 4), cex=0.75) # make plot easier to read
  plot(ape::consensus(iw.consensus))
  text(0.5, 1.4, paste0("k = ", k, "; IW Score: ", signif(score, 5), "-", signif(score * (1 + suboptFraction), 5)), pos=4)
  text(0.5, 0.5, Sys.time(), pos=4)  
  plot(iw.tree)
  text(0.5, 1.4, paste0("k = ", k, "; IW Score: ", signif(score, 5)), pos=4)
  text(0.5, 0.5, Sys.time(), pos=4)  
  dev.off()
}

