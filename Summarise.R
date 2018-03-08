library('TreeSearch'); library('ape'); library('phangorn')
cd <- "../Hyoliths/"
files <- list.files(cd, pattern='mbank_.*\\.nex', full.names=TRUE)
filename <- files[which.max(file.mtime(files))]
cat("Reading data from", filename)
my_data <- ReadAsPhyDat(filename)


tntFiles <- list.files(paste0(cd, 'TNT/'), pattern='xpiwe.*\\.tre', full=TRUE)
filename <- tntFiles[1]
tnTrees <- lapply(tntFiles, ReadTntTree)

plot(ConsensusWithout(tnTrees, c('Clupeafumosus_socialis', 'Micrina', 'Mickwitzia_muralensis')))

plot(ConsensusWithout(tnTrees, c('Clupeafumosus_socialis', 'Micrina', 'Mickwitzia_muralensis', 'Paterimitra', 'Heliomedusa_orienta', 'Tomteluva_perturbata', 'Yuganotheca_elegans', 'Salanygolina', 'Gasconsia', 'Mummpikia_nuda')))

plot(ConsensusWithout(tnTrees, c('Clupeafumosus_socialis', 'Micrina', 'Mickwitzia_muralensis', 'Paterimitra', 'Heliomedusa_orienta', 'Tomteluva_perturbata', 'Yuganotheca_elegans', 'Salanygolina', 'Gasconsia', 'Mummpikia_nuda', 'Kutorgina_chengjiangensis', 'Nisusia_sulcata')))
plot(ConsensusWithout(tnTrees, c('Clupeafumosus_socialis', 'Micrina', 'Mickwitzia_muralensis', 'Paterimitra', 'Heliomedusa_orienta', 'Tomteluva_perturbata', 'Yuganotheca_elegans', 'Salanygolina', 'Gasconsia', 'Mummpikia_nuda', 'Kutorgina_chengjiangensis', 'Nisusia_sulcata')))
