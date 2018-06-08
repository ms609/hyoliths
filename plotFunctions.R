#pal <- TreeSearch::brewer[[11]]
pal <- c("#F07894", "#E18758", "#C49800", "#99A700", "#55B23A", "#00B97B", "#00BAAE", "#00B3D7", "#52A2EF", "#AF8BF2", "#DE77DD") # <- colorspace::rainbow_hcl(11, 80, 65, 0, 308)

taxonColour <- c(
  Namacalathus                    = pal[1],
  Flustra                         = pal[1],
  Amathia                         = pal[1],
  Loxosomella                     = pal[1],
  Serpula                         = pal[11],
  Sipunculus                      = pal[11],
  Dentalium                       = pal[11],
  Tonicella                       = pal[11],
  Halkieria                       = 'grey',
  Halkieria_evangelista           = 'grey',
  Wiwaxia_corrugata               = pal[11],
  Wiwaxia                         = pal[11],
  Cotyledion_tylodes              = 'grey',
  Cotyledion                      = pal[1],
  Novocrania                      = pal[4],
  Craniops                        = pal[4],
  Ussunia                         = pal[4],
  Gasconsia                       = pal[4],
  Heliomedusa_orienta             = pal[3],
  Heliomedusa                     = pal[3],
  Micrina                         = pal[3],
  Mickwitzia_muralensis           = pal[3],
  Mickwitzia                      = pal[3],
  Micromitra                      = pal[8],
  Askepasma_toddense              = pal[8],
  Askepasma                       = pal[8],
  Pelagodiscus_atlanticus         = pal[6],
  Pelagodiscus                    = pal[6],
  Mummpikia_nuda                  = pal[5],
  Mummpikia                       = pal[5],
  Lingula                         = pal[5],
  Eoobolus                        = pal[5],
  Botsfordia                      = pal[5],
  Siphonobolus_priscus            = pal[5],
  Siphonobolus                    = pal[5],
  Clupeafumosus_socialis          = pal[5],
  Clupeafumosus                   = pal[5],
  Lingulellotreta_malongensis     = pal[5],
  Lingulellotreta                 = pal[5],
  Acanthotretella_spinosa         = pal[5],
  Acanthotretella                 = pal[5],
  Lingulosacculus                 = 'grey',
  Phoronis                        = pal[2],
  Dailyatia                       = pal[2],
  Eccentrotheca                   = pal[2],
  Yuganotheca_elegans             = 'grey',
  Yuganotheca                     = 'grey',
  Longtancunella_chengjiangensis  = 'grey',
  Longtancunella                  = 'grey',
  Paterimitra                     = pal[1],
  Pedunculotheca_diania           = 'black',
  Pedunculotheca                  = 'black',
  Haplophrentis_carinatus         = 'black',
  Haplophrentis                   = 'black',
  Tomteluva_perturbata            = 'grey',
  Tomteluva                       = 'grey',
  Salanygolina                    = pal[9],
  Coolinia_pecten                 = pal[9],
  Coolinia                        = pal[9],
  Antigonambonites_planus         = pal[9],
  Antigonambonites                = pal[9],
  Kutorgina_chengjiangensis       = pal[10],
  Kutorgina                       = pal[10],
  Nisusia_sulcata                 = pal[10],
  Nisusia                         = pal[10],
  Alisina                         = pal[10],
  Glyptoria                       = pal[10],
  Orthis                          = pal[10],
  Terebratulina                   = pal[10])


ColMissing <- function (omit) {
  MarkMissing(omit, text.font=3, cex=0.8, text.col=taxonColour[omit])
}

NodeColour <- function (support) {
  ifelse(is.na(support), 'red', divergingScale[(support * 100) + 1L])
}

# continuousScale <- rev(colorspace::heat_hcl(101, h=c(300, 75), c.=c(35, 95), l=c(15, 90), power=c(0.8, 1.2))) # Viridis prefered
divergingScale <- rev(colorspace::diverge_hcl(101, h=c(260, 0), c=100, l=c(50, 90), power=1.0))


WriteNumber <- function (n) {
  if (n == 0) {
    'zero'
  } else if (0 < n && n < 10) {
    c('one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine')[n]
  } else {
    n
  }
}

UnitEdges <- function (tree) {
  tree$edge.length <- rep(1, dim(tree$edge)[1])
  # Return:
  tree
}

SetPar <- function() par(mar=rep(0.2, 4), cex=0.8)

ColPlot <- function (tree) {
  SetPar()
  plot(tree, tip.color = taxonColour[tree$tip.label], lwd=2)
}

# Plot results for each value of k
PlotPanel <- function (treeList, i) {
  if (i > length(treeList) || length(treeList[[i]]) == 0) {
    cat("\n > Results not available for panel", i)
  } else {
    ColPlot(RootTree(consensus(treeList[[i]]), rootingTips))
    text(1, 1, paste0('k = ', kValues[i]), pos=4)
  }
}
