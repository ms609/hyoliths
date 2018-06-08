Italicize <- function (string) {
  string <- gsub(paste0("\\b(",
                        paste0(gsub("_", "|", fixed = TRUE,
                                    c(taxa_names,
                                      'Bethia', 'Bowerbankia', 'Calloria',
                                      'Disoma', 'Discinisca',
                                      'Eichwaldia', 'Electra', 'Enchytraeus',
                                      'Equisetum', 'Eoorthis', 'Fredericella',
                                      'Galeolaria', 'Glottidia', 'Gryphus',
                                      'Hadrotreta', 'Harmothoe', 'Hydroides',
                                      'Ischnochiton', 'Kraussina',
                                      'Lochkothele', 'Loxosoma', 'Loxosomatoides',
                                      'Magelona', 'Membranipora', 'Mopalia',
                                      'Nereis', 'Neocrania', 'Notosaria',
                                      'Oikozetetes',
                                      'Paracraniops', 'Paramicrocornus',
                                      'Phascolion', 'Phoronopsis', 'Scolelepis',
                                      'Terebratalia', 'Themiste', 'Tubulipora',

                                      'atkinsae', 'comleyensis', 'decaius',
                                      'minuta', 'murmanica',
                                      'schucherti', 'tenuis', 'worsleyi',

                                      'vascula', 'genitalia','lateralia',
                                      'media', 'myaria', 'terminalia',
                                      'levator ani',
                                      'et al\\.', 'et al', 'sensu')),
                               collapse='|'),
                        ")\\b"), "_\\1_", string, perl=TRUE, ignore.case=FALSE)
  gsub("\\b([A-Z]\\.) _([a-z])", "_\\1 \\2", string, perl=TRUE)
}

PrintStates <- function (states) {
  states <- gsub("^'(.*)'$", "\\1", states)
  tokens <- seq_along(states) - 1L
  if (states[1] == "") {
    states <- states[-1]
    tokens <- tokens[-1]
    transformational <- TRUE
  } else {
    transformational <- FALSE
  }
  cat(paste0(" > ", tokens, ": ", states, "  \n"))
  cat("> ", if (transformational) "Transformational" else "Neomorphic", "character.  \n>\n")
}

PrintNaughtyInapplicables <- function (states) {
  if (any(states == '-'))
    cat("  \n Oh dear! <mark>**You included the inapplicable token in a neomorphic character!**</mark>",
        "  \n That's really very naughty, as @Brazeau2018 will tell you.",
        "  \n Unless you are very sure that you understand the consequences, you should re-code  \n\n - ", paste(names(states[states == '-']), collapse="  \n - "))
}
