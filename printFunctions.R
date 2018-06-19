Italicize <- function (string) {
  string <- gsub(paste0("\\b(",
                        paste0(gsub("_", "|", fixed = TRUE,
                                    c(taxa_names,
                                      'Bactrotheca', 'Bethia', 'Bowerbankia',
                                      'Calloria',
                                      'Disoma', 'Discinisca',
                                      'Eichwaldia', 'Electra', 'Enchytraeus',
                                      'Equisetum', 'Eoorthis', 'Fredericella',
                                      'Galeolaria', 'Glottidia', 'Gompholites',
                                      'Gryphus',
                                      'Hadrotreta', 'Harmothoe', 'Hydroides',
                                      'Hyolithes',
                                      'Ischnochiton', 'Kraussina', 'Lingulella',
                                      'Lochkothele', 'Loxosoma', 'Loxosomatoides',
                                      'Magelona', 'Membranipora', 'Mopalia',
                                      'Nereis', 'Neocrania', 'Notosaria',
                                      'Neoancistrocrania', 'Oikozetetes',
                                      'Paracraniops', 'Paramicrocornus',
                                      'Phascolion', 'Phoronopsis', 'Recilites',
                                      'Scolelepis',
                                      'Terebratalia', 'Themiste', 'Tubulipora',

                                      'atkinsae', 'comleyensis', 'cyrene', 'decaius',
                                      'deleta', 'minuta', 'murmanica',
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
        "  \n Unless you are very sure that you understand the consequences, ",
        "you should mark the character as Transformational by setting State 0 to",
        "`[Transformational character]`, or re-code: \n\n - ",
        paste(names(states[states == '-']), collapse="  \n - "))
}

GitLink <- function (gitSuffix, alt=NULL) {
  paste0(" [",
         if (!knitr::is_html_output() || is.null(alt)) {
           paste0(gsub("https://", "", rawGit, fixed=TRUE), gitSuffix)
         } else alt, "](",
         rawGit, gitSuffix, ")")
}

morphoBankRefereeAccess <- "<mark>[This dataset will be released on publication
of the paper. Referee access is availble by
[logging in to MorphoBank](https://morphobank.org/index.php/LoginReg/form)
using the e-mail address and password given in the manuscript.]</mark>"
