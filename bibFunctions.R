BibEntry <- function (ref) {
  BibLine <- function (bibKey, key, text=ref[key]) {
    if (is.na(ref[key]) || ref[key] == "") {
      ""
    } else {
      ret <- paste0("  ", bibKey, " = {",
                    if (bibKey == 'journal') "{" else "",
                    gsub(" & ", " \\& ",
                         gsub("\u2013", "--",  # en dash
                              gsub("\u2014", "---",  # em dash
                                   gsub("ä", '{\\"a}',
                                        gsub("é", "{\\'e}",
                                             gsub("ö", '{\\"o}',
                                                  trimws(text),
                                                  fixed=TRUE),
                                             fixed=TRUE),
                                        fixed=TRUE),
                                   fixed=TRUE),
                              fixed=TRUE),
                         fixed=TRUE),
                    if (bibKey == 'journal') "}" else "",
                    "},\n")
      ret
    }
  }
  refType <- tolower(ref['Reference.Type'])
  refType <- if (refType %in% c('generic', 'journal', 'journal article')) 'article' else
    # if (refType == 'book section') 'inbook' else  # incollection would be nicer but throws errors?
    if (refType == 'book section') 'incollection' else
      if (refType == 'conference proceedings') 'incollection' else refType

  ref <- RefExpand(ref)
  paste0("@", refType, '{',
         RefKey(ref), ",\n",
         BibLine('author', 'Author', FormatAuthors(ref['Author'])),
         BibLine('doi', 'DOI'),
         BibLine('editor', 'Secondary.Author', FormatAuthors(ref['Secondary.Author'])),
         BibLine(if (refType == 'article') 'journal' else
           if (refType == 'inbook') 'title' else
             'booktitle', 'Secondary.Title'),
         #"  score = {", ref['score'], '; ', ref['relScore'], "},\n",
         BibLine(if (refType == 'inbook') 'chapter' else 'title', 'Title',
                 paste0('{', gsub("<i>", "\\\\textit{",
                                  gsub("</i>", "}", ref["Title"],
                                       perl=TRUE, ignore.case = TRUE),
                                  perl=TRUE, ignore.case=TRUE), '}')),
         BibLine('number', 'Number'),
         BibLine('pages', 'Pages', gsub("\\-+", "--", ref["Pages"], perl=TRUE)),
         BibLine('publisher', 'Publisher'),
         BibLine('volume', "Volume"),
         BibLine('year', 'Year'),
         "}\n"
  )
}

RefExpand <- function (ref) {

  # Manual overrides based on title cannot assume DOI present, so go first.
  title <- ref['Title']
  ref['Secondary.Title'] <-
    if (grepl("interpreting evolutionary relationships within early Rhynchonelliformea",
              title, fixed=TRUE) ||
        grepl("Dalarna, Sweden", title, fixed=TRUE)) {
      "Fossils and Strata"
    } else if (grepl("Lower Cambrian brachiopods from the rift valley",
                     title, ignore.case=TRUE) ||
               grepl("an early articulate brachiopod", title, fixed=TRUE)) {
      "Journal of Paleontology"
    } else if (grepl("Ultrastructure of spermatids and spermatozoa in",
                     title, fixed=TRUE)) {
      "Journal of Submicroscopic Cytology"
    } else if (grepl("palaeoecology of the trimerellid brachiopod", title, fixed=TRUE) ||
               grepl("shell structure in the classification of the orthotetidine",
                     title, fixed=TRUE) ||
               title == "The diversity and phylogeny of the paterinate brachiopods") {
      "Palaeontology"
    } else if (grepl("The early Cambrian tommotiid genus",
                     title, fixed=TRUE)) {
      "Memoirs of the Association of Australasian Palaeontologists"
    } else if (grepl("Early Cambrian (Botomian)", title, fixed=TRUE)) {
      "Acta Palaeontologica Polonica"
    } else if (grepl("in the Lower Ordovician sequence of Bohemia",
                     title, fixed=TRUE)) {
      "Sbornik geolgickych ved, Paleontologie"
    } else if (grepl("musculature and vascular systems of two species of Cambrian Paterinide",
                     title, fixed=TRUE)) {
      "Bureau of Mineral Resources Journal of Australian Geology and Geophysics"
    } else if (title == "The Cambrian radiation of brachiopods") {
      ref['Secondary.Author'] <- "Lipps, J. H.; Signor, Philip W."
      ref['Publisher'] <- "Pergamon"
      "Origin and Early Evolution of Metazoa"
    } else if (title == "Phoronida") {
      ref['Publisher'] <- "Wiley-Blackwell"
      ref['Secondary.Author'] <- "Harrison, F. W.; Woollacott, R. M."
      "Microscopic Anatomy of Invertebrates, 13: Lophophorates, Entoprocta, and Cycliophora"
    }  else ref['Secondary.Title']

  # Publisher:
  ref['Publisher'] <- if (grepl('Fish Evolution and Systematics', title, fixed=TRUE)) {
    "Cambridge University Press"
  } else ref['Publisher']


  CleanTitle <- function (title) {
    title <-gsub('</?i>|[.,;:\\-]', '', title, perl=TRUE)
    # Return:
    if (title=="") NULL else trimws(title)
  }
  RefAttribute <- function (name, pref='', suff='') {
    if (is.na(ref[name])) NULL else paste0(pref, ref[name], suff)
  }

  cr_details <- cached_dois(ref['DOI'])
  if (is.null(cr_details) || is.null(cr_details$data)) {
    flq <- c(query.bibliographic = paste(
      RefAttribute('Author'),
      RefAttribute('Year', '(', ')'),
      CleanTitle(ref['Title']),
      RefAttribute('Journal', '', ','),
      RefAttribute('Volume', ':'),
      RefAttribute('Pages'))
    )

    cr_details <-
      cached_works(query=CleanTitle(ref['Title']), flq=flq, limit=2,
                   select='container-title,DOI,title,issue,volume,score')
    if (cr_details$meta[, 'total_results'] > 0) {
      topResult <- cr_details$data[1, ]
      ref['score'] <- as.numeric(topResult$score)  ## TODO delete
      if (cr_details$meta[, 'total_results'] == 1 &&
          as.numeric(topResult$score) < 60) {
        return(ref)
      }
      relativeScore <- as.numeric(topResult$score) /
        as.numeric(cr_details$data[2, 'score'])
      ref['relScore'] <- round(relativeScore, 3)
      if (relativeScore < 1.3) {
        return (ref)
      }
    } else {
      # No results found.
      return(ref)
    }
  } else {
    topResult <- cr_details$data
  }

  AddIfBlank <- function (ref, cr_attr_name, attr) {
    if (
      (is.na(ref[attr]) || ref[attr] == '') &&
      cr_attr_name %in% names(topResult)
    ) {
      ref[attr] <- gsub("(\\w,)(\\w)", "\\1 \\2",
                        gsub("\\.$", "",
                             as.character(topResult[, cr_attr_name]), perl=TRUE),
                        perl=TRUE)
    }
    # Return:
    ref
  }
  ref <- AddIfBlank(ref, 'container.title', 'Secondary.Title')
  ref <- AddIfBlank(ref, 'doi', 'DOI')
  ref <- AddIfBlank(ref, 'issue', 'Issue')
  ref <- AddIfBlank(ref, 'title', 'Title')
  ref <- AddIfBlank(ref, 'volume', 'Volume')

  # Manual overrides
  if (ref['Title'] == "Invertebrate zoology: a functional evolutionary approach") {
    ref['Publisher'] <- "Thompson Learning"
    ref['DOI'] <- ref['volume'] <- ref['Secondary.Title'] <- '' # False positive
  }
  if (!is.null(ref['DOI']) && ref['DOI'] != "") {
    doi <- ref['DOI']
    ref['Secondary.Title'] <-
      if (grepl('10.1002/9781118896372', doi, fixed=TRUE)) {
        ref['Publisher'] <- "Blackwell"
        "The Cambrian Fossils of Chengjiang, China: The Flowering of Early Animal Life"
      } else if (grepl('10.17161/dt.v0i0', doi, fixed=TRUE)) {
        ref['Publisher'] <- "Geological Society of America & Paleontological Institute"
        "Treatise on Invertebrate Paleontology, Part H, Brachiopoda (Revised)"
      } else if (grepl('10.1201/9780203210437', doi, fixed=TRUE)) {
        ref['Secondary.Author'] <- "Brunton, H.; Cocks, R. M.; Long, S. L."
        ref['Publisher'] <- 'Taylor & Francis'
        "Brachiopods, Past and Present"
      } else if (grepl("setae and 3D MicroCT", title, ignore.case=TRUE)) {
        ref['DOI'] <- ''
        ref['Publisher'] <- "The Palaeontological Association"
        "Palaeontological Association Annual Meeting"
      } else ref['Secondary.Title']
  }
  # Return:
  ref
}

FirstWord <- function (string) {
  gsub("^(\\w+).*$", "\\1", trimws(string), perl=TRUE)
}

FirstWords <- function (string) {
  gsub("^\\W*(\\w+?)\\b\\W*(\\w+)?.*", "\\1\\2", trimws(string), perl=TRUE)
}

FormatAuthors <- function (string) {
  authors <- trimws(strsplit(string, ';')[[1]])
  names.pattern <- "(\\w+,)\\s(.*)"
  surnames <-  gsub(names.pattern, "\\1 ", authors)
  forenames <- gsub(names.pattern, "\\2", authors)
  forenames <- gsub("(\\w)\\w+\\.?", "\\1.", forenames)
  forenames <- gsub("\\.(\\w)", ". \\1", forenames)

  # Return:
  paste0(surnames, forenames, collapse=' and ')
}

RefKey <- function (ref) paste0(FirstWord(ref['Author']),
                                trimws(ref['Year']),
                                FirstWords(ref['Title']))

cached_works <- R.cache::addMemoization(rcrossref::cr_works)
cached_dois <- R.cache::addMemoization(function(doi)
  if (doi == '') NULL else rcrossref::cr_works(dois=doi))


StrReplace <- stringi::stri_replace_all_fixed

LookUpReference <- function (keys, original, prefix='', suffix='') {
  keyParts <- "^~~~(.+?)(\\b.*)~~~$"
  keySuffixes <- gsub(keyParts, "\\2", keys, perl=TRUE)
  keys <- gsub(keyParts, "\\1", keys, perl=TRUE)
  if (length(keys) == 0) {
    return (character(0))
  } else if (length(keys) == 1) {
    matches <- grepl(keys, refKeys, fixed=TRUE)
    if (sum(matches) == 1) {
      paste0(prefix, refKeys[matches], keySuffixes, suffix)
    } else {
      paste0('<mark>', original, '</mark>')
    }
  } else { # length(keys) > 1
    # initialize
    ret <- paste0('<mark>', original, '</mark>')
    if (length(suffix) == 1) suffix <- rep(suffix, length(keys))
    if (length(prefix) == 1) prefix <- rep(prefix, length(keys))

    matches <- vapply(keys, FUN=grepl, FUN.VALUE=logical(length(refKeys)), refKeys, fixed=TRUE)
    uniques <- colSums(matches) == 1
    ret[uniques] <- vapply(which(uniques), function (i)
      paste0(prefix[i], refKeys[matches[, i]], keySuffixes[i], suffix[i])
      , character(1))

    # Return:
    ret
  }
}

ConvertReferences <- function (string) {
  noParentheses <- "([A-Z][A-z]+)( et al\\.| ?(?:&|and) ?\\w+)?,? (\\d\\d\\d\\d\\w?)\\b"
  ReplNoPar <- function (keyedText, oldKey, original) {
    StrReplace(keyedText, oldKey,
               LookUpReference(oldKey, original, prefix='@'),
               vectorize_all=FALSE)
  }

  parentheses <- "([A-Z][A-z]+)( et al\\.| ?(?:&|and) ?\\w+)?,? \\((\\d\\d\\d\\d\\w?)\\b([^\\)]*)\\)"
  ReplPar <- function (keyedText, oldKey, original) {
    StrReplace(keyedText, oldKey,
               LookUpReference(oldKey, original,
                               prefix=paste0(gsub(parentheses, "\\1\\2", original, perl=TRUE), ' [-@'),
                               suffix=']'
               ), vectorize_all=FALSE)

  }

  ReplaceRefs <- function (pattern, strings, Replace) {
    execResults <- gregexpr(pattern, strings, perl=TRUE)
    matchText <- lapply(regmatches(strings, execResults), unique)
    keyed <- gsub(pattern, "~~~\\1\\3\\4~~~", strings, perl=TRUE)
    keys <- lapply(matchText, function (citation) {
      ret <- gsub(pattern, "~~~\\1\\3\\4~~~", citation, perl=TRUE)
      duplicates <- ret[duplicated(ret)]
      substr(ret[duplicated(ret)], nchar(duplicates)-6, nchar(duplicates)-3) <- '!x2!'
      ret
    })
    toEdit <- which(vapply(keys, length, 0) > 0)
    strings[toEdit] <- vapply(toEdit, function (i) {
      Replace(keyed[i], keys[[i]], matchText[[i]])
    }, character(1))
    # Return:
    strings
  }

  matches <- grepl(noParentheses, string)
  string[matches] <- ReplaceRefs(noParentheses, string[matches], ReplNoPar)

  matches <- grepl(parentheses, string)
  string[matches] <- ReplaceRefs(parentheses, string[matches], ReplPar)

  string <- gsub("\\(([^\\)]*?@\\S+?[^\\)]*)\\)", "[\\1]", string, perl=TRUE)

  # Return:
  string
}

MorphoBankDecode <- function (string) {
  string <- gsub("^n", "  \n", string, fixed=TRUE)
  string <- gsub("''", "'", string, fixed=TRUE)
  string <- gsub(" - ", " -- ", string, fixed=TRUE)
  string <- gsub("(\\d)\\-(\\d)", "\\1--\\2", string, perl=TRUE)
  string <- gsub("(\\d) ?um\\b", "\\1 µm", string, perl=TRUE)
  string <- gsub(" et al ", " et al. ", string, fixed=TRUE)
  string <- gsub(" et alia", " et al.", string, fixed=TRUE)

  # Return:
  ConvertReferences(string)
}
