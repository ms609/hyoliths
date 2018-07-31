# This file is copied from https://github.com/ms609/iotuba/blob/master/FormatRefs.R
# Please propogate any changes there.
refCon <- file('References.bib')
on.exit(close(refCon))
refLines <- readLines(refCon)
refLines <- gsub("^(\\w+\\s*=\\s*\\{)", "  \\1", refLines, perl=TRUE)
refLines <- gsub("\\{\\\\~\\{\\}\\}(\\w)", "\\\\textit{\\1", refLines, perl=TRUE)
refLines <- gsub("(\\w)\\{\\\\~\\{\\}\\}(\\W)", "\\1}\\2", refLines, perl=TRUE)
trimmedLines <- trimws(refLines)
remove <- substr(trimmedLines, 1, 9) %in% c('abstract ', 'keywords ') |
  substr(trimmedLines, 1, 7) %in% c('annote ') |
  substr(trimmedLines, 1, 6) %in% c('month ') |
  substr(trimmedLines, 1, 5) %in% c('file ', 'issn ', 'isbn ') |
  substr(trimmedLines, 1, 4) == 'url '
refLines <- refLines[!remove]

authLines <- substr(trimws(refLines), 1, 7) == 'author '
refLines[authLines] <- gsub(",(\\s+\\w)(?:[\\w']|\\{\\\\.\\{[\\w']\\}\\})*\\.?(?: ?(\\-?\\w)[\\w'\\.]*)?(?: ?(\\-?\\w)[\\w\\.]*)?(?: ?(\\-?\\w)[\\w'\\.]*)?( and|\\s*\\})", ",\\1. \\2. \\3. \\4. \\5", refLines[authLines], perl=TRUE)
refLines[authLines] <- gsub("\\.[\\s\\.]+", ". ", refLines[authLines], perl=TRUE)
refLines[authLines] <- gsub(". }", ".}", refLines[authLines], fixed=TRUE)
refLines[authLines] <- gsub(". -", ".-", refLines[authLines], fixed=TRUE)

firstChar <- substr(refLines, 1, 1)
keyLines <- firstChar == '@'
blankLines <- firstChar == ''
endLines <- firstChar == '}'

alphOrder <- order(gsub("^.+\\{(\\w+.*),", "\\1", refLines[keyLines], perl=TRUE))
whichKeys <- which(keyLines)
whichEnds <- which(endLines)
aBlank <- which(blankLines)[1]
sortedRefs <- refLines[unlist(lapply(alphOrder, function (i) c(whichKeys[i]:whichEnds[i], aBlank)))]

writeLines(sortedRefs, refCon)
close(refCon)
