NexusTime <- function (filename, format='double') {
  FILE <- file(filename)
  open(FILE)
  comment <- readLines(FILE, n=3)[3]
  close(FILE)
  comment <- sub("\\-(\\d) ", "-0\\1 ", comment) # Morphobank do odd things with times!
  if (format == 'double') {
    as.double(sub('.*(\\d{4})\\-(\\d{2})\\-(\\d{2})\\s(\\d{2})\\.(\\d{2}\\.\\d{2}).*', "\\1\\2\\3\\4\\5", comment, perl=TRUE))
  } else {
    sub('.*(\\d{4}\\-\\d{2}\\-\\d{2}\\s\\d{2})\\.(\\d{2})\\.(\\d{2}).*', "\\1:\\2:\\3", comment, perl=TRUE)
  }
}

NewickTree <- function(tree) gsub('_', ' ', write.tree(tree), fixed=TRUE)

MatrixData <- function (states_matrix, state.labels) {
  tree <- states_matrix$tree
  regions <- states_matrix$regions
  changes <- states_matrix$changes
  n_tip <- states_matrix$n_tip
  plot.convert.state <- function(character, missing = FALSE) {
    plot.convert.inappli <- function(X) {
      return(ifelse(X == -1, "-", X))
    }
    plot.convert.missing <- function(X, all_states) {
      if (length(all_states) > 1 && length(X) == length(all_states) &&
          all(sort(X) == sort(all_states))) {
        return("?")
      }
      else {
        return(X)
      }
    }
    if (missing) {
      all_states <- unique(unlist(character))
      character <- lapply(character, plot.convert.missing,
                          all_states)
    }
    character <- lapply(character, plot.convert.inappli)
    return(unlist(lapply(character, function(X) paste(as.character(X),
                                                      collapse = ""))))
  }
  get.NA.edges <- function(states_matrix, tree, pass = 4) {
    check.applicable <- function(nodes, states_matrix, pass) {
      node1 <- states_matrix[[pass + 1]][nodes[1]][[1]]
      node2 <- states_matrix[[pass + 1]][nodes[2]][[1]]
      all_char <- sort(unique(unlist(states_matrix$Char)))
      options(warn = -1)
      node2 <- ifelse(all(node2 == all_char), node1, node2)
      options(warn = 0)
      return(ifelse(all(c(node1, node2) == -1), 0, 1))
    }
    return(apply(tree$edge, 1, check.applicable, states_matrix,
                 pass))
  }

  edge_col <- "black"
  tips_labels <- plot.convert.state(states_matrix[[1]][1:n_tip], missing = TRUE)

  tips_colours <- tips_labels
  tips_colours[nchar(tips_labels) > 1] <- "?"
  max_colour <- max(as.integer(tips_colours[tips_colours %in%
                                              0:9]))
  state_colours <- c(Inapp::brewer[[max_colour + 1]],
                     "grey")
  names(state_colours) <- c(0:max_colour, "?")
  if ("-" %in% tips_labels) state_colours <- c(state_colours, `-` = "lightgrey")
  edge_palette <- state_colours
  edge_palette["?"] <- "darkgrey"

  if (!is.null(unlist(states_matrix$Up2))) {
    na_edges <- get.NA.edges(states_matrix, tree, pass = 4) ==
      1
    edge_final <- ifelse(na_edges, "0", "-")
    edge_col <- ifelse(na_edges, "black", "grey")
  } else {
    edge_final = 0
  }
  #if (!is.null(unlist(states_matrix$Up1))) {
  if (!is.null(unlist(states_matrix$Up2))) {
    final_state <- states_matrix$Up2
  } else {
    final_state <- states_matrix$Up1
  }
  all_states <- -1:max_colour
  col_states <- c("-", 0:max_colour)
  colour.edge <- function(edge) {
    parent <- all_states %in% final_state[[edge[1]]]
    child <- all_states %in% final_state[[edge[2]]]
    common <- parent & child
    if (sum(common) == 1) {
      col_states[common]
    }
    else if (sum(child) == 1) {
      col_states[child]
    }
    else if (sum(parent) == 1 && !identical(parent,
                                            (col_states == "-"))) {
      col_states[parent]
    }
    else "?"
  }
  edge_final <- apply(tree$edge, 1, colour.edge)
  edge_col <- as.character(edge_palette[edge_final])
  #}

  if (length(state.labels) == length(edge_palette) - 2) {
    state.labels <- c(state.labels, "Ambiguous", "Inapplicable")
  } else if (length(state.labels) == length(edge_palette) - 1) {
    state.labels <- c(state.labels, "Ambiguous")
  } else if (length(state.labels) != length(edge_palette)) {
    warning("State labels do not seem to match states.  You need to label all states from 0 to the maximum observed.")
  }

  edge_col_array <- paste0('["', paste0(edge_col, collapse='", "'), '"]')

  state_labels <- paste(names(edge_palette), gsub("^['\"]|['\"]$",
                                                  "", state.labels), sep = ": ")
  observed <- names(edge_palette) %in% edge_final
  list (edge_col = edge_col,
        edge_col_array = edge_col_array,
        legend = state_labels[observed],
        legend_col = edge_palette[observed],
        tips_labels = tips_labels,
        tips_colours = as.character(state_colours[tips_colours]))
}


ReadNotes <- function (filepath) {
  lines <- enc2utf8(readLines(filepath))
  charNote.pattern <- "^\\s+TEXT\\s+CHARACTER=(\\d+)\\s+TEXT='(.*)';\\s*$"
  stateNote.pattern <- "^\\s+TEXT\\s+TAXON=(\\d+)\\s+CHARACTER=(\\d+)\\s+TEXT='(.*)';\\s*$"
  upperLines <- toupper(lines)
  notesStart <- which(upperLines == "BEGIN NOTES;")
  endBlocks <- which(upperLines == "ENDBLOCK;")
  if (length(notesStart) == 0) {
    return(list("NOTES block not found in Nexus file."))
  } else if (length(notesStart) > 1) {
    return(list("Multiple NOTES blocks found in Nexus file."))
  } else {
    notesEnd <- endBlocks[endBlocks > notesStart][1] - 1L
    notesLines <- lines[(notesStart + 1):notesEnd]
    charNote.matches <- grepl(charNote.pattern, notesLines, perl=TRUE)
    charNotes <- gsub(charNote.pattern, "\\2",
                      notesLines[charNote.matches], perl=TRUE)
    charNotes <- EndSentence(MorphoBankDecode(charNotes))
    charNumbers <- gsub(charNote.pattern, "\\1",
                        notesLines[charNote.matches], perl=TRUE)

    stateNote.matches <- grepl(stateNote.pattern, notesLines, perl=TRUE)
    stateNotes <- gsub(stateNote.pattern, "\\3",
                       notesLines[stateNote.matches], perl=TRUE)
    stateNotes <- EndSentence(MorphoBankDecode(stateNotes))
    stateTaxon <- gsub(stateNote.pattern, "\\1",
                       notesLines[stateNote.matches], perl=TRUE)
    stateChar  <- gsub(stateNote.pattern, "\\2",
                       notesLines[stateNote.matches], perl=TRUE)

    seqAlongNotes <- seq_len(max(as.integer(c(stateChar, charNumbers))))
    charNotes <- lapply(seqAlongNotes, function (i) {
      ret <- list(
        charNotes[charNumbers == i],
        stateNotes[stateChar == i])
      names(ret[[2]]) <- stateTaxon[stateChar == i]

      # Return:
      ret
    })
    names(charNotes) <- seqAlongNotes

    # Return:
    charNotes
  }
}

IsTransformational <- function (states) {
  gsub("^'(.*)'$", "\\1", states)[1] == ""
}

