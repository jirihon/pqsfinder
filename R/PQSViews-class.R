###
## PQSViews class
##
## Author: Jiri Hon <jiri.hon@gmail.com>
## Date: 2016/01/17
## Package: pqsfinder
##


#' An S4 class to represent potential quadruplex forming sequences.
#'
#' Represents potential quadruplex forming sequences found by
#' \code{\link{pqsfinder}} function. This is a subclass of
#' \code{\link{XStringViews-class}} class and adds one more slot.
#'
#' @slot density Numbers of PQS (potential quadruplex forming sequences)
#'               overlapping at each position in input sequence.
#'
.PQSViews <- setClass(
  "PQSViews",
  contains = "XStringViews",
  slots = c(
    density = "numeric"
  ),
  validity = function(object) {
    if (length(object@subject) != length(object@density)) {
      return("Length of density vector is not equal to length of subject.")
    }
    return(TRUE)
  }
)


#' PQSViews class constructor.
#'
#' User friendly constructor for PQSViews class representing potential
#' quadruplex forming sequences (PQS). PQSViews is a subclass of
#' \code{\link{XStringViews}} class and adds one more slot to store
#' PQS density.
#'
#' @param subject DNAString object.
#' @param start Vector of PQS start positions.
#' @param width Vector of PQS lengths.
#' @param strand Vector of PQS strand specifications.
#' @param score Vector of PQS scores.
#' @param density Numbers of PQS overlapping at each position in \code{subject}.
#' @return PQSViews object
#'
#' @examples
#' pv <- PQSViews(DNAString("CGGGCGGGGC"), 1:2, 2:3, "+", 10:11, 1:10)
#' start(pv)
#' width(pv)
#' strand(pv)
#' score(pv)
#' density(pv)
#'
PQSViews <- function(
  subject, start, width, strand, score, density)
{
  ix <- order(start)

  .PQSViews(subject = subject, ranges = IRanges(start = start[ix], width = width[ix]),
            elementMetadata = DataFrame(strand = strand[ix], score = score[ix]),
            density = density)
}


#' Get PQS score vector
#'
#' @param x PQSViews object
#' @return Score vector
#' @examples
#' pqs <- pqsfinder(DNAString("CCCCCCGGGTGGGTGGGTGGGAAAA"))
#' score(pqs)
#'
setMethod("score", "PQSViews", function(x) mcols(x)$score)

#' Get PQS strand vector
#'
#' @param x PQSViews object
#' @return Strand vector
#' @examples
#' pqs <- pqsfinder(DNAString("CCCCCCGGGTGGGTGGGTGGGAAAA"))
#' strand(pqs)
#'
setMethod("strand", "PQSViews", function(x) mcols(x)$strand)

#' Get density vector
#'
#' Desity vector represents numbers of PQS (potential quadruplex forming
#' sequences) overlapping at each position in input sequence.
#'
#' @param x PQSViews object
#' @return Density vector
#' @examples
#' pqs <- pqsfinder(DNAString("CCCCCCGGGTGGGTGGGTGGGAAAA"))
#' density(pqs)
#'
setMethod("density", "PQSViews", function(x) x@density)


## The 2 helper functions below convert a given view on an XString object
## into a character-string.
##
## Both assume that 'start' <= 'end' (so they don't check it) and
## padd the result with spaces to produce the "margin effect"
## if 'start' or 'end' are out of limits.
##
## NOTE: Heavily based on Biostrings package, file XStringViews-class.R
##

## nchar(get_view(x, start, end)) is always end-start+1
##
.get_view <- function(x, start, end)
{
  as.character(subseq(x, start, end))
}

## nchar(get_snippet(x, start, end, snippetWidth)) is <= snippetWidth
##
.get_snippet <- function(x, start, end, snippetWidth, strand)
{
  if (snippetWidth < 7)
    snippetWidth <- 7
  width <- end - start + 1
  if (width <= snippetWidth) {
    .get_view(x, start, end)
  } else {
    w1 <- (snippetWidth - 2) %/% 2
    w2 <- (snippetWidth - 3) %/% 2
    paste(.get_view(x, start, start+w1-1),
          "...",
          .get_view(x, end-w2+1, end), sep="")
  }
}


## Show header of output table
##
.show_vframe_header <- function(iW, cols)
{
  cat(format("", width=iW)) # Print padding

  for (col in cols)
  {# Print column names
    cat(" ")
    cat(format(col$nm, width=col$width, justify="right"))
  }
  cat("\n")
}


## Show row of output table
##
.show_vframe_line <- function(x, i, iW, cols)
{
  # Print triplex index
  cat(format(paste("[", i,"]", sep=""), width=iW, justify="right"))

  colW <- 0 # Sum of all column width

  for (col in cols)
  {# Print column values
    cat(" ")
    value <- do.call(col$fn, list(x))[i]
    cat(do.call("format", c(value, col[3:length(col)], justify="right")))

    colW = colW + col$width
  }
  snippetW <- getOption("width") - iW - colW - length(cols) - 3
  cat(" [",
      .get_snippet(subject(x), start(x)[i], end(x)[i], snippetW),
      "]\n", sep="")
}


## Shot dots in place of hidden table rows
##
.show_vframe_line_dots <- function(x, iW, cols)
{
  cat(format("...", width=iW, justify="right"))

  for (col in cols) {
    cat(" ")
    cat(format("...", width=col$width, justify="right"))
  }
  cat("\n")
}


## Show all output table rows
## 'half_nrow' must be >= 1
##
.show_vframe <- function(x, half_nrow=9L)
{
  ## Column definitions
  ## nm = Column header
  ## fn = Function name to get row values
  ## ... other parameters passed to format function
  ##
  cols <- list(
    list(nm="start",  fn="start" ),
    list(nm="width",  fn="width" ),
    list(nm="score",  fn="score" ),
    list(nm="strand", fn="strand")
    #list(nm="pvalue", fn="pvalue", scientific=TRUE, digits=2),
  )

  i <- 1
  for (col in cols)
  {# Calculate column widths
    col_max <- max(do.call(col$fn, list(x)))
    col_maxstr <- do.call("format", c(col_max, col[3:length(col)]))
    cols[[i]] <- c(col, width=max(nchar(col_maxstr), nchar(col$nm)))
    i = i + 1
  }

  lx <- length(x)
  iW <- nchar(format(lx)) + 2 # Two extra for square brackets

  .show_vframe_header(iW, cols)

  if (lx <= 2*half_nrow + 1)
  {# Show all
    for (i in seq_len(lx))
      .show_vframe_line(x, i, iW, cols)
  }
  else
  {# Show first and last views
    for (i in 1:half_nrow)
      .show_vframe_line(x, i, iW, cols)

    .show_vframe_line_dots(x, iW, cols)

    for (i in (lx-half_nrow+1L):lx)
      .show_vframe_line(x, i, iW, cols)
  }
}


#' Show method
#'
#' @param object PQSViews object.
#' @return PQSViews object printed.
#'
setMethod("show", "PQSViews", function(object)
{
  subject <- subject(object)
  lsub <- length(subject)

  cat("  PQS views on a ", lsub, "-letter ", class(subject),
      " subject", sep="")
  cat("\nsubject:", .get_snippet(subject, 1, lsub, getOption("width")-9, "+"))
  cat("\nquadruplexes:")

  if (length(object) == 0) {
    cat(" NONE\n")
  }
  else {
    cat("\n")
    .show_vframe(object)
  }
})


###
## Coerce TriplexViews to DNAStringSet
##
setAs("PQSViews", "DNAStringSet", function(from)
{
  s <- DNAStringSet(subject(from), start(from), end(from))

  for (i in 1:length(s))
  {# Set proper names
    names(s)[i] <- paste(
      "start=", start(from)[i], ";",
      "end=", end(from)[i], ";",
      "strand=", strand(from)[i], ";",
      "score=", score(from)[i], ";",
      sep=""
    )
  }
  return(s)
})


###
## Coerce TriplexViews to GRanges
##
setAs("PQSViews", "GRanges", function(from)
{
  seqlen <- length(subject(from))
  names(seqlen) <- "chr1"

  GRanges(
    "chr1",
    IRanges(start(from), end(from)),
    strand(from),
    score = score(from),
    seqlengths = seqlen
  )
})


#' Coerce to character vector
#'
#' @param x PQSViews object.
#' @return Character vector representing PQS.
#'
setMethod("as.character", "PQSViews", function(x)
{
  s <- as(x, "DNAStringSet")
  as.character(s)
})


###
## Convert to printable string
##
setMethod("toString", "PQSViews", function(x, ...)
{
  toString(as.character(x), ...)
})
