#' Add a LaTeX-typeset legend to a plotFun or plotPoints plot
#'
#' Adds a lattice-style legend/key to a trellis plot, such as the output of
#' `plotFun()` or `plotPoints()`. This function is intended to be called once
#' after the plot has been built. If the plot already has a legend, it is
#' removed and replaced by the new legend.
#'
#' Legend labels are processed with `latex2exp::TeX()`.
#'
#' Legend entries are passed through `...`. The name of each entry is the
#' legend label. Each entry may be either a color string or a named vector
#' with entries `col`, `lty`, `lwd`, and/or `pch`.
#'
#' @param ... Legend entries.
#' @param space Legend position. Common values are `"right"`, `"top"`,
#'   `"bottom"`, and `"left"`. Defaults to `"right"`.
#' @param columns Number of legend columns. If `NULL`, the default is `1`
#'   for legends on the left or right, and the number of legend entries for
#'   legends on the top or bottom.
#' @param background Background color for the legend. Defaults to `"white"`.
#' @param default_lty Default line type. Defaults to `1`.
#' @param default_lwd Default line width. Defaults to `1`.
#' @param default_pch Default plotting character. Defaults to `NA`, meaning no point.
#' @param plot A trellis plot object. Defaults to the most recent trellis plot.
#'
#' @return The updated trellis plot object.
#'
#' @examples
#' p <- plotFun(x^2 ~ x, xlim = c(-3, 3), col = "red", lty = 2, lwd = 3)
#' p <- plotFun(x^3 ~ x, add = TRUE, plot = p, col = "blue")
#'
#' add.legend(
#'   "$y=x^2$" = c(col = "red", lty = 2, lwd = 3),
#'   "$y=x^3$" = "blue",
#'   plot = p
#' )
#'
#' add.legend("$y$" = "maroon")
#'
#' @export
add.legend <- function(
    ...,
    space = "right",
    columns = NULL,
    background = "white",
    default_lty = 1,
    default_lwd = 1,
    default_pch = NA,
    plot = lattice::trellis.last.object()
) {
  entries <- unlist(list(...), use.names = TRUE)

  parsed <- parse.legend.entries(
    entries,
    default_lty = default_lty,
    default_lwd = default_lwd,
    default_pch = default_pch
  )

  if (is.null(columns)) {
    columns <- if (space %in% c("top", "bottom")) {
      length(parsed$label)
    } else {
      1
    }
  }

  key <- list(
    space = space,
    columns = columns,
    background = background,
    text = list(latex2exp::TeX(parsed$label))
  )

  if (any(!is.na(parsed$lty) & parsed$lty != 0)) {
    key$lines <- list(
      col = parsed$col,
      lty = parsed$lty,
      lwd = parsed$lwd
    )
  }

  if (any(!is.na(parsed$pch))) {
    key$points <- list(
      col = parsed$col,
      pch = parsed$pch
    )
  }

  # add.legend() is not incremental.
  # Any existing legend/key is removed before the new one is added.
  plot$legend <- NULL

  update(plot, key = key)
}


parse.legend.entries <- function(
    entries,
    default_lty = 1,
    default_lwd = 1,
    default_pch = NA
) {
  nms <- names(entries)

  attr_pattern <- "\\.(col|color|lty|lwd|pch)$"
  has_attr <- grepl(attr_pattern, nms)

  labels <- ifelse(
    has_attr,
    sub(attr_pattern, "", nms),
    nms
  )

  attrs <- ifelse(
    has_attr,
    sub("^.*\\.(col|color|lty|lwd|pch)$", "\\1", nms),
    "col"
  )

  attrs[attrs == "color"] <- "col"

  label_order <- unique(labels)

  col <- setNames(rep(NA_character_, length(label_order)), label_order)
  lty <- setNames(rep(default_lty, length(label_order)), label_order)
  lwd <- setNames(rep(default_lwd, length(label_order)), label_order)
  pch <- setNames(rep(default_pch, length(label_order)), label_order)

  for (i in seq_along(entries)) {
    label <- labels[i]
    attr <- attrs[i]
    value <- unname(entries[i])

    if (attr == "col") {
      col[label] <- value
    }

    if (attr == "lty") {
      lty[label] <- as.numeric(value)
    }

    if (attr == "lwd") {
      lwd[label] <- as.numeric(value)
    }

    if (attr == "pch") {
      pch[label] <- as.numeric(value)
    }
  }

  list(
    label = unname(label_order),
    col = unname(col),
    lty = unname(lty),
    lwd = unname(lwd),
    pch = unname(pch)
  )
}
