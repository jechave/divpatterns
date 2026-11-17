#' View a family in the browser
#'
#' Writes an interactive page for a fitted family. A matrix of the distance
#' from each homolog to the reference (homologs by alignment columns) is linked
#' to the 3D superposition and to the alignment of the chains shown: pointing
#' at a cell marks that residue in 3D, and clicking a row shows or hides that
#' homolog.
#'
#' The matrix is coloured by lrmsd, the log of the distance to the reference
#' (saturated at 0.25 and 4 Angstrom), or by nlrmsd, lrmsd minus the homolog's
#' mean lrmsd over the reference sites (saturated at a quarter and 4 times the
#' mean). Both use the viridis scale, from dark purple (low) to yellow (high).
#' Its columns are the whole alignment, or the reference sites: the columns
#' where the reference has a residue, numbered as the reference's residues. Rows
#' are ordered by average-linkage clustering on the root-mean-square difference
#' of their lrmsd over the reference sites (saturated as the colours are;
#' distances below 0.001 Angstrom count as 0.001). Cartoons need N, C and O
#' atoms; chains without them are drawn as Calpha traces.
#'
#' The page is one file with its data and the 3Dmol.js library inside, so it
#' works offline and can be sent as it is.
#'
#' @param fit A fit from [fit_to_reference()].
#' @param file Path of the HTML file to write. By default, a temporary file.
#' @param title Title of the page.
#' @param open Whether to open the page in the browser.
#'
#' @return The path of the file, invisibly.
#' @seealso [read_family()], [fit_to_reference()]
#' @export
#' @examples
#' \dontrun{
#' family <- read_family("family.fasta", files)
#' fit <- fit_to_reference(family, reference = "1abc_A")
#' view_family(fit, title = "My family")
#' }
view_family <- function(fit, file = NULL, title = "Family", open = interactive()) {
  if (!inherits(fit, "divpatterns_fit")) {
    stop("`fit` must come from `fit_to_reference()`.", call. = FALSE)
  }
  if (!is.character(title) || length(title) != 1) {
    stop("`title` must be one string.", call. = FALSE)
  }
  if (is.null(file)) file <- tempfile("family-", fileext = ".html")

  json <- jsonlite::toJSON(page_data(fit), na = "null", null = "null", digits = NA)
  # "</" in the data would end the script element that holds it; "<\/" is the same text in JSON
  json <- gsub("</", "<\\/", json, fixed = TRUE)

  viewer_dir <- system.file("viewer", package = "divpatterns", mustWork = TRUE)
  page <- readLines(file.path(viewer_dir, "template.html"), encoding = "UTF-8", warn = FALSE)
  library_code <- readLines(file.path(viewer_dir, "3Dmol-min.js"), encoding = "UTF-8", warn = FALSE)

  # the template has a line "__DATA__" for the data, a line "__3DMOL__" for the library, and "__TITLE__"
  if (sum(page == "__DATA__") != 1 || sum(page == "__3DMOL__") != 1) {
    stop("The page template is broken.", call. = FALSE)
  }
  if (any(grepl("</script", library_code, ignore.case = TRUE))) {
    stop("3Dmol-min.js cannot be placed inside the page.", call. = FALSE)
  }
  page[page == "__DATA__"] <- json
  page[page == "__3DMOL__"] <- paste(library_code, collapse = "\n")
  page <- gsub("__TITLE__", html_escape(title), page, fixed = TRUE)

  con <- file(file, open = "w", encoding = "UTF-8")
  writeLines(page, con, useBytes = TRUE)
  close(con)
  if (open) utils::browseURL(normalizePath(file))
  invisible(file)
}

# Everything the page shows, as a list to be written in JSON. Positions in the page's arrays start at 0.
page_data <- function(fit) {
  family <- fit$family
  ref <- fit$reference
  n_structures <- length(family$id)
  n_columns <- ncol(family$alignment)
  homologs <- setdiff(seq_len(n_structures), ref)
  site_columns <- which(!is.na(family$resno[ref, ]))
  cap <- log(4) # colours saturate at 4 times and a quarter of 1 Angstrom

  # lrmsd at the reference sites: one row per homolog. Distances below 0.001 Angstrom, the precision of
  # coordinates in PDB files, count as 0.001, so that the log is finite.
  lrmsd <- 0.5 * log(pmax(fit$d2[homologs, site_columns, drop = FALSE], 1e-6))
  rownames(lrmsd) <- family$id[homologs]
  mean_lrmsd <- rowMeans(lrmsd, na.rm = TRUE)
  order <- homologs[cluster_order(lrmsd, cap)]

  structures <- list()
  for (k in seq_len(n_structures)) {
    s <- list(
      id = jsonlite::unbox(family$id[k]),
      letter = family$alignment[k, ],
      resno = family$resno[k, ],
      insert = family$insert[k, ],
      resname = family$resname[k, ],
      ca = coordinates_by_column(family$ca[k, , ]),
      n = coordinates_by_column(family$n[k, , ]),
      c = coordinates_by_column(family$c[k, , ]),
      o = coordinates_by_column(family$o[k, , ]),
      distance = round(sqrt(fit$d2[k, ]), 3)
    )
    if (k != ref) {
      both <- family$alignment[ref, ] != "-" & family$alignment[k, ] != "-"
      s$identity <- jsonlite::unbox(round(mean(family$alignment[ref, both] == family$alignment[k, both]), 3))
      s$core <- jsonlite::unbox(length(fit$core[[k]]))
      s$shared <- jsonlite::unbox(sum(both))
      s$mean_lrmsd <- jsonlite::unbox(round(mean_lrmsd[[family$id[k]]], 4))
    }
    structures[[k]] <- s
  }

  list(
    cap = jsonlite::unbox(cap),
    colours = substr(viridisLite::viridis(256), 1, 7), # the colour scale from -cap to cap, as #rrggbb
    n_columns = jsonlite::unbox(n_columns),
    reference = jsonlite::unbox(ref - 1),
    order = order - 1,
    structures = structures
  )
}

# Coordinates (3 x columns) as a list with one entry per column: c(x, y, z), or NULL where there is no atom.
coordinates_by_column <- function(xyz) {
  out <- vector("list", ncol(xyz))
  for (a in seq_len(ncol(xyz))) {
    if (!is.na(xyz[1, a])) out[[a]] <- round(xyz[, a], 3)
  }
  out
}

# Order of the rows of `lrmsd` (homologs x reference sites): average-linkage clustering on the
# root-mean-square difference of two rows, over the sites both have, with values saturated at +-cap; within
# each merge, the branch with the lower mean lrmsd comes first.
cluster_order <- function(lrmsd, cap) {
  n <- nrow(lrmsd)
  if (n == 1) {
    return(1)
  }
  x <- pmax(pmin(lrmsd, cap), -cap)
  distance <- matrix(0, n, n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      both <- !is.na(x[i, ]) & !is.na(x[j, ])
      if (sum(both) < 10) {
        stop(sprintf("%s and %s share fewer than 10 reference sites.", rownames(lrmsd)[i], rownames(lrmsd)[j]), call. = FALSE)
      }
      distance[i, j] <- sqrt(mean((x[i, both] - x[j, both])^2))
      distance[j, i] <- distance[i, j]
    }
  }
  tree <- stats::as.dendrogram(stats::hclust(stats::as.dist(distance), method = "average"))
  tree <- stats::reorder(tree, wts = rowMeans(lrmsd, na.rm = TRUE), agglo.FUN = mean)
  stats::order.dendrogram(tree)
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}
