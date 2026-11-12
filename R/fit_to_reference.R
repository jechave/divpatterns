#' Superpose each homolog on a reference
#'
#' Fits each structure of a family onto the reference on its own pair core,
#' and gives the distance from each homolog's Calpha to the reference's Calpha
#' in the same alignment column.
#'
#' The pair core of a homolog: fit it onto the reference on all the residues
#' both have, drop the residue furthest from its counterpart, and refit, until
#' every remaining residue lies within `cutoff` or `min_core` residues remain.
#' The homolog is then fitted on its core, and distances are given at every
#' residue both have. Fits are least-squares superpositions of Calpha atoms,
#' without weights.
#'
#' @param family A family from [read_family()].
#' @param reference The reference structure, by its position in the alignment
#'   or by its identifier.
#' @param cutoff Core cutoff, in Angstrom.
#' @param min_core Smallest core, in residues.
#'
#' @return A list of class `divpatterns_fit`. Entries with one value per
#'   structure are in the order of the alignment, and are `NA` (or `NULL`) for
#'   the reference.
#'   - `family`: the family, with the coordinates of each homolog moved by its
#'     fit. The reference does not move.
#'   - `reference`: the position of the reference in the alignment.
#'   - `d2`: squared distance from each structure's Calpha to the reference's,
#'     in square Angstrom: one row per structure, one column per alignment
#'     column, `NA` where either has a gap, and in the reference's row.
#'   - `core`: for each structure, the alignment columns of its core.
#'   - `reached`: for each structure, whether every core residue lies within
#'     `cutoff`.
#'   - `cutoff`, `min_core`: as given.
#'
#' @seealso [read_family()], [view_family()]
#' @export
#' @examples
#' \dontrun{
#' fit <- fit_to_reference(family, reference = "1abc_A")
#' }
fit_to_reference <- function(family, reference = 1, cutoff = 1.2, min_core = 15) {
  if (!inherits(family, "divpatterns_family")) {
    stop("`family` must come from `read_family()`.", call. = FALSE)
  }
  n_structures <- length(family$id)
  if (n_structures < 2) {
    stop("The family must hold at least two structures.", call. = FALSE)
  }
  if (is.character(reference)) reference <- match(reference, family$id)
  if (length(reference) != 1 || is.na(reference) || !(reference %in% seq_len(n_structures))) {
    stop("`reference` must be one structure of the family, by position or identifier.", call. = FALSE)
  }
  if (!is.numeric(cutoff) || length(cutoff) != 1 || !(cutoff > 0)) {
    stop("`cutoff` must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(min_core) || length(min_core) != 1 || !(min_core >= 3)) {
    stop("`min_core` must be at least 3.", call. = FALSE)
  }

  n_columns <- ncol(family$alignment)
  fixed <- family$ca[reference, , ] # 3 x columns
  d2 <- matrix(NA_real_, n_structures, n_columns)
  core <- vector("list", n_structures)
  reached <- rep(NA, n_structures)
  fitted <- family

  for (k in seq_len(n_structures)) {
    if (k == reference) next
    mobile <- family$ca[k, , ]
    shared <- which(!is.na(fixed[1, ]) & !is.na(mobile[1, ]))
    if (length(shared) <= min_core) {
      stop(sprintf("%s shares %d residues with the reference; more than %d are needed.", family$id[k], length(shared), min_core), call. = FALSE)
    }

    # the core, and the fit on it
    keep <- pair_core(fixed[, shared], mobile[, shared], cutoff, min_core)
    core[[k]] <- shared[keep]
    fit <- kabsch(fixed[, core[[k]]], mobile[, core[[k]]])

    # move every atom of the homolog by that fit
    for (atom in c("ca", "n", "c", "o")) {
      xyz <- family[[atom]][k, , ]
      present <- which(!is.na(xyz[1, ]))
      xyz[, present] <- apply_fit(fit, xyz[, present, drop = FALSE])
      fitted[[atom]][k, , ] <- xyz
    }

    moved <- fitted$ca[k, , ]
    d2[k, shared] <- colSums((fixed[, shared] - moved[, shared])^2)
    reached[k] <- max(sqrt(d2[k, core[[k]]])) <= cutoff
  }

  structure(
    list(family = fitted, reference = reference, d2 = d2, core = core, reached = reached, cutoff = cutoff, min_core = min_core),
    class = "divpatterns_fit"
  )
}

#' @export
print.divpatterns_fit <- function(x, ...) {
  homologs <- setdiff(seq_along(x$family$id), x$reference)
  cat(sprintf(
    "<divpatterns_fit> %d homologs fitted onto %s (pair cores within %g A, at least %d residues)\n",
    length(homologs), x$family$id[x$reference], x$cutoff, x$min_core
  ))
  cat(sprintf("  core sizes: %s\n", paste(range(lengths(x$core[homologs])), collapse = " to ")))
  missed <- homologs[!x$reached[homologs]]
  if (length(missed) > 0) {
    cat(sprintf("  cutoff not reached by: %s\n", paste(x$family$id[missed], collapse = ", ")))
  }
  invisible(x)
}

# Least-squares rotation of `mobile` onto `fixed` (3 x n each), with the centres to apply it.
kabsch <- function(fixed, mobile) {
  fixed_centre <- rowMeans(fixed)
  mobile_centre <- rowMeans(mobile)
  s <- svd((mobile - mobile_centre) %*% t(fixed - fixed_centre))
  d <- sign(det(s$v %*% t(s$u)))
  rotation <- s$v %*% diag(c(1, 1, d)) %*% t(s$u)
  list(rotation = rotation, fixed_centre = fixed_centre, mobile_centre = mobile_centre)
}

# Move coordinates `x` (3 x n) by a fit from kabsch().
apply_fit <- function(fit, x) {
  fit$rotation %*% (x - fit$mobile_centre) + fit$fixed_centre
}

# The pair core: the positions (columns of `fixed` and `mobile`) kept after dropping the furthest residue,
# one at a time, until all lie within `cutoff` or `min_core` remain.
pair_core <- function(fixed, mobile, cutoff, min_core) {
  keep <- seq_len(ncol(fixed))
  repeat {
    fit <- kabsch(fixed[, keep], mobile[, keep])
    distance <- sqrt(colSums((fixed[, keep] - apply_fit(fit, mobile[, keep]))^2))
    if (max(distance) <= cutoff || length(keep) <= min_core) {
      return(keep)
    }
    keep <- keep[-which.max(distance)]
  }
}
