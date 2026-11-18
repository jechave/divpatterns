#' Read a family of homologous structures
#'
#' Reads a sequence alignment and one structure file per aligned sequence, and
#' places every residue of every structure in its alignment column. A table
#' says which file holds the structure of each sequence.
#'
#' Each sequence of the alignment, without its gaps, must be the sequence of
#' the residues in its file, in the order of the file. The first residue of
#' the file goes to the column of the first letter of its sequence, the second
#' residue to the column of the second letter, and so on. A residue is an
#' amino acid with a Calpha atom; other records, such as water or ligands, are
#' ignored. Where a file gives alternate locations for an atom, the first is
#' used.
#'
#' @param alignment Path to a FASTA alignment.
#' @param files A data frame with columns `id` and `file`: one row per
#'   sequence of the alignment, with its identifier, as in the alignment, and
#'   the path to its structure file, in PDB format. Each file must hold one
#'   chain. The rows can be in any order; other columns are ignored.
#' @param properties Optional. Values of residue properties, such as pLDDT or
#'   B-factors, shown in [view_family()] (one named `plddt` can mask residues
#'   of low confidence): a data frame with columns
#'   `id`, `resno` and `insert`, which identify a residue by the id of its
#'   sequence and its residue number and insertion code in its file (`NA` or
#'   `""` for none), and then one numeric column per property, named by the
#'   property. Each row must be a residue of the family; residues without a
#'   row, and `NA` values, have no value. The rows can be in any order.
#'
#' @return A list of class `divpatterns_family`. Its matrices have one row
#'   per structure and one column per alignment column:
#'   - `id`: the sequence identifiers of the alignment.
#'   - `file`: the structure files, in the order of the alignment.
#'   - `alignment`: the aligned sequences, one letter per residue, `"-"` at
#'     gaps.
#'   - `resno`, `insert`, `resname`: residue number, insertion code and
#'     residue name, as in the files; `NA` at gaps, and `insert` is `NA` where
#'     a residue has no insertion code.
#'   - `ca`, `n`, `c`, `o`: coordinates of the Calpha, N, C and O atoms,
#'     arrays indexed by structure, coordinate (x, y, z) and column; `NA` at
#'     gaps and where a file lacks the atom.
#'   - `properties`: one matrix per property, named as its column in
#'     `properties`; `NA` at gaps and where a residue has no value. An empty
#'     list without `properties`.
#'
#' @seealso [fit_to_reference()], [view_family()]
#' @export
#' @examples
#' \dontrun{
#' files <- data.frame(
#'   id = c("1abc_A", "2def_B", "3ghi_A"),
#'   file = c("1abc_A.pdb", "2def_B.pdb", "3ghi_A.pdb")
#' )
#' family <- read_family("family.fasta", files)
#' }
read_family <- function(alignment, files, properties = NULL) {
  aln <- bio3d::read.fasta(alignment, rm.dup = FALSE) # bio3d would drop a repeated name; it is an error below
  ids <- aln$id
  letters <- toupper(unname(aln$ali))
  n_structures <- length(ids)
  n_columns <- ncol(letters)

  repeated <- unique(ids[duplicated(ids)])
  if (length(repeated) > 0) {
    stop("The alignment has more than one sequence named: ", paste(repeated, collapse = ", "), call. = FALSE)
  }
  if (!is.data.frame(files) || !all(c("id", "file") %in% names(files))) {
    stop("`files` must be a data frame with columns `id` and `file`.", call. = FALSE)
  }
  table_ids <- as.character(files$id)
  repeated <- unique(table_ids[duplicated(table_ids)])
  if (length(repeated) > 0) {
    stop("`files` has more than one row for: ", paste(repeated, collapse = ", "), call. = FALSE)
  }
  without_file <- setdiff(ids, table_ids)
  if (length(without_file) > 0) {
    stop("`files` has no row for these sequences of the alignment: ", paste(without_file, collapse = ", "), call. = FALSE)
  }
  not_aligned <- setdiff(table_ids, ids)
  if (length(not_aligned) > 0) {
    stop("`files` has rows for ids that are not in the alignment: ", paste(not_aligned, collapse = ", "), call. = FALSE)
  }
  # the file of each sequence, in the order of the alignment
  files <- as.character(files$file)[match(ids, table_ids)]
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    stop("Files not found: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  resno <- matrix(NA_integer_, n_structures, n_columns)
  insert <- matrix(NA_character_, n_structures, n_columns)
  resname <- matrix(NA_character_, n_structures, n_columns)
  # coordinates: structure, x/y/z, column
  ca_xyz <- array(NA_real_, c(n_structures, 3, n_columns))
  n_xyz <- ca_xyz
  c_xyz <- ca_xyz
  o_xyz <- ca_xyz

  for (k in seq_len(n_structures)) {
    residues <- read_residues(files[k])
    columns <- which(letters[k, ] != "-")
    if (!identical(residues$letter, letters[k, columns])) {
      stop(sprintf(
        "%s: its residues do not match sequence %s of the alignment. File: %s. Alignment: %s.",
        files[k], ids[k], paste(residues$letter, collapse = ""), paste(letters[k, columns], collapse = "")
      ), call. = FALSE)
    }
    resno[k, columns] <- residues$resno
    insert[k, columns] <- residues$insert
    resname[k, columns] <- residues$resname
    ca_xyz[k, , columns] <- residues$ca
    n_xyz[k, , columns] <- residues$n
    c_xyz[k, , columns] <- residues$c
    o_xyz[k, , columns] <- residues$o
  }

  structure(
    list(
      id = ids, file = files, alignment = letters, resno = resno, insert = insert, resname = resname,
      ca = ca_xyz, n = n_xyz, c = c_xyz, o = o_xyz,
      properties = place_properties(properties, ids, resno, insert)
    ),
    class = "divpatterns_family"
  )
}

# The properties of read_family(): each property column of the table `properties` as a matrix like `resno`
# (structures x alignment columns), each value in the cell of its residue. An empty list for NULL.
place_properties <- function(properties, ids, resno, insert) {
  if (is.null(properties)) {
    return(list())
  }
  if (!is.data.frame(properties) || !all(c("id", "resno", "insert") %in% names(properties))) {
    stop("`properties` must be a data frame with columns `id`, `resno` and `insert`, and one column per property.", call. = FALSE)
  }
  property_names <- setdiff(names(properties), c("id", "resno", "insert"))
  if (length(property_names) == 0) {
    stop("`properties` has no property columns.", call. = FALSE)
  }
  for (name in property_names) {
    if (!is.numeric(properties[[name]])) {
      stop(sprintf("Property `%s` must be numeric.", name), call. = FALSE)
    }
    if (all(is.na(properties[[name]]))) {
      stop(sprintf("Property `%s` has no values.", name), call. = FALSE)
    }
  }
  if (!is.numeric(properties$resno) || anyNA(properties$resno) || any(properties$resno != round(properties$resno))) {
    stop("`resno` in `properties` must be whole numbers.", call. = FALSE)
  }

  # each row's residue, and each residue of the family, as "id resno" with the insertion code after resno
  code <- as.character(properties$insert)
  code[is.na(code) | trimws(code) == ""] <- ""
  row_residue <- paste0(properties$id, " ", as.integer(properties$resno), code) # "1abc_A 10", "1abc_A 10A"
  repeated <- unique(row_residue[duplicated(row_residue)])
  if (length(repeated) > 0) {
    stop("`properties` has more than one row for: ", paste(utils::head(repeated, 5), collapse = ", "), call. = FALSE)
  }
  cells <- which(!is.na(resno), arr.ind = TRUE) # row: structure, col: alignment column
  cell_insert <- insert[cells]
  cell_insert[is.na(cell_insert)] <- ""
  cell_residue <- paste0(ids[cells[, 1]], " ", resno[cells], cell_insert)
  where <- match(row_residue, cell_residue)
  if (anyNA(where)) {
    stop("`properties` has rows for residues that are not in the family: ", paste(utils::head(row_residue[is.na(where)], 5), collapse = ", "), call. = FALSE)
  }

  placed <- list()
  for (name in property_names) {
    values <- matrix(NA_real_, nrow(resno), ncol(resno))
    values[cells[where, , drop = FALSE]] <- properties[[name]]
    placed[[name]] <- values
  }
  placed
}

#' @export
print.divpatterns_family <- function(x, ...) {
  cat(sprintf("<divpatterns_family> %d structures, %d alignment columns\n", length(x$id), ncol(x$alignment)))
  invisible(x)
}

# The residues of a PDB file, in the order of the file: the amino acids with a Calpha atom. Returns their
# numbers, insertion codes, names and one-letter codes, and the coordinates of their Calpha, N, C and O
# atoms as matrices with one column per residue (NA where the file lacks an atom).
read_residues <- function(file) {
  # read.pdb() prints a note where a file has alternate locations (it keeps the first); not shown here
  utils::capture.output(pdb <- bio3d::read.pdb(file, verbose = FALSE))
  atoms <- pdb$atom
  ca_rows <- bio3d::atom.select(pdb, "calpha", verbose = FALSE)$atom
  if (length(ca_rows) == 0) {
    stop(sprintf("%s: no Calpha atoms.", file), call. = FALSE)
  }
  if (length(unique(atoms$chain[ca_rows])) > 1) {
    stop(sprintf("%s: the file must hold one chain.", file), call. = FALSE)
  }

  # an atom belongs to the residue with the same chain, number, insertion code and name
  residue_of_atom <- paste(atoms$chain, atoms$resno, atoms$insert, atoms$resid)
  residues <- residue_of_atom[ca_rows]
  if (anyDuplicated(residues)) {
    stop(sprintf("%s: a residue has more than one Calpha atom.", file), call. = FALSE)
  }

  coordinates <- list()
  for (atom_name in c("CA", "N", "C", "O")) {
    rows <- which(atoms$elety == atom_name & residue_of_atom %in% residues)
    if (anyDuplicated(residue_of_atom[rows])) {
      stop(sprintf("%s: a residue has more than one %s atom.", file, atom_name), call. = FALSE)
    }
    row_of_residue <- rows[match(residues, residue_of_atom[rows])]
    coordinates[[atom_name]] <- t(as.matrix(atoms[row_of_residue, c("x", "y", "z")]))
  }

  list(
    resno = atoms$resno[ca_rows],
    insert = atoms$insert[ca_rows],
    resname = atoms$resid[ca_rows],
    letter = unname(bio3d::aa321(atoms$resid[ca_rows])),
    ca = unname(coordinates$CA),
    n = unname(coordinates$N),
    c = unname(coordinates$C),
    o = unname(coordinates$O)
  )
}
