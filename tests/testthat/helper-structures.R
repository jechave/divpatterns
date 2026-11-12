# A helical Calpha trace, consecutive Calpha 3.8 Angstrom apart (3 x n).
helix <- function(n) {
  i <- seq_len(n) - 1
  rbind(2.3 * cos(100 * pi / 180 * i), 2.3 * sin(100 * pi / 180 * i), 1.5 * i)
}

# A proper rotation from three angles (radians).
rotation <- function(a, b, c) {
  rx <- matrix(c(1, 0, 0, 0, cos(a), sin(a), 0, -sin(a), cos(a)), 3)
  ry <- matrix(c(cos(b), 0, -sin(b), 0, 1, 0, sin(b), 0, cos(b)), 3)
  rz <- matrix(c(cos(c), sin(c), 0, -sin(c), cos(c), 0, 0, 0, 1), 3)
  rz %*% ry %*% rx
}

# A family made directly from Calpha coordinates, without files: one 3 x columns matrix per structure, NA
# columns being gaps. All residues are alanines; there are no backbone atoms.
as_family <- function(coords, ids = paste0("s", seq_along(coords))) {
  n_columns <- ncol(coords[[1]])
  ca <- array(NA_real_, c(length(coords), 3, n_columns))
  for (k in seq_along(coords)) ca[k, , ] <- coords[[k]]
  present <- !is.na(ca[, 1, , drop = FALSE])[, 1, ]
  structure(
    list(
      id = ids, file = rep(NA_character_, length(ids)),
      alignment = ifelse(present, "A", "-"),
      resno = ifelse(present, col(present), NA_integer_),
      insert = matrix(NA_character_, length(ids), n_columns),
      resname = ifelse(present, "ALA", NA_character_),
      ca = ca, n = ca * NA, c = ca * NA, o = ca * NA
    ),
    class = "divpatterns_family"
  )
}

# A PDB file of one chain: Calpha at `ca` (3 x n), residues named from the one-letter `seq`, numbered
# `resno` with insertion codes `insert`, and, with `backbone`, N, C and O at bond lengths from their Calpha.
write_pdb <- function(path, ca, seq, backbone = TRUE, resno = seq_len(ncol(ca)), insert = rep(" ", ncol(ca)), chain = "A") {
  res <- toupper(bio3d::aa123(strsplit(seq, "")[[1]]))
  stopifnot(length(res) == ncol(ca))
  chain <- rep_len(chain, ncol(ca))
  lines <- character()
  serial <- 0
  for (i in seq_len(ncol(ca))) {
    atoms <- list(N = ca[, i] + c(0, 0, -1.46), CA = ca[, i], C = ca[, i] + c(1.52, 0, 0), O = ca[, i] + c(1.52, 1.23, 0))
    if (!backbone) atoms <- atoms["CA"]
    for (name in names(atoms)) {
      serial <- serial + 1
      p <- atoms[[name]]
      lines <- c(lines, sprintf(
        "ATOM  %5d %-4s %3s %s%4d%s   %8.3f%8.3f%8.3f  1.00 90.00          %2s",
        serial, paste0(" ", name), res[i], chain[i], resno[i], insert[i], p[1], p[2], p[3], substr(name, 1, 1)
      ))
    }
  }
  writeLines(c(lines, "END"), path)
}

# A family of three in `dir`: the reference (30 residues), a rigid copy missing residues 11-13, and a
# rigid copy, Calpha only, with two residues inserted after residue 20 and residues 26-30 moved 4 A.
# Returns the alignment file and the structure files.
toy_family <- function(dir) {
  seq <- "MKTAYIAKQRQISFVKSHFSRQLEERLGLI"
  ref <- helix(30)
  move <- function(x, a, t) rotation(a, a / 2, -a) %*% x + t
  h1 <- move(ref[, -(11:13)], 0.7, c(5, -3, 2))
  h2 <- ref
  h2[1, 26:30] <- h2[1, 26:30] + 4
  h2 <- cbind(h2[, 1:20], ref[, 20] + c(3, 3, 0), ref[, 20] + c(5, 6, 0), h2[, 21:30])
  h2 <- move(h2, -1.1, c(-7, 1, 4))
  files <- file.path(dir, c("ref.pdb", "h1.pdb", "h2.pdb"))
  write_pdb(files[1], ref, seq)
  write_pdb(files[2], h1, paste0(substr(seq, 1, 10), substr(seq, 14, 30)))
  write_pdb(files[3], h2, paste0(substr(seq, 1, 20), "GG", substr(seq, 21, 30)), backbone = FALSE)
  aln <- file.path(dir, "family.fasta")
  writeLines(c(
    ">ref", paste0(substr(seq, 1, 20), "--", substr(seq, 21, 30)),
    ">h1", paste0(substr(seq, 1, 10), "---", substr(seq, 14, 20), "--", substr(seq, 21, 30)),
    ">h2", paste0(substr(seq, 1, 20), "GG", substr(seq, 21, 30))
  ), aln)
  list(alignment = aln, files = files)
}

# A JSON array with nulls as a numeric vector with NA.
num <- function(x) vapply(x, \(v) if (is.null(v)) NA_real_ else v, 1)

# The data written into a page by view_family().
read_page_data <- function(file) {
  lines <- readLines(file, encoding = "UTF-8")
  i <- grep('<script id="data" type="application/json">', lines, fixed = TRUE)
  jsonlite::fromJSON(lines[i + 1], simplifyVector = FALSE)
}
