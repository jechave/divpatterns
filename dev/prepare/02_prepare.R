# Build dev/data from data-raw. For every family of Mathilde's three sets, one folder
# dev/data/<set>/<family>/ with:
#
#   alignment.fasta   her alignment (<family>.fasta.ali), with her reference moved to the first row and its _REF
#                     suffix removed: the reference is the first sequence
#   <id>.pdb          one file per sequence, named by its id in alignment.fasta: one protein chain, all its atoms
#   mathilde_rmsd_core.tsv   her per-site table (<family>.RMSD_core), unchanged: one row per alignment column, with
#                     her site profile and the M-CSA reference numbering of each column. Its protein names keep
#                     her _REF suffix.
#
# and, for each set, dev/data/<set>/mathilde_families.tsv: her Families.tsv, unchanged.
#
# Sets and where their structures come from:
#
#   crystals         the reference: Mathilde's file, which already has all its atoms; the others: the chain cut
#                    out of its PDB entry, downloaded by 01_download.R (PDB format, or mmCIF for the entries too
#                    large for it)
#   crystal_models   Mathilde's files, which already have all their atoms
#   uniprot_models   the AlphaFold model downloaded by 01_download.R
#
# Check, for every family: read_family() on the new folder must give exactly the alignment and residues that it
# gives on Mathilde's files, and the same Calpha coordinates up to 0.001 Angstrom (the last decimal of a PDB
# file). A family that fails the check is reported and left out.
#
#   Rscript dev/prepare/02_prepare.R        (from the package folder, after 01_download.R)

devtools::load_all(quiet = TRUE)

source_dir <- "data-raw/afdata20261104"
rcsb_dir <- "data-raw/downloads/rcsb"
afdb_dir <- "data-raw/downloads/afdb"
out_dir <- "dev/data"
sets <- c(crystals = "EXP_ALL_MCSA", crystal_models = "AF_MDL_EXP_ALL_MCSA", uniprot_models = "UNIPROT_ALL_MCSA")

# Mathilde's file for one sequence (id as in her alignment, with _REF if it is her reference). Her crystal and
# UniProt files are named by id; her models of the crystals by UniProt accession ("P62942_1d6o_A_REF" is in
# AF-P62942-F1-model_v6.pdb).
mathilde_file <- function(set, family, id) {
  folder <- file.path(source_dir, sets[[set]], "data/pdb", family)
  if (set == "crystal_models") {
    accession <- sub("_.*$", "", id)
    file <- list.files(folder, paste0("^AF-", accession, "-F1[-.]"), full.names = TRUE)
  } else {
    file <- file.path(folder, paste0(id, ".pdb"))
  }
  if (length(file) != 1 || !file.exists(file)) stop("no single file of Mathilde's for ", id, " in ", folder)
  file
}

# Write the protein atoms of one chain of a downloaded PDB entry ("3ni6_B": chain B of entry 3ni6) to file:
# its ATOM records only, as in Mathilde's files, whose sequences (and so her alignments) leave out the residues
# written as HETATM (selenomethionine, modified residues, amino acids bound as ligands).
# bio3d keeps the first model and, at alternate locations, the first one. Its mmCIF reader warns that it is a
# beta version; the check against Mathilde's files is what tells whether the chain was read right.
write_crystal_chain <- function(id, file) {
  entry <- sub("_.*$", "", id)
  chain <- sub("^.*_", "", id)
  pdb_file <- file.path(rcsb_dir, paste0(entry, ".pdb"))
  cif_file <- file.path(rcsb_dir, paste0(entry, ".cif"))
  if (file.exists(pdb_file)) {
    utils::capture.output(structure <- bio3d::read.pdb(pdb_file, verbose = FALSE))
  } else if (file.exists(cif_file)) {
    utils::capture.output(structure <- suppressWarnings(bio3d::read.cif(cif_file, verbose = FALSE)))
  } else {
    stop("entry ", entry, " was not downloaded")
  }
  atoms <- bio3d::atom.select(structure, "protein", type = "ATOM", chain = chain, verbose = FALSE)
  if (length(atoms$atom) == 0) stop("entry ", entry, " has no protein chain ", chain)
  one_chain <- bio3d::trim.pdb(structure, atoms)
  # The PDB format has one character for the chain. Longer chain ids (mmCIF only, e.g. "AAA") are written as
  # their first character, as in Mathilde's files.
  one_chain$atom$chain <- substr(one_chain$atom$chain, 1, 1)
  bio3d::write.pdb(one_chain, file = file)
}

# Copy the downloaded AlphaFold model of one sequence ("AF-P62942-F1") to file.
copy_model <- function(id, file) {
  model <- list.files(afdb_dir, paste0("^", id, "-model_v[0-9]+[.]pdb$"), full.names = TRUE)
  if (length(model) != 1) stop("no single downloaded model for ", id)
  file.copy(model, file, overwrite = TRUE)
}

one_family <- function(set, alignment_file) {
  family <- sub("[.]fasta[.]ali$", "", basename(alignment_file))
  folder <- file.path(out_dir, set, family)
  unlink(folder, recursive = TRUE)
  dir.create(folder, recursive = TRUE)

  # the alignment, reference first
  alignment <- bio3d::read.fasta(alignment_file)
  ref <- grep("_REF$", alignment$id)
  if (length(ref) != 1) stop("family ", family, ": not exactly one sequence ends in _REF")
  rows <- c(ref, setdiff(seq_along(alignment$id), ref))
  original_ids <- alignment$id[rows]
  ids <- sub("_REF$", "", original_ids)
  bio3d::write.fasta(ids = ids, seqs = alignment$ali[rows, , drop = FALSE], file = file.path(folder, "alignment.fasta"))

  # the structures
  files <- file.path(folder, paste0(ids, ".pdb"))
  for (i in seq_along(ids)) {
    if (set == "crystals" && i == 1) file.copy(mathilde_file(set, family, original_ids[i]), files[i], overwrite = TRUE)
    if (set == "crystals" && i > 1) write_crystal_chain(ids[i], files[i])
    if (set == "crystal_models") file.copy(mathilde_file(set, family, original_ids[i]), files[i], overwrite = TRUE)
    if (set == "uniprot_models") copy_model(ids[i], files[i])
  }

  # her per-site table (family 2's file alone is named .RMSD_core.tsv)
  table <- Sys.glob(file.path(source_dir, sets[[set]], "data/rmsd", paste0(family, ".RMSD_core*")))
  if (length(table) != 1) stop("family ", family, ": not exactly one RMSD_core table")
  if (!file.copy(table, file.path(folder, "mathilde_rmsd_core.tsv"))) stop("family ", family, ": could not copy its RMSD_core table")

  # the check
  new <- read_family(file.path(folder, "alignment.fasta"), data.frame(id = ids, file = files))
  old_files <- character(length(original_ids))
  for (i in seq_along(original_ids)) old_files[i] <- mathilde_file(set, family, original_ids[i])
  old <- read_family(alignment_file, data.frame(id = original_ids, file = old_files))
  for (part in c("alignment", "resno", "insert", "resname")) {
    if (!identical(unname(new[[part]]), unname(old[[part]][rows, , drop = FALSE]))) stop("family ", family, ": ", part, " differs from Mathilde's files")
  }
  # Calpha coordinates: the same atoms, and the same values up to 0.001 Angstrom, the last decimal of a PDB file
  # (some downloaded PDB files round a few last digits differently from Mathilde's files). The limit is 0.0015
  # because the difference of two numbers with 3 decimals is not exactly 0.001 in floating point.
  old_ca <- old$ca[rows, , , drop = FALSE]
  if (!identical(is.na(unname(new$ca)), is.na(unname(old_ca)))) stop("family ", family, ": Calpha atoms differ from Mathilde's files")
  largest <- max(abs(new$ca - old_ca), na.rm = TRUE)
  if (largest > 0.0015) stop("family ", family, ": Calpha coordinates differ from Mathilde's files by up to ", largest, " Angstrom")
  data.frame(set = set, family = family, structures = length(ids), largest_difference = largest,
             with_backbone = sum(apply(!is.na(new$n[, 1, , drop = FALSE]), 1, any)))
}

done <- list()
failed <- character()
for (set in names(sets)) {
  alignments <- list.files(file.path(source_dir, sets[[set]], "data/sequences"), "^[0-9]+[.]fasta[.]ali$", full.names = TRUE)
  cat(set, ":", length(alignments), "families\n")
  dir.create(file.path(out_dir, set), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(file.path(source_dir, sets[[set]], "Families.tsv"), file.path(out_dir, set, "mathilde_families.tsv"), overwrite = TRUE)) stop("could not copy Families.tsv of ", set)
  for (alignment_file in alignments) {
    result <- tryCatch(one_family(set, alignment_file), error = function(e) e)
    if (inherits(result, "error")) {
      family <- sub("[.]fasta[.]ali$", "", basename(alignment_file))
      unlink(file.path(out_dir, set, family), recursive = TRUE)
      failed <- c(failed, paste0(set, " ", family, ": ", conditionMessage(result)))
    } else {
      done[[length(done) + 1]] <- result
    }
  }
}

summary <- do.call(rbind, done)
for (set in names(sets)) {
  s <- summary[summary$set == set, ]
  cat(sprintf("%s: %d families, %d structures, %d with backbone atoms; %d families with Calpha differences, largest %.4f Angstrom\n",
              set, nrow(s), sum(s$structures), sum(s$with_backbone), sum(s$largest_difference > 0), max(s$largest_difference)))
}
if (length(failed) > 0) {
  cat("\nleft out:\n")
  cat(failed, sep = "\n")
  quit(status = 1)
}
