# Download the full-atom structures of the proteins in Mathilde's alignments (data-raw/afdata20261104), as
# they are served, into data-raw/downloads:
#
#   rcsb/   one file per PDB entry of the crystal alignments: <entry>.pdb, or <entry>.cif for the entries
#           too large for the PDB format
#   afdb/   one AlphaFold model per UniProt accession of the UniProt alignments, under the name the AlphaFold
#           Database gives it (e.g. AF-P62942-F1-model_v6.pdb)
#
# Each folder has downloaded.tsv: file, address, date. Files already downloaded are skipped, so the script can
# be run again after an interruption. The models of the crystals need no download: Mathilde's files have all
# their atoms.
#
#   Rscript dev/prepare/01_download.R        (from the package folder)

source_dir <- "data-raw/afdata20261104"
rcsb_dir <- "data-raw/downloads/rcsb"
afdb_dir <- "data-raw/downloads/afdb"
dir.create(rcsb_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(afdb_dir, recursive = TRUE, showWarnings = FALSE)
headers <- c("User-Agent" = "divpatterns")

# The sequence ids of all the main alignments (<family>.fasta.ali) of one of Mathilde's folders, without the
# _REF suffix that marks her reference.
alignment_ids <- function(folder) {
  files <- list.files(file.path(source_dir, folder, "data/sequences"), "^[0-9]+[.]fasta[.]ali$", full.names = TRUE)
  ids <- character()
  for (f in files) ids <- c(ids, bio3d::read.fasta(f)$id)
  sub("_REF$", "", ids)
}

# Download address to file. Returns TRUE if it worked, FALSE if the server refused or failed (the file is then
# not written). Writes to a temporary name first, so that an interrupted download leaves no partial file.
fetch <- function(address, file) {
  part <- paste0(file, ".part")
  ok <- tryCatch(download.file(address, part, quiet = TRUE, headers = headers) == 0, error = function(e) FALSE, warning = function(w) FALSE)
  if (ok) {
    file.rename(part, file)
  } else {
    unlink(part)
  }
  ok
}

record <- function(folder, file, address) {
  log <- file.path(folder, "downloaded.tsv")
  if (!file.exists(log)) cat("file\taddress\tdate\n", file = log)
  cat(file, "\t", address, "\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = log, append = TRUE, sep = "")
}

failed <- character()

## Crystals: PDB entries --------------------------------------------------------------------------------------

crystal_ids <- alignment_ids("EXP_ALL_MCSA")                       # e.g. "1d6o_A", "3ni6_B"
entries <- unique(sub("_.*$", "", crystal_ids))                     # e.g. "1d6o", "3ni6"
cat(length(entries), "PDB entries\n")
for (entry in entries) {
  pdb_file <- file.path(rcsb_dir, paste0(entry, ".pdb"))
  cif_file <- file.path(rcsb_dir, paste0(entry, ".cif"))
  if (file.exists(pdb_file) || file.exists(cif_file)) next
  pdb_address <- sprintf("https://files.rcsb.org/download/%s.pdb", toupper(entry))
  cif_address <- sprintf("https://files.rcsb.org/download/%s.cif", toupper(entry))
  if (fetch(pdb_address, pdb_file)) {
    record(rcsb_dir, basename(pdb_file), pdb_address)
  } else if (fetch(cif_address, cif_file)) {
    record(rcsb_dir, basename(cif_file), cif_address)
  } else {
    failed <- c(failed, entry)
  }
}

## UniProt models: AlphaFold Database ---------------------------------------------------------------------------

model_ids <- unique(alignment_ids("UNIPROT_ALL_MCSA"))             # e.g. "AF-P62942-F1"
accessions <- sub("^AF-(.*)-F1$", "\\1", model_ids)                # e.g. "P62942"
stopifnot(paste0("AF-", accessions, "-F1") == model_ids)
already <- list.files(afdb_dir, "^AF-.*[.]pdb$")
cat(length(accessions), "AlphaFold models\n")
for (accession in accessions) {
  if (any(startsWith(already, paste0("AF-", accession, "-F1-")))) next
  # the AlphaFold Database's description of the model gives the address of its current PDB file
  api <- sprintf("https://alphafold.ebi.ac.uk/api/prediction/%s", accession)
  answer <- tryCatch(jsonlite::fromJSON(api), error = function(e) NULL)
  if (is.null(answer)) {
    failed <- c(failed, accession)
    next
  }
  address <- answer$pdbUrl[answer$entryId == paste0("AF-", accession, "-F1")]
  file <- file.path(afdb_dir, basename(address))
  if (length(address) == 1 && fetch(address, file)) {
    record(afdb_dir, basename(file), address)
  } else {
    failed <- c(failed, accession)
  }
}

cat("downloaded files:", length(list.files(rcsb_dir, "[.](pdb|cif)$")), "from the RCSB,", length(list.files(afdb_dir, "[.]pdb$")), "from the AlphaFold Database\n")
if (length(failed) > 0) {
  cat("failed (run again to retry):", failed, "\n")
  quit(status = 1)
}
