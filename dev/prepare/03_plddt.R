# Write the pLDDT of every AlphaFold model in dev/data (built by 02_prepare.R) as a property table, for
# read_family(properties = ...): dev/data/<set>/<family>/properties.csv, for the two sets of AlphaFold models
# (crystal_models and uniprot_models), with columns id, resno, insert and plddt, one row per residue.
#
# AlphaFold writes each residue's pLDDT in the B-factor column of its atoms (the same value for all the atoms
# of a residue); the value of the Calpha atom is taken. The crystals have no pLDDT and get no table.
#
#   Rscript dev/prepare/03_plddt.R        (from the package folder, after 02_prepare.R)

library(parallel)

sets <- c("crystal_models", "uniprot_models")

one_family <- function(folder) {
  ids <- bio3d::read.fasta(file.path(folder, "alignment.fasta"))$id
  rows <- list()
  for (id in ids) {
    utils::capture.output(structure <- bio3d::read.pdb(file.path(folder, paste0(id, ".pdb")), verbose = FALSE))
    ca <- structure$atom[structure$atom$elety == "CA" & structure$atom$type == "ATOM", ]
    if (anyNA(ca$b) || any(ca$b < 0 | ca$b > 100)) stop(id, ": B-factor column is not a pLDDT (values outside 0-100)")
    rows[[id]] <- data.frame(id = id, resno = ca$resno, insert = ifelse(is.na(ca$insert), "", ca$insert), plddt = ca$b)
  }
  table <- do.call(rbind, rows)
  utils::write.csv(table, file.path(folder, "properties.csv"), row.names = FALSE)
  nrow(table)
}

folders <- unlist(lapply(sets, function(set) file.path("dev/data", set, list.files(file.path("dev/data", set), "^[0-9]+$"))))
results <- mclapply(folders, function(folder) tryCatch(one_family(folder), error = function(e) e), mc.cores = detectCores(), mc.preschedule = FALSE)
failed <- vapply(results, inherits, TRUE, what = "error")
cat(sum(!failed), "property tables written,", sum(unlist(results[!failed])), "residues\n")
if (any(failed)) {
  for (i in which(failed)) cat("failed:", folders[i], ":", conditionMessage(results[[i]]), "\n")
  quit(status = 1)
}
