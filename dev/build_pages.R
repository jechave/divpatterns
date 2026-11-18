# Build the viewer page of every family of every set in dev/data (built by dev/prepare/02_prepare.R):
# dev/pages/<set>/<family>.html, one self-contained HTML file each, with the reference as the first sequence
# of the alignment. Each structure is in <id>.pdb in the family's folder. The AlphaFold models also have
# properties.csv (pLDDT, written by dev/prepare/03_plddt.R), so their pages can be coloured by pLDDT; the
# crystals have none.
#
#   Rscript dev/build_pages.R        (from the package folder)

devtools::load_all(quiet = TRUE)
library(parallel)

set_names <- c(crystals = "crystals", crystal_models = "models of the crystals", uniprot_models = "UniProt models")

jobs <- data.frame(set = character(), family = character())
for (set in names(set_names)) {
  families <- list.files(file.path("dev/data", set), "^[0-9]+$")
  jobs <- rbind(jobs, data.frame(set = set, family = families))
  dir.create(file.path("dev/pages", set), recursive = TRUE, showWarnings = FALSE)
}

build_page <- function(set, family) {
  folder <- file.path("dev/data", set, family)
  alignment <- file.path(folder, "alignment.fasta")
  ids <- bio3d::read.fasta(alignment)$id
  files <- data.frame(id = ids, file = file.path(folder, paste0(ids, ".pdb")))
  properties <- NULL
  if (set != "crystals") properties <- utils::read.csv(file.path(folder, "properties.csv"), colClasses = c(insert = "character"))
  fit <- fit_to_reference(read_family(alignment, files, properties))
  view_family(fit, file = file.path("dev/pages", set, paste0(family, ".html")),
              title = paste0("Family ", family, ", ", set_names[[set]]), open = FALSE)
}

results <- mclapply(seq_len(nrow(jobs)), function(i) {
  tryCatch(build_page(jobs$set[i], jobs$family[i]), error = function(e) e)
}, mc.cores = detectCores(), mc.preschedule = FALSE)

failed <- vapply(results, inherits, TRUE, what = "error")
cat(sum(!failed), "pages written to dev/pages\n")
if (any(failed)) {
  for (i in which(failed)) cat("failed:", jobs$set[i], jobs$family[i], ":", conditionMessage(results[[i]]), "\n")
  quit(status = 1)
}
