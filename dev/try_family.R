# Try the divpatterns package on one family: one page for each of the three sets in dev/data (built by
# dev/prepare/02_prepare.R). In each folder, alignment.fasta has the reference as its first sequence, and each
# sequence's structure is in <id>.pdb; the table `files` says so to read_family(). The two sets of AlphaFold
# models also have their pLDDT in properties.csv (written by dev/prepare/03_plddt.R).
#
# Install the package first (R CMD INSTALL . in the package folder), then run this from the package folder,
# line by line or with source(). The pages are written to temporary files and opened in the browser.

library(divpatterns)

family <- "362"

## Crystals ------------------------------------------------------------------------------------------------

folder <- file.path("dev/data/crystals", family)
ids <- bio3d::read.fasta(file.path(folder, "alignment.fasta"))$id
files <- data.frame(id = ids, file = file.path(folder, paste0(ids, ".pdb")))
crystals <- read_family(file.path(folder, "alignment.fasta"), files)
fit_crystals <- fit_to_reference(crystals)
view_family(fit_crystals, title = paste("Family", family, "- crystals"))

## Models of the crystals ------------------------------------------------------------------------------------

folder <- file.path("dev/data/crystal_models", family)
ids <- bio3d::read.fasta(file.path(folder, "alignment.fasta"))$id
files <- data.frame(id = ids, file = file.path(folder, paste0(ids, ".pdb")))
plddt <- read.csv(file.path(folder, "properties.csv"), colClasses = c(insert = "character"))
crystal_models <- read_family(file.path(folder, "alignment.fasta"), files, plddt)
fit_crystal_models <- fit_to_reference(crystal_models)
view_family(fit_crystal_models, title = paste("Family", family, "- models of the crystals"))

## UniProt models --------------------------------------------------------------------------------------------

folder <- file.path("dev/data/uniprot_models", family)
ids <- bio3d::read.fasta(file.path(folder, "alignment.fasta"))$id
files <- data.frame(id = ids, file = file.path(folder, paste0(ids, ".pdb")))
plddt <- read.csv(file.path(folder, "properties.csv"), colClasses = c(insert = "character"))
uniprot_models <- read_family(file.path(folder, "alignment.fasta"), files, plddt)
fit_uniprot_models <- fit_to_reference(uniprot_models)
view_family(fit_uniprot_models, title = paste("Family", family, "- UniProt models"))
