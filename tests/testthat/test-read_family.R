test_that("residues are placed in their alignment columns", {
  toy <- toy_family(withr::local_tempdir())
  family <- read_family(toy$alignment, toy$files)
  expect_equal(family$id, c("ref", "h1", "h2"))
  expect_equal(family$file, toy$files$file)
  expect_equal(which(is.na(family$resno[1, ])), 21:22)
  expect_equal(which(is.na(family$resno[2, ])), c(11:13, 21:22))
  expect_equal(family$resno[3, ], 1:32)
  expect_equal(family$resno[2, 14], 11L)
  expect_equal(family$resname[2, 14], "PHE")
  expect_equal(paste(family$alignment[2, ], collapse = ""), "MKTAYIAKQR---FVKSHFS--RQLEERLGLI")
  expect_equal(family$ca[1, , c(1:20, 23:32)], helix(30), tolerance = 1e-3)
  expect_equal(is.na(family$ca[2, 1, ]), is.na(family$resno[2, ]))
})

test_that("the rows of `files` can be in any order", {
  toy <- toy_family(withr::local_tempdir())
  expect_identical(read_family(toy$alignment, toy$files[c(3, 1, 2), ]), read_family(toy$alignment, toy$files))
})

test_that("backbone atoms are read where files have them", {
  toy <- toy_family(withr::local_tempdir())
  family <- read_family(toy$alignment, toy$files)
  columns <- which(!is.na(family$resno[1, ]))
  n_to_ca <- sqrt(colSums((family$n[1, , columns] - family$ca[1, , columns])^2))
  c_to_ca <- sqrt(colSums((family$c[1, , columns] - family$ca[1, , columns])^2))
  expect_equal(n_to_ca, rep(1.46, 30), tolerance = 1e-3)
  expect_equal(c_to_ca, rep(1.52, 30), tolerance = 1e-3)
  expect_equal(sum(!is.na(family$o[3, , ])), 0)
})

test_that("residues with insertion codes are told apart", {
  dir <- withr::local_tempdir()
  seq <- "MKTAYIAKQLLLSFVKSHFS" # residues 10, 10A and 10B: all LEU
  files <- file.path(dir, c("a.pdb", "b.pdb"))
  write_pdb(files[1], helix(20), seq, resno = c(1:10, 10, 10, 11:18), insert = c(rep(" ", 10), "A", "B", rep(" ", 8)))
  write_pdb(files[2], helix(20) + 1, seq)
  aln <- file.path(dir, "a.fasta")
  writeLines(c(">a", seq, ">b", seq), aln)
  family <- read_family(aln, data.frame(id = c("a", "b"), file = files))
  expect_equal(family$resno[1, 10:12], c(10L, 10L, 10L))
  expect_equal(family$insert[1, 10:12], c(NA, "A", "B"))
  expect_equal(family$ca[1, , ], helix(20), tolerance = 1e-3)
  n_to_ca <- sqrt(colSums((family$n[1, , ] - family$ca[1, , ])^2))
  expect_equal(n_to_ca, rep(1.46, 20), tolerance = 1e-3)
})

test_that("a family read from files fits as expected", {
  toy <- toy_family(withr::local_tempdir())
  fit <- fit_to_reference(read_family(toy$alignment, toy$files))
  # h1 is a rigid copy, up to the 3 decimals of PDB coordinates
  expect_equal(sum(!is.na(fit$d2[2, ])), 27)
  expect_lt(max(fit$d2[2, ], na.rm = TRUE), 1e-5)
  # h2: residues 26-30 moved 4 A, outside its core
  expect_equal(fit$core[[3]], setdiff(which(!is.na(fit$d2[3, ])), 28:32))
  expect_equal(fit$d2[3, 28:32], rep(16, 5), tolerance = 1e-3)
})

test_that("files that do not match the alignment are errors", {
  dir <- withr::local_tempdir()
  toy <- toy_family(dir)
  files <- toy$files
  expect_snapshot(read_family(toy$alignment, files$file), error = TRUE)
  expect_snapshot(read_family(toy$alignment, files[1:2, ]), error = TRUE)
  expect_snapshot(read_family(toy$alignment, rbind(files, files[2, ])), error = TRUE)
  expect_snapshot(read_family(toy$alignment, rbind(files, data.frame(id = "x", file = files$file[1]))), error = TRUE)
  missing <- files
  missing$file[3] <- "no/such/file.pdb"
  expect_snapshot(read_family(toy$alignment, missing), error = TRUE)
  swapped <- files
  swapped$file <- files$file[c(1, 3, 2)]
  expect_snapshot(read_family(toy$alignment, swapped), error = TRUE, transform = \(x) gsub(dir, "DIR", x, fixed = TRUE))
  twice <- file.path(dir, "twice.fasta")
  writeLines(c(readLines(toy$alignment), ">h1", "MKTAYIAKQR---FVKSHFS--RQLEERLGLI"), twice)
  expect_snapshot(read_family(twice, files), error = TRUE)
})

test_that("a file with two chains is an error", {
  dir <- withr::local_tempdir()
  seq <- "MKTAYIAKQRQISFVKSHFS"
  file <- file.path(dir, "two.pdb")
  write_pdb(file, helix(20), seq, chain = rep(c("A", "B"), each = 10))
  aln <- file.path(dir, "a.fasta")
  writeLines(c(">a", seq), aln)
  expect_snapshot(read_family(aln, data.frame(id = "a", file = file)), error = TRUE, transform = \(x) gsub(dir, "DIR", x, fixed = TRUE))
})

test_that("properties are placed in the columns of their residues", {
  toy <- toy_family(withr::local_tempdir())
  properties <- data.frame(
    id = c("h2", "ref", "ref", "h2"),
    resno = c(32, 1, 30, 21),
    insert = c("", NA, NA, ""),
    plddt = c(41.5, 90, NA, 77),
    other = c(1, 2, 3, 4)
  )
  family <- read_family(toy$alignment, toy$files, properties)
  expect_equal(names(family$properties), c("plddt", "other"))
  plddt <- family$properties$plddt
  expect_equal(dim(plddt), c(3, 32))
  expect_equal(plddt[1, 1], 90)
  expect_equal(plddt[3, 32], 41.5)
  expect_equal(plddt[3, 21], 77) # h2's residue 21 is in column 21, a gap of the reference
  expect_equal(sum(!is.na(plddt)), 3)
  expect_equal(family$properties$other[1, 32], 3) # the reference's residue 30 is in column 32
  expect_equal(read_family(toy$alignment, toy$files)$properties, list())
})

test_that("properties match residues by insertion code", {
  dir <- withr::local_tempdir()
  seq <- "MKTAYIAKQLLLSFVKSHFS"
  file <- file.path(dir, "a.pdb")
  write_pdb(file, helix(20), seq, resno = c(1:10, 10, 10, 11:18), insert = c(rep(" ", 10), "A", "B", rep(" ", 8)))
  aln <- file.path(dir, "a.fasta")
  writeLines(c(">a", seq), aln)
  properties <- data.frame(id = "a", resno = 10, insert = c("B", "", "A"), x = c(3, 1, 2))
  family <- read_family(aln, data.frame(id = "a", file = file), properties)
  expect_equal(family$properties$x[1, 10:12], c(1, 2, 3))
})

test_that("property tables that do not match the family are errors", {
  toy <- toy_family(withr::local_tempdir())
  good <- data.frame(id = "ref", resno = 1:3, insert = NA, plddt = c(90, 80, 70))
  expect_snapshot(read_family(toy$alignment, toy$files, good[, c("id", "resno", "plddt")]), error = TRUE)
  expect_snapshot(read_family(toy$alignment, toy$files, good[, c("id", "resno", "insert")]), error = TRUE)
  expect_snapshot(read_family(toy$alignment, toy$files, transform(good, plddt = as.character(plddt))), error = TRUE)
  expect_snapshot(read_family(toy$alignment, toy$files, transform(good, plddt = NA_real_)), error = TRUE)
  expect_snapshot(read_family(toy$alignment, toy$files, transform(good, resno = c(1, 2.5, 3))), error = TRUE)
  expect_snapshot(read_family(toy$alignment, toy$files, rbind(good, good[2, ])), error = TRUE)
  expect_snapshot(read_family(toy$alignment, toy$files, rbind(good, data.frame(id = c("h1", "x"), resno = c(28, 1), insert = NA, plddt = 1))), error = TRUE)
})
