test_that("residues are placed in their alignment columns", {
  toy <- toy_family(withr::local_tempdir())
  family <- read_family(toy$alignment, toy$files)
  expect_equal(family$id, c("ref", "h1", "h2"))
  expect_equal(family$file, toy$files)
  expect_equal(which(is.na(family$resno[1, ])), 21:22)
  expect_equal(which(is.na(family$resno[2, ])), c(11:13, 21:22))
  expect_equal(family$resno[3, ], 1:32)
  expect_equal(family$resno[2, 14], 11L)
  expect_equal(family$resname[2, 14], "PHE")
  expect_equal(paste(family$alignment[2, ], collapse = ""), "MKTAYIAKQR---FVKSHFS--RQLEERLGLI")
  expect_equal(family$ca[1, , c(1:20, 23:32)], helix(30), tolerance = 1e-3)
  expect_equal(is.na(family$ca[2, 1, ]), is.na(family$resno[2, ]))
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
  family <- read_family(aln, files)
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
  expect_snapshot(read_family(toy$alignment, toy$files[1:2]), error = TRUE, transform = \(x) gsub(dir, "DIR", x, fixed = TRUE))
  expect_snapshot(read_family(toy$alignment, c(toy$files[1:2], "no/such/file.pdb")), error = TRUE)
  expect_snapshot(read_family(toy$alignment, toy$files[c(1, 3, 2)]), error = TRUE, transform = \(x) gsub(dir, "DIR", x, fixed = TRUE))
})

test_that("a file with two chains is an error", {
  dir <- withr::local_tempdir()
  seq <- "MKTAYIAKQRQISFVKSHFS"
  file <- file.path(dir, "two.pdb")
  write_pdb(file, helix(20), seq, chain = rep(c("A", "B"), each = 10))
  aln <- file.path(dir, "a.fasta")
  writeLines(c(">a", seq), aln)
  expect_snapshot(read_family(aln, file), error = TRUE, transform = \(x) gsub(dir, "DIR", x, fixed = TRUE))
})
