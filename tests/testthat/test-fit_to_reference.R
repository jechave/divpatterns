test_that("a rigid copy superposes exactly, on every residue", {
  ref <- helix(40)
  copy <- rotation(0.4, -1.2, 2) %*% ref + c(10, -4, 7)
  fit <- fit_to_reference(as_family(list(ref, copy)))
  expect_lt(max(fit$d2[2, ]), 1e-20)
  expect_equal(fit$core[[2]], 1:40)
  expect_equal(fit$reached[2], TRUE)
  expect_equal(fit$family$ca[2, , ], ref, tolerance = 1e-10)
  expect_identical(fit$family$ca[1, , ], ref)
})

test_that("the reference has no distances, core or flag", {
  ref <- helix(40)
  fit <- fit_to_reference(as_family(list(ref, ref + 1, ref - 1)), reference = 2)
  expect_equal(is.na(fit$d2[2, ]), rep(TRUE, 40))
  expect_null(fit$core[[2]])
  expect_equal(fit$reached, c(TRUE, NA, TRUE))
})

test_that("a displaced segment is left out of the core and keeps its displacement", {
  ref <- helix(40)
  moved <- ref
  moved[1, 36:40] <- moved[1, 36:40] + 5
  homolog <- rotation(1, 0.3, -0.5) %*% moved + c(-2, 3, 1)
  fit <- fit_to_reference(as_family(list(ref, homolog)))
  expect_equal(fit$core[[2]], 1:35)
  expect_equal(fit$d2[2, 1:35], rep(0, 35), tolerance = 1e-10)
  expect_equal(fit$d2[2, 36:40], rep(25, 5), tolerance = 1e-10)
})

test_that("distances are NA exactly where either structure has a gap", {
  ref <- helix(40)
  homolog <- rotation(0.2, 0.2, 0.2) %*% ref
  homolog[, 5:7] <- NA
  ref[, 39:40] <- NA
  fit <- fit_to_reference(as_family(list(ref, homolog)))
  expect_equal(which(is.na(fit$d2[2, ])), c(5:7, 39:40))
})

test_that("the reference can be given by identifier or by any position", {
  ref <- helix(40)
  moved <- ref
  moved[2, 1:8] <- moved[2, 1:8] - 3
  homolog <- rotation(0.5, 0.1, 1.3) %*% moved
  a <- fit_to_reference(as_family(list(ref, homolog), c("ref", "hom")))
  b <- fit_to_reference(as_family(list(homolog, ref), c("hom", "ref")), reference = "ref")
  expect_equal(b$reference, 2)
  expect_equal(b$d2[1, ], a$d2[2, ])
  expect_equal(b$core[[1]], a$core[[2]])
})

test_that("when no core reaches the cutoff, the smallest core is kept and flagged", {
  set.seed(1)
  ref <- helix(40)
  noisy <- ref + matrix(rnorm(120, sd = 3), 3)
  fit <- fit_to_reference(as_family(list(ref, noisy)), min_core = 15)
  expect_length(fit$core[[2]], 15)
  expect_equal(fit$reached[2], FALSE)
})

test_that("backbone atoms move with their Calpha", {
  ref <- helix(20)
  r <- rotation(0.9, -0.4, 0.6)
  homolog <- r %*% ref + c(3, 3, 3)
  family <- as_family(list(ref, homolog))
  n_ref <- ref + c(0, 0, -1.46)
  family$n[1, , ] <- n_ref
  family$n[2, , ] <- r %*% n_ref + c(3, 3, 3)
  fit <- fit_to_reference(family)
  expect_equal(fit$family$n[2, , ], n_ref, tolerance = 1e-10)
  expect_equal(is.na(fit$family$c), is.na(family$c))
})

test_that("homologs sharing too few residues with the reference are an error", {
  ref <- helix(40)
  short <- ref
  short[, 11:40] <- NA
  expect_snapshot(fit_to_reference(as_family(list(ref, ref, short), c("ref", "full", "short"))), error = TRUE)
})

test_that("bad arguments are errors", {
  family <- as_family(list(helix(20), helix(20)))
  expect_snapshot(fit_to_reference(list(ca = 1)), error = TRUE)
  expect_snapshot(fit_to_reference(family, reference = "nope"), error = TRUE)
  expect_snapshot(fit_to_reference(family, reference = 3), error = TRUE)
  expect_snapshot(fit_to_reference(family, cutoff = -1), error = TRUE)
})
