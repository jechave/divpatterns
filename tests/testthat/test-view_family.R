test_that("the page holds the structures, their distances and the row order", {
  dir <- withr::local_tempdir()
  toy <- toy_family(dir)
  fit <- fit_to_reference(read_family(toy$alignment, toy$files))
  file <- view_family(fit, file = file.path(dir, "page.html"), open = FALSE)
  expect_equal(file, file.path(dir, "page.html"))
  d <- read_page_data(file)
  expect_equal(d$n_columns, 32)
  expect_equal(d$reference, 0)
  expect_equal(unlist(d$colours)[c(1, 128, 256)], c("#440154", "#21908D", "#FDE725"))
  expect_setequal(unlist(d$order), 1:2)
  expect_equal(vapply(d$structures, \(s) s$id, ""), c("ref", "h1", "h2"))

  ref <- d$structures[[1]]
  expect_null(ref$core)
  expect_equal(unlist(ref$letter), fit$family$alignment[1, ])
  expect_equal(num(ref$resno), fit$family$resno[1, ])

  h2 <- d$structures[[3]]
  expect_equal(num(h2$distance), round(sqrt(fit$d2[3, ]), 3))
  expect_equal(num(h2$distance)[28:32], rep(4, 5), tolerance = 1e-3)
  expect_null(h2$distance[[21]])
  expect_equal(h2$core, length(fit$core[[3]]))
  expect_equal(h2$shared, 30)
  expect_equal(sum(!vapply(h2$n, is.null, TRUE)), 0)
  expect_equal(num(h2$ca[[32]]), round(fit$family$ca[3, , 32], 3))

  h1 <- d$structures[[2]]
  expect_equal(sum(!vapply(h1$n, is.null, TRUE)), 27)
  expect_null(h1$ca[[12]])
})

test_that("rows are ordered by clustering, the lower mean lrmsd first", {
  # rows 1 and 2 merge first (closest), rows 3 and 4 next; the pair with the lower mean goes first
  lrmsd <- rbind(rep(1, 20), rep(1.05, 20), rep(-1, 20), rep(-0.5, 20))
  expect_equal(cluster_order(lrmsd, log(4)), c(3, 4, 1, 2))
})

test_that("the page needs nothing from outside", {
  dir <- withr::local_tempdir()
  toy <- toy_family(dir)
  fit <- fit_to_reference(read_family(toy$alignment, toy$files))
  page <- readLines(view_family(fit, file = file.path(dir, "p.html"), title = "A <b> & c", open = FALSE), encoding = "UTF-8")
  lib <- readLines(system.file("viewer", "3Dmol-min.js", package = "divpatterns"), encoding = "UTF-8", warn = FALSE)
  expect_equal(setdiff(lib, page), character())
  expect_equal(grep("<script src|<link", page, value = TRUE), character())
  expect_equal(grep("__DATA__|__3DMOL__|__TITLE__", page, value = TRUE), character())
  expect_length(grep("<h1>A &lt;b&gt; &amp; c</h1>", page, fixed = TRUE), 1)
})

test_that("bad arguments are errors", {
  dir <- withr::local_tempdir()
  toy <- toy_family(dir)
  fit <- fit_to_reference(read_family(toy$alignment, toy$files))
  expect_snapshot(view_family(list(), open = FALSE), error = TRUE)
  expect_snapshot(view_family(fit, title = c("a", "b"), open = FALSE), error = TRUE)
})
