# files that do not match the alignment are errors

    Code
      read_family(toy$alignment, files$file)
    Condition
      Error:
      ! `files` must be a data frame with columns `id` and `file`.

---

    Code
      read_family(toy$alignment, files[1:2, ])
    Condition
      Error:
      ! `files` has no row for these sequences of the alignment: h2

---

    Code
      read_family(toy$alignment, rbind(files, files[2, ]))
    Condition
      Error:
      ! `files` has more than one row for: h1

---

    Code
      read_family(toy$alignment, rbind(files, data.frame(id = "x", file = files$file[
        1])))
    Condition
      Error:
      ! `files` has rows for ids that are not in the alignment: x

---

    Code
      read_family(toy$alignment, missing)
    Condition
      Error:
      ! Files not found: no/such/file.pdb

---

    Code
      read_family(toy$alignment, swapped)
    Condition
      Error:
      ! DIR/h2.pdb: its residues do not match sequence h1 of the alignment. File: MKTAYIAKQRQISFVKSHFSGGRQLEERLGLI. Alignment: MKTAYIAKQRFVKSHFSRQLEERLGLI.

---

    Code
      read_family(twice, files)
    Condition
      Error:
      ! The alignment has more than one sequence named: h1

# a file with two chains is an error

    Code
      read_family(aln, data.frame(id = "a", file = file))
    Condition
      Error:
      ! DIR/two.pdb: the file must hold one chain.

# property tables that do not match the family are errors

    Code
      read_family(toy$alignment, toy$files, good[, c("id", "resno", "plddt")])
    Condition
      Error:
      ! `properties` must be a data frame with columns `id`, `resno` and `insert`, and one column per property.

---

    Code
      read_family(toy$alignment, toy$files, good[, c("id", "resno", "insert")])
    Condition
      Error:
      ! `properties` has no property columns.

---

    Code
      read_family(toy$alignment, toy$files, transform(good, plddt = as.character(
        plddt)))
    Condition
      Error:
      ! Property `plddt` must be numeric.

---

    Code
      read_family(toy$alignment, toy$files, transform(good, plddt = NA_real_))
    Condition
      Error:
      ! Property `plddt` has no values.

---

    Code
      read_family(toy$alignment, toy$files, transform(good, resno = c(1, 2.5, 3)))
    Condition
      Error:
      ! `resno` in `properties` must be whole numbers.

---

    Code
      read_family(toy$alignment, toy$files, rbind(good, good[2, ]))
    Condition
      Error:
      ! `properties` has more than one row for: ref 2

---

    Code
      read_family(toy$alignment, toy$files, rbind(good, data.frame(id = c("h1", "x"),
      resno = c(28, 1), insert = NA, plddt = 1)))
    Condition
      Error:
      ! `properties` has rows for residues that are not in the family: h1 28, x 1

