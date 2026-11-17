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

