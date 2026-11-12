# files that do not match the alignment are errors

    Code
      read_family(toy$alignment, toy$files[1:2])
    Condition
      Error:
      ! The alignment has 3 sequences, so `files` must be 3 file paths.

---

    Code
      read_family(toy$alignment, c(toy$files[1:2], "no/such/file.pdb"))
    Condition
      Error:
      ! Files not found: no/such/file.pdb

---

    Code
      read_family(toy$alignment, toy$files[c(1, 3, 2)])
    Condition
      Error:
      ! DIR/h2.pdb: its residues do not match sequence h1 of the alignment. File: MKTAYIAKQRQISFVKSHFSGGRQLEERLGLI. Alignment: MKTAYIAKQRFVKSHFSRQLEERLGLI.

# a file with two chains is an error

    Code
      read_family(aln, file)
    Condition
      Error:
      ! DIR/two.pdb: the file must hold one chain.

