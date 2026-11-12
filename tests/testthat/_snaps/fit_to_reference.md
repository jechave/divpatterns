# homologs sharing too few residues with the reference are an error

    Code
      fit_to_reference(as_family(list(ref, ref, short), c("ref", "full", "short")))
    Condition
      Error:
      ! short shares 10 residues with the reference; more than 15 are needed.

# bad arguments are errors

    Code
      fit_to_reference(list(ca = 1))
    Condition
      Error:
      ! `family` must come from `read_family()`.

---

    Code
      fit_to_reference(family, reference = "nope")
    Condition
      Error:
      ! `reference` must be one structure of the family, by position or identifier.

---

    Code
      fit_to_reference(family, reference = 3)
    Condition
      Error:
      ! `reference` must be one structure of the family, by position or identifier.

---

    Code
      fit_to_reference(family, cutoff = -1)
    Condition
      Error:
      ! `cutoff` must be a positive number.

