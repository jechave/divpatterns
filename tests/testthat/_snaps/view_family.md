# bad arguments are errors

    Code
      view_family(list(), open = FALSE)
    Condition
      Error:
      ! `fit` must come from `fit_to_reference()`.

---

    Code
      view_family(fit, title = c("a", "b"), open = FALSE)
    Condition
      Error:
      ! `title` must be one string.

