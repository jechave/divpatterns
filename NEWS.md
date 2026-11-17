# divpatterns 0.1.0

* First version: view a family of homologous structures in the browser.
* `fit_to_reference()` superposes each homolog on a reference on its own pair
  core, and gives each homolog's distance to the reference at every alignment
  column.
* `read_family()` reads a sequence alignment and one structure file per
  sequence, given by a table of ids and files, and places every residue in its
  alignment column.
* `view_family()` writes an interactive page, one file that works offline: the
  matrix of distances to the reference, linked to the 3D superposition and to
  the alignment.
