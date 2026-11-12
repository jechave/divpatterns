# divpatterns

<!-- badges: start -->
[![R-CMD-check](https://github.com/jechave/divpatterns/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jechave/divpatterns/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

divpatterns is a library of functions to calculate and analyse patterns of
protein evolutionary divergence at the sequence, structure and motion levels,
from homologous proteins, given their structures and a sequence alignment.

The package is in early development. Version 0.1.0 shows a family of
homologous structures in the browser.

## Installation

``` r
# install.packages("pak")
pak::pak("jechave/divpatterns")
```

## Viewing a family

You need a sequence alignment (FASTA) and one structure file (PDB format, one
chain) per aligned sequence. Each sequence, without its gaps, must be the
sequence of the residues in its file.

``` r
library(divpatterns)

family <- read_family("family.fasta", c("1abc_A.pdb", "2def_B.pdb", "3ghi_A.pdb"))
fit <- fit_to_reference(family, reference = "1abc_A")
view_family(fit, title = "My family")
```

`read_family()` places every residue of every structure in its alignment
column.

`fit_to_reference()` superposes each homolog on the reference on its own pair
core (residues within 1.2 Å after pruning) and gives each homolog's distance to
the reference at every alignment column.

`view_family()` writes one HTML file and opens it. The matrix of distances to
the reference (homologs by alignment columns) is linked to the 3D superposition
and to the alignment of the chains shown. Point at a cell to mark that residue
in 3D; click a row to show or hide that homolog.

The file holds its data and the 3Dmol.js library, so it works offline and can
be sent as it is. Cartoons need N, C and O atoms; chains with Calpha only are
drawn as traces.
