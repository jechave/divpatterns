# CLAUDE.md

## What this package is

`divpatterns` calculates observed patterns of divergence for real homologous
proteins, given their structures and a sequence alignment. These are the observed
patterns to which simulated evolution models are fitted and compared. Simulating
evolution is outside this package.

## Plan

Purpose: calculate observed divergence patterns from real homologous proteins,
for fitting and comparing simulated evolution models.

Principle: ENMs are built with penm. Comparing proteins is done here, because
alignments, indels and averaging over pairs have no counterpart in penm.
Integrating the two packages can be reconsidered later.

1. Site metrics, two proteins: structure and motion comparisons that use no
   normal modes.
2. Site metrics, family: explore how to combine proteins with different indels,
   then settle it. This includes site `lrmsd` and `nlrmsd`.
3. Mode metrics: explore how to compare two proteins, then a family. One idea
   to test is using the eigenvectors of the conserved covariance block.
   This includes mode `lrmsd` and `nlrmsd`.
4. Sequence metrics: a future version.

Each exploration step ends with a decision, recorded here, before its
functions become part of the package.

## What to calculate

`dev/metrics_to_calculate.md` is the list of metrics.

## Comparing two proteins

Never trim an ENM. For mode comparisons, using the block of the covariance
matrix for the aligned residues without gaps is an idea under exploration, not
a decision.

## Comparing a family

Undecided. It may be the common core, or an average over pairs.
