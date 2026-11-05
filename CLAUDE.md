# CLAUDE.md

## What this package is

`divpatterns` calculates observed patterns of divergence for real homologous
proteins, given their structures and a sequence alignment. These are the observed
patterns to which simulated evolution models are fitted and compared. Simulating
evolution is outside this package.

## Scope

ENMs are built with penm. Comparing proteins is done here, because alignments,
indels and averaging over pairs have no counterpart in penm.

Topics, in rough order:

- site metrics
- mode metrics
- sequence metrics, in a future version

## What to calculate

`dev/metrics_to_calculate.md` is the list of metrics.

## Comparing two proteins

Never trim an ENM. For mode comparisons, using the block of the covariance
matrix for the aligned residues without gaps is an idea under exploration, not
a decision.

## Comparing a family

Undecided. It may be the common core, or an average over pairs.
