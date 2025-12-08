# QTL Hotspots

QTL hotspots, groups of traits co-mapping to the same genomic location,
are a common feature of genetical genomics studies.
Genomic locations associated with many traits are biologically interesting
since they may harbor influential regulators.
Nonetheless, non-genetic mechanisms, uncontrolled environmental factors
and unmeasured variables are capable of inducing a strong correlation structure
among clusters of transcripts, and as a consequence,
whenever a transcript shows a spurious linkage,
many correlated transcripts will likely map to the same locus,
creating a spurious QTL hotspot.
Permutation approaches that do not take into account
the phenotypic correlation tend to underestimate
the size of the hotspots that might appear by change in these situations
(Breitling et al. 2008).

This issue motivated the development of permutation tests that preserve
the correlation structure of the phenotypes in order to determine
the significance of QTL hotspots
(Breitling et al. 2008, Chaibub Neto et al. 2012).
In this tutorial we present software tools implementing the _NL_-method
(Chaibub Neto et al. 2012), the _N_-method (Breitling et al. 2008),
and the _Q_-method (West et al. 2007, Wu et al. 2008) permutation approaches.

## Installation

```
devtools::install_github("byandell/qtlhot", build_vignettes=TRUE)
```

## Issues

Modify to use names element of highlod object?
Fix so call with length(pheno2) == 1 and > 1 give same results; see JoinTestOutputs.


