# tPC - Causal discovery with temporal background knowledge

This package implements the tPC algorithm for causal discovery. The 't' stands for 'temporal' or 'tiers' and
indicates that background knowledge in the form of a partial node/variable ordering is available.
Our implementation is a modified version of `pc` from the `pcalg` package (Kalisch et al. 2012) with the following
additional options:

Using the `tiers` argument, the user can allocate each node/variable to a tier. Specifying tiers has two effects:
First, conditional independence testing is restricted such that the variables in the conditioning set do not
lie in the future of the variables whose independence is being tested. This reduces the number of
unnecessary conditional independence tests and thus makes the algorithm more reliable. Second, edges between
nodes in different tiers are oriented from the earlier tier to the later tier. This usually results in a more
informative output. Both modifications were suggested in Spirtes et al. (2000), p. 93.

Additionally, further directed edges may be blacklisted using the `forbEdges` argument. In contrast to `pcalg`,
this allows the user to forbid one direction of an edge, but allow the other one. The arguments `context.all` and
`context.tier` function as whitelists. Variables in `context.all` are glocal context variables; as such, they are
parents of all other non-context nodes in the graphs (examples are variables encoding batch effect in gene expression data, or 'sex'
and 'country' in a cohort study). Variables in `context.tier` are tier-specific context variables,
which are parents of all non-context nodes in the same tier (e.g. 'calender year' if the tiers encode different years).

The package also includes a function called `ida_invalid`, which determines possibly valid adjustment sets from a graph
that is not a valid CPDAG or MPDAG.

## Install
To install and load this package in R from GitHub, make sure that the `devtools` package is installed and run the following commands:

```R
devtools::install_github("bips-hb/tpc")
library(tpc)
```

## References
  Markus Kalisch, Martin Maechler, Diego Colombo, Marloes H. Maathuis, Peter Buehlmann (2012). Causal
  Inference Using Graphical Models with the R Package pcalg. Journal of Statistical Software, 47(11), 1-26.
  URL https://www.jstatsoft.org/article/view/v047i11.
  
  Peter Spirtes, Clark Glymour, Richard Scheines (2000). Causation, Prediction, and Search. Second Edition. MIT Press, Cambridge, Massachusetts, USA.
