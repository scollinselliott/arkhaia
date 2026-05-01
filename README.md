
<!-- README.md is generated from README.Rmd. Please edit that file -->

# arkhaia: Archaeological and Historical Analysis

<!-- badges: start -->

<!-- badges: end -->

The R package `arkhaia` contains functions related to research on
economic relationships via archaeological or otherwise historical data.
The main focus is on evaluation of changes in long-term “integrating” or
“integrated” relationships over time between communities, primarily
through the material evidence of artifact assemblages. The problem of
how to detect similar patterns over time in a given behavior is a broad
one that requires specificity in order to address mathematically, and so
this package provides the tools necessary for handling counts of
artifacts, as well as measurement data (such as prices).

The package relies on `Rcpp` and `RcppArmadillo` (Eddelbuettel and
Sanderson 2014; Eddelbuettel and Balamuta 2018).

## Installation

To obtain the development version of the package, `devtools` can be used
to install `arkhaia` from GitHub:

``` r
library(devtools)
install_github("scollinselliott/arkhaia", dependencies = TRUE, build_vignettes = TRUE) 
```

## Methods

Vignettes in the package illustrate

### Homogeneity via Effect Size

Establishing a measure of practical signficance for the comparison of
finds assemblages in their depositional contexts. Whether or not a
context is “representative” of another is assessed on the basis of the
homogeneity of the distribution of finds. The relevant paper has been
reviewed and is under revision (Collins-Elliott Under Review).

- Cressie-Read power-divergence statistic to estimate $\chi^2$ (Cressie
  and Read 1984; Read and Cressie 1988)
- Bergsma’s bias-corrected Cramér’s $V$ (Bergsma 2013)
- Effect size as a measure of homogeneity between “related” and
  “unrelated” archaeological contexts, evaluated on the basis of count
  data and presence/absence data.
- Leave-one-out (LOO) validation of effect sizes.
- Inference based on sign error and practical signficance.

### Random Right-Censored Count Data

Random right-censoring of archaeological count data. Assemblages of
artifacts that do *not* belong to primary contexts (i.e., “random”
secondary or tertiary contexts), comprise a minimum amount of finds that
were “is use” in a given locality. From a contingency table of those
minimum counts, it is possible to generate distributions of counts in
use from which to evaluate This paper is currently under review
(Collins-Elliott Under Review).

- Estimating the rate of a Poisson distribution based on a contingency
  table of minimum counts.
- Resampling routines to generate contingency tables accoriding to a
  truncated Poisson distribution, representing counts “in use” as
  opposed to those deposied.

### Least Squares Sepctral Analysis

Least squares spectral analysis (LSSA), with an implementation of
fitting by lowest frequency iteratively (LSSA-LFI), to fit sparse
time-indexed observations and then evaluate whether there exists linear
dependence in their data-generating process via model comparison. The
paper applying this method to Babylonian price data is under review
(Collins-Elliott Under Review).

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-bergsma_bias-correction_2013" class="csl-entry">

Bergsma, W. 2013. “A Bias-Correction for Cramér’s \$V\$ and Tschuprow’s
\$T\$.” *Journal of the Korean Statistical Society* 42: 323–28.
<https://doi.org/10.1016/j.jkss.2012.10.002>.

</div>

<div id="ref-collins-elliott_evaluating_nodate" class="csl-entry">

Collins-Elliott, S. A. Under Review. “Evaluating the Relationship
Between Surface, Subsurface, and Stratigraphic Assemblages,” Under
Review.

</div>

<div id="ref-collins-elliott_random_nodate" class="csl-entry">

———. Under Review. “Random Right Censoring of Archaeological Count
Data,” Under Review.

</div>

<div id="ref-collins-elliott_revisiting_nodate" class="csl-entry">

———. Under Review. “Revisiting Babylonian Prices: Long-Run and
Variable-Length Equilibria, <span class="nocase">ca.</span> 400-80 BCE,”
Under Review.

</div>

<div id="ref-cressie_multinomial_1984" class="csl-entry">

Cressie, N. A. C., and T. R. C. Read. 1984. “Multinomial Goodness-of-Fit
Tests.” *Journal of the Royal Statistical Society. Series B
(Methodological)* 46: 440–64.
<https://doi.org/10.1111/j.2517-6161.1984.tb01318.x>.

</div>

<div id="ref-eddelbuettel_extending_2018" class="csl-entry">

Eddelbuettel, D., and J. J. Balamuta. 2018. “Extending R with C++: A
Brief Introduction to Rcpp.” *The American Statistician* 72: 28–36.
<https://doi.org/10.1080/00031305.2017.1375990>.

</div>

<div id="ref-eddelbuettel_rcpparmadillo_2014" class="csl-entry">

Eddelbuettel, D., and C. Sanderson. 2014. “RcppArmadillo: Accelerating R
<span class="nocase">with</span>
<span class="nocase">high</span>-<span class="nocase">performance</span>
C++ <span class="nocase">linear</span>
<span class="nocase">algebra</span>.” *Computational Statistics and Data
Analysis* 71: 1054–63. <https://doi.org/10.1016/j.csda.2013.02.005>.

</div>

<div id="ref-read_goodnessfit_1988" class="csl-entry">

Read, T. R. C., and N. A. C. Cressie. 1988. *Goodness-of-Fit Statistics
for Discrete Multivariate Data*. New York: Springer.

</div>

</div>
