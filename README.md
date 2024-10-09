# BayesianPairedSamples

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`BayesianPairedSamples` is an R package designed to perform Bayesian analysis for pretest-posttest binary outcomes, with adaptive significance levels. This analysis is useful in situations where pairs of measurements come from the same individual before and after treatment, like in pretest-posttest studies. The adaptive significance levels are computed objectively from the tests.

## Installation

You can install the latest version of the package directly from GitHub using `remotes` or `devtools`:

```r
# Install the remotes package if you don't have it already
install.packages("remotes")

# Install the package from GitHub
remotes::install_github("torodriguezt/BayesianPairedSamples")
