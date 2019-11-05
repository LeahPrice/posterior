<!-- badges: start -->
[![Build Status](https://travis-ci.org/jgabry/posterior.svg?branch=master)](https://travis-ci.org/jgabry/posterior)
[![Coverage Status](https://codecov.io/github/jgabry/posterior/coverage.svg?branch=master)](https://codecov.io/github/jgabry/posterior?branch=master)
<!-- badges: end -->

# posterior

The **posterior** R package is intended to provide useful tools for both users
and developers of packages for fitting Bayesian models or working with output
from Bayesian models. The primary goals of the package are to:

* Efficiently convert between many different useful formats of
  draws (samples) from posterior or prior distributions.
* Provide various summaries of draws in convenient formats.
* Provide lightweight implementations of state of the art MCMC diagnostics.


### Installation

The package is not released yet. Currently you can install the development
version from GitHub, but expect frequent changes until an official release.

```r
# install.packages("devtools")
devtools::install_github("jgabry/posterior")
```
