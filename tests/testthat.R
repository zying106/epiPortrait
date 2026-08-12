library(testthat)
library(epiPortrait)

# Ensure data is available
if (!exists("example_se")) {
  data(example_se)
}

test_check("epiPortrait")
