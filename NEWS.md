# wINEQ 


# *News*

# wINEQ 1.1.0 _(2023-01-13)_

## Changes

* R Documentation is enriched with formula for each inequality measure, more details and more references.
* `Leti` function works for numerical and now also for ordered factor variable. 
* A warning was added for `Kolm` formula that the result is scale dependent. Example provided. `Kolm` function is enriched with some standardization methods.

## New features

* Added new data `Well_being` - data from sample survey on quality of life conducted on Polish-Ukrainian border. Modified due to SDC.
* Added new data `Tourism` - data from sample survey on trips conducted in Polish households. Modified due to SDC.
* Added examples for new datasets.

## Bugs fixed:

* Function `ineq.weighted.boot` did not take a default `bounds` parameter for `calib` function. Now `bounds` parameter is added with default values.
* `Leti` function incorrectly calculated a cumulative distribution. Fixed.
* `Leti` function standardized its results to unit interval without option to change it. Additional parameter `norm` provided.



# wINEQ 1.0.1 _(2022-02-11)_

## First CRAN Release

* Initial CRAN Release



