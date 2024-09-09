# wINEQ 


# *News*

# wINEQ 1.2.1 _(2024-08-03)_

## New features

* Package has been enriched with additional functions, that is:
  + `Quantile` computes quantile derived for the given probability taking into account weights
  + `LowerSum` computes weighted sum of values lower then a quantile derived for the given probability. It is a part of very well-know inequality indices. 
* Function `Theil_T`, `Theil_L` and `Entropy` have an option to deal with zeroes in a data vector.
  
## Bugs fixed

* Functions `Palma` and `Prop20_20` incorrectly determined quantiles. 

## Changes

* Function `Gini` have an alternative computing algorithm: slower but memory saving. 

## Future works

* new inequality measures
* vignette presenting an example workflow with our package
* enhanced dataset with tourists' expenditures on trips




# wINEQ 1.2.0 _(2023-04-22)_

## New features

* Package has been enriched with new inequality measures for variables on an ordinal scale, that is:
  + Apouey index
  + Abul Naga and Yalcin index
  + Blair and Lacy index
* Functions `ineq.weighted` and `ineq.weighted.boot` now start with assessing if variable of interest is numeric whether ordered factor. Then the inequality measures appropriate for the given variable are returned.
* New function `medianF` returns median for an ordered factor.


## Changes

* Improvement of handling with missing data or wrong data format.
* Allison and Foster index formula is not suitable for an ordered factor. Hence, function `AF` converts an ordered factor variable to numeric variable and then performs calculation.  




# wINEQ 1.1.2 _(2023-02-16)_

## Bugs fixed

* Function `ineq.weighted.boot` worked incorrectly for calibrated bootstrap storing results in an object of wrong dimensions. 

## Changes

* Function `ineq.weighted.boot` returns named list to access its items easily. 

## Future works

* more inequality measures for ordinal data


# wINEQ 1.1.1 _(2023-01-13)_

## Changes

* R Documentation is enriched with formula for each inequality measure, more details and more references.
* `Leti` function works for numerical and now also for ordered factor variable. 
* A warning was added for `Kolm` formula that the result is scale dependent. Example provided. `Kolm` function is enriched with some standardization methods.

## New features

* Added new data `Well_being` - data from sample survey on quality of life conducted on Polish-Ukrainian border. Modified due to SDC.
* Added new data `Tourism` - data from sample survey on trips conducted in Polish households. Modified due to SDC.
* Added examples for new datasets.

## Bugs fixed

* Function `ineq.weighted.boot` did not take a default `bounds` parameter for `calib` function. Now `bounds` parameter is added with default values.
* `Leti` function incorrectly calculated a cumulative distribution. Fixed.
* `Leti` function standardized its results to unit interval without option to change it. Additional parameter `norm` provided.



# wINEQ 1.0.1 _(2022-02-11)_

## First CRAN Release

* Initial CRAN Release



