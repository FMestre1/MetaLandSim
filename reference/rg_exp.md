# List with range.expansion output

Output of
[`range_expansion`](https://fmestre1.github.io/MetaLandSim/reference/range_expansion.md).
Object of class 'expansion'.

## Usage

``` r
data(rg_exp)
```

## Format

Data frame with the probability of occupations at several distances from
the closest occupied landscape mosaic. The data frame has the following
columns:

- DISTANCE - Distance (mapsize x number of landscapes).

- OCCUPATION - How many times did the landscape at this distance got
  occupied by the species (from a total of 'iter' repetitions).

- PROPORTION - Proportion of occupation for the landscape at this
  distance (OCCUPATION/iter).

- TIME STEP - The average time step during which a given distance is
  reached.

## Examples

``` r

data(rg_exp)
```
