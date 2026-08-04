# Simulate species occupancy in one dynamic landscape

Simulates the species' occupation on a landscape sequence, resorting to
the [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md)
function.

## Usage

``` r
simulate_graph(rl, rlist, simulate.start, method, parm, nsew="none", succ="none", 
param_df, kern, conn, colnz, ext, beta1, b, c1, c2, z, R)
```

## Arguments

- rl:

  Object of class 'landscape' or 'metapopulation'.

- rlist:

  List delivered by
  [`span.graph`](https://fmestre1.github.io/MetaLandSim/reference/span.graph.md).

- simulate.start:

  TRUE (rl is of class 'landscape') or FALSE (rl is of class
  'metapopulation')

- method:

  One of the following: **click** - individually select the patches with
  occurrence of the species by clicking on the map. Use only for
  individual landscape simulations. However, this option should not be
  used with iterate.graph. **percentage** - percentage of the patches to
  by occupied by the species. **number** - number of patches to be
  occupied by the species. To be internally passed to
  [`species.graph`](https://fmestre1.github.io/MetaLandSim/reference/species.graph.md).

- parm:

  Parameter to specify the species occurrence - either percentage of
  occupied patches or number of occupied patches, depending on the
  method chosen. To be internally passed to
  [`species.graph`](https://fmestre1.github.io/MetaLandSim/reference/species.graph.md).

- nsew:

  'N', 'S', 'E', 'W' or none - point of entry of the species in the
  landscape. By default set to "none". To be internally passed to
  [`species.graph`](https://fmestre1.github.io/MetaLandSim/reference/species.graph.md).

- succ:

  Set the preference of the species for patch successional stage:
  'none', 'early', 'mid' and 'late'.

- param_df:

  Parameter data frame delivered by
  [`parameter.estimate`](https://fmestre1.github.io/MetaLandSim/reference/parameter.estimate.md),
  including:

  - alpha - Parameter relating extinction with distance.

  - y - Parameter y in the colonization probability.

  - e - Parameter defining the extinction probability in a patch of unit
    area.

  - x - Parameter scaling extinction risk with patch area.

  To be internally passed to `simulate_graph`.

- kern:

  'op1' or 'op2'. Dispersal kernel. See details in the
  [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md)
  function. To be internally passed to
  [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md).

- conn:

  'op1' or 'op2'. Connectivity function. See details in the
  [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md)
  function. To be internally passed to
  [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md).

- colnz:

  'op1', 'op2' or 'op3'. Colonization function. See details in the
  [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md)
  function. To be internally passed to
  [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md).

- ext:

  'op1', 'op2' or 'op3'. Extinction function. See details in the
  [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md)
  function. To be internally passed to
  [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md).

- beta1:

  Parameter afecting long distance dispersal probability (if the
  Kern='op2'). To be internally passed to
  [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md).

- b:

  Parameter scaling emigration with patch area (if conn='op1' or 'op2').
  To be internally passed to
  [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md).

- c1:

  Parameter scaling immigration with the focal patch area (if
  conn='op2'). To be internally passed to
  [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md).

- c2:

  Parameter c in the option 3 of the colonization probability (if
  colnz='op3'). To be internally passed to
  [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md).

- z:

  Parameter giving the strength of the Allee effect (if colnz='op3'). To
  be internally passed to
  [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md).

- R:

  Parameter giving the strength of the Rescue effect (if ext='op3'). To
  be internally passed to
  [`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md).

## Value

Returns a list of occupied landscapes, representing the same occupied
landscape at different time steps.

## Author

Frederico Mestre and Fernando Canovas

## See also

[`spom`](https://fmestre1.github.io/MetaLandSim/reference/spom.md),
[`span.graph`](https://fmestre1.github.io/MetaLandSim/reference/span.graph.md),
[`rland.graph`](https://fmestre1.github.io/MetaLandSim/reference/rland.graph.md),
[`iterate.graph`](https://fmestre1.github.io/MetaLandSim/reference/iterate.graph.md)

## Examples

``` r

data(rland)
data(landscape_change)
data(param1)

sim1 <- simulate_graph(rl=rland, 
      rlist=landscape_change, 
      simulate.start=TRUE, 
      method="percentage", 
      parm=50, 
      nsew="none",
      succ = "none", 
      param_df=param1, 
      kern="op1", 
      conn="op1", 
      colnz="op1", 
      ext="op1", 
      beta1=NULL, 
      b=1, 
      c1=NULL, 
      c2=NULL, 
      z=NULL, 
      R=NULL
      )

```
