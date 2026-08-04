# Set of simulated data to use with the IFM parameter estimation functions. The data were generated using the code provided in "details".

This dataset loads several objects:

'sim.area', 'sim.det.20', 'sim.distance', 'z.sim', 'z.sim.20' and
'z.sim.20.fa'.

## Format

- `sim.area`:

  Vector of the areas for each site; here, 100 sites.

- `sim.det.20`:

  100 x 10 x 3 array corresponding to nsites x nyears x nvisits. Data
  simulated with year-specific detection probabilities equal to
  0.4,0.6,0.2,0.9,0.3,0.4,0.6,0.2,0.9,0.3.

- `sim.distance`:

  nsite x nsite distance matrix.

- `z.sim`:

  nyear x nsite occupancy data generated with perfect detection.

- `z.sim.20`:

  nyear x nsite occupancy data generated with perfect detection with
  approximately 20% of data missing at random.

- `z.sim.20.fa`:

  nyear x nsite occupancy data containing false absences, which can be
  used to explore the bias of ifm.missing.MCMC and ifm.naive.MCMC when
  there is imperfect detection.

## Details

These datasets were created using the code included in **Examples**.

## Examples
