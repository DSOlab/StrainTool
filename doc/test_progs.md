## Test Programs

### `zweights_test.py`

This program is used to check the computation of **Z** weights, according to <sup>[1](#shen2015)</sup>.
Run as `zweights_test.py -i ../tmp/input2.vel`, it will write all computation
info to stdout and print a map with the stations.

The program is based on the function `pystrain::z_weights`; it will:
 1. Parse all stations from the input list (`-i` switch)
 2. Transform all station coordinates to UTM (after computing the UTM zone)
 3. Compute the barycenter
 4. Compute the z-weights
    * Compute the azimouth of the line connecting the barycenter to each of the points
    * Sort the azimouths (in ascending order)
    * Compute the θ angle of each point, as an azimouth difference
    * Compute the individual z-weights, using z=θ*n/4π
    * Restore the z weight list in the order the stations where passed in
 5. Print results and plot map

All intermediate results are printed to check the algorithm.

### `lweights_test.py`

This program is used to check the computation of **L** weights, according to <sup>[1](#shen2015)</sup>.
The L (spatial) weights, are computed using one of the two functions:
* L(i) = exp(-ΔR(i)<sup>2</sup>/D<sup>2</sup>) -- __Gaussian__, or
* L(i) = 1/(1+ΔR(i)<sup>2</sup>/D<sup>2</sup>) -- __Quadratic__
Note that to compute the weights, we need the parameter **D**.
Run as `lweights_test.py -i ../tmp/input2.vel [-t gaussian|quadratic] [-D value_of_D_param]`

The program is based on the function `pystrain::l_weights`; it will:
 1. Parse all stations from the input list (`-i` switch)
 2. Transform all station coordinates to UTM (after computing the UTM zone)
 3. Compute the barycenter
 4. Compute z-weights
 5. Compute the l-weights
    * Compute distance from each point to the barycenter
    * if parameter D is given, then just compute the l-weights and return them
    * if parameter D is not given, then iterate through the range `[dmin, dmax]`
      with a step of `dstep`; for each value, compute the l-weights and the
      quantity `W=Σ(l*z)` (aka sum of z- and l-weights). If `W >= Wt`, then this
      value is the optimal D, and the l-weights are returned.
 6. Print results and plot map

All intermediate results are printed to check the algorithm.

## References

<a name="shen2015">1</a>: _Optimal Interpolation of Spatially Discretized Geodetic Data, Shen et al, 2015_
