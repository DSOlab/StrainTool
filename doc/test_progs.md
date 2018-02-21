## Test Programs

* `zweights_test.py`

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

* `lweights_test.py`

## References

<a name="shen2015">1</a>: _Optimal Interpolation of Spatially Discretized Geodetic Data, Shen et al, 2015_
