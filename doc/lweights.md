## Test Programs

* `lweights_test.py`
Reference: _Optimal Interpolation of Spatially Discretized Geodetic Data, Shen et al, 2015_

[...]We reconstruct the covariance matrix `C` by multiplying a
weighting function to each of its diagonal terms `C<sub>i</sub>` , and the
weighting is given as `C<sub>i</sub> <- C<sub>i</sub> G^1^<sub>i</sub>`. The weighting function
`G<sub>i</sub>=L<sub>i</sub>×Z<sub>i</sub>`, in which `L<sub>i</sub>` and `Z<sub>i</sub>` are functions of distance
and spatial coverage dependent, respectively. For distance-
dependent weighting, `L<sub>i</sub>` is assumed to be in the form of
```L<sub>i</sub> = exp(−ΔR<sup>2</sup><sub>i</sub>/D^2^)```
or
```L<sub>i</sub> = 1/(1+ΔR<sup>2</sup><sub>i</sub>/D^2^)```
in which a spatial smoothing parameter `D` is introduced. Both
functions allow reduced weight of the data as distance in creases.
