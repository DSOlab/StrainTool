## Test Programs

* `lweights_test.py`
Reference: _Optimal Interpolation of Spatially Discretized Geodetic Data, Shen et al, 2015_

[...]We reconstruct the covariance matrix `C` by multiplying a
weighting function to each of its diagonal terms C~i~ , and the
weighting is given as C~i~ <- C~i~ G^1^~i~ . The weighting function
G~i~=L~i~×Z~i~, in which L~i~ and Z~i~ are functions of distance
and spatial coverage dependent, respectively. For distance-
dependent weighting, L~i~ is assumed to be in the form of
L~i~ = exp(−ΔR^2^~i~/D^2^) 
L~i~ = 1/(1+ΔR^2^~i~/D^2^)
in which a spatial smoothing parameter `D` is introduced. Both
functions allow reduced weight of the data as distance in-
creases.
