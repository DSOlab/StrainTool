# Basic Setup

## StrainTool WebApp (todo)
In the following, `root_dir` is the current dir, i.e. where this readme file 
is located. Go to `root_dir`.

```
$> virtualenv flask
$> cd flask
## activate the virtualenv
$> source bin/activate
## now you should see (flask) on the left of the command line. Let's install flask
$> pip install flask
$> deactivate
```

To start the application, run the `run.py` script (in the `root_folder`). You should see something
like:
>  * Running on http://127.0.0.1:5000/ (Press CTRL+C to quit)
>  * Restarting with stat
>  * Debugger is active!
>  * Debugger PIN: 201-460-604

Go to `http://127.0.0.1:5000/` to see it!

## StrainTool Stand-Alone Program

### Input Files

To estimate the strain parameters, the user needs to feed the program an input file (via the `-i` 
or `--input-file` switch). This input file, should contain station information in the following
format:
> station-name longtitude latitude  Ve     Vn     SigmaVe SigmaVn   Sne   time-span
>                 deg.       deg.  mm/yr  mm/yr    mm/yr   mm/yr     /    dec. years
Station coordinates are provided in longtitude/latitude pairs in *decimal degrees*. Velocities and
velocity standard deviations are provided in *mm per years (mm/yr)*. `Sne` is the correlation
coefficient between East and North velocity components.

**Note that at his point the last two columns (aka Sne and time-span) are not used, so they could have random values**.

### Weighting

Let `C` be the covariance matrix of the velocity data. Reminder:
*Cofactor Matrix Q* Q = C / σ<sub>0</sub><sup><2</sup>
*Weight Matrix* P = σ<sub>0</sub><sup><2</sup> * C<sup>-1</sup>
σ<sub>0</sub><sup><2</sup> is called variance factor or variance of unit weight or apriori variance factor
We reconstruct the covariance matrix C by multiplying a weighting function to 
each of its diagonal terms C<sub>i</sub> and the weighting is given as
C<sub>i</sub> <- C<sub>i</sub> * G<sup>-1</sup>. The weighting function
G<sub>i</sub> = L<sub>i</sub> * Z<sub>i</sub>, in which L<sub>i</sub> and 
Z<sub>i</sub> are functions of distance and spatial ccoverage dependent, respectively.
```
Wx = (1e-3/sta.se)*zw[idx]*lw[idx]
Wy = (1e-3/sta.sn)*zw[idx]*lw[idx]
```

