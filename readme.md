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
