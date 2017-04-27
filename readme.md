Basic Setup
==================

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

To start the application, run the `run.py` script (in the `root_folder`). You should see something like:
>  * Running on http://127.0.0.1:5000/ (Press CTRL+C to quit)
>  * Restarting with stat
>  * Debugger is active!
>  * Debugger PIN: 201-460-604

Go to `http://127.0.0.1:5000/` to see it!
