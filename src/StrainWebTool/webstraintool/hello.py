import os.path
import webstraintool

from webstraintool import app
from flask import flash, render_template, request, redirect, url_for
from werkzeug.utils import secure_filename

## for python version
import platform
print platform.python_version()

@app.route('/')
def index():
    return 'This is the index!'

@app.route('/hello')
def hello_world():
    return render_template('hello.html')

@app.route('/upload', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        # check if the post request has the file part
        if 'file' not in request.files:
            flash('No file part', 'error')
            return redirect(request.url)
        file = request.files['file']
        # if user does not select file, browser also
        # submit a empty part without filename
        if file.filename == '':
            flash('No selected file', 'error')
            return redirect(request.url)
        if file:
            filename = secure_filename(file.filename) + '.tmp'
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            try:
                stations = parse_ascii_input(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            except:
                #print '[DEBUG] Error parsing input file!'
                flash('Failed to parse input file', 'error')
                return redirect(request.url)
        return render_template('hello.html', stations=stations)
    return render_template('hello.html')
