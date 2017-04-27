import os.path
from app import app
from flask import flash, render_template, request, redirect, url_for
from werkzeug.utils import secure_filename

g_stations = [];

@app.route('/')
def index():
    return 'This is the index!'

@app.route('/hello')
def hello_world():
    return render_template('hello.html')

class Station:
    def __init__(self, input_line):
        ''' parse from standard ascii input line '''
        l = input_line.split()
        self.name = l[0]
        self.lon  = float(l[1])
        self.lat  = float(l[2])
        self.ve   = float(l[3])
        self.vn   = float(l[4])
        self.se   = float(l[5])
        self.sn   = float(l[6])
        self.rho  = float(l[7])
        self.t    = float(l[8])

def parse_ascii_input(filename):
    stations = []
    with open(filename) as fin:
        for line in fin.readlines():
            stations.append(Station(line))
    print '[DEBUG] Read {} stations.'.format(len(stations))
    if len(stations):
        return stations
    else:
        return None

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
                print '[DEBUG] Error parsing input file!'
                flash('Failed to parse input file', 'error')
                return redirect(request.url)
        return render_template('hello.html', stations=stations)
    return render_template('hello.html')
