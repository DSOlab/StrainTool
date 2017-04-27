from app import app
from flask import render_template

@app.route('/')
def index():
    return 'This is the index!'

@app.route('/hello')
def hello_world():
    return render_template('hello.html')
