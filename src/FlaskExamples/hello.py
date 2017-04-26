from flask import Flask
app = Flask(__name__)

@app.route('/')
def index():
    return 'This is the index!'

@app.route('/hello')
def hello_world():
    return 'Hello, World!'
