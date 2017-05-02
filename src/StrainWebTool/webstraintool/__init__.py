# Import flask and template operators
import os.path
from flask import Flask

# Define the WSGI application object
app = Flask(__name__)
from webstraintool import hello
# Configurations
app.config.from_object('config')
# Temporary uploading area (this should somehow be in the config file)
app.config["UPLOAD_FOLDER"] = os.path.join(os.path.abspath(os.path.dirname('.')), 'tmp')
