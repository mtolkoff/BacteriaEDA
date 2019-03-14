from flask import render_template
from COloNY_functions.FLASK_functions import app

@app.route('/')
def main():
    return render_template('index.html')

