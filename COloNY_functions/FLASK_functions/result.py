import pandas as pd
import numpy as np

from flask import request
from COloNY_functions.helper_functions.species_list import *
from COloNY_functions.helper_functions.getHTML import *
from COloNY_functions.FLASK_functions import app

data = pd.read_csv('picrust_Named.csv', index_col = 0)



#handles results and different requests
@app.route('/result', methods = ['POST', 'GET'])
def result():
    #preselected groups/orders
    if request.form["submit"] == 'Compute':
        if request.form["grouping"] == '300 Random Species':
            analyzer = data.sample(500, random_state = 1234)
        elif request.form["grouping"] == "Bacteroidales":
            analyzer = data.iloc[Bacteroidales, ]
        elif request.form["grouping"] == "Actinomycetales":
            analyzer = data.iloc[Actinomycetales, ]
        elif request.form["grouping"] == "Bacillales":
            analyzer = data.iloc[Bacillales, ]
        elif request.form["grouping"] == "Clostridiales":
            analyzer = data.iloc[Clostridiales, ]
        elif request.form["grouping"] == "500 Random Species":
            analyzer = data.sample(500, random_state = 12345)
    #handles random selection
    elif request.form["submit"] == "Random Selection":
        analyzer = data.sample(int(request.form["random"]))

    #handles uploaded file    
    elif request.form["submit"] == "Upload":
        analyzeList = []
        lines = request.files["uploader"].readlines()
        for line in lines:
            analyzeList.append(line.decode("utf-8").rstrip('\n'))
        analyzer = data[np.array(data['taxa_name'].isin(analyzeList), dtype=bool)]
    return getHTML(analyzer)

