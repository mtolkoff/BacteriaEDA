import ast
from flask import request

from COloNY_functions.FLASK_functions import app
from COloNY_functions.helper_functions.getPlot import *

#Graph the TSNE
@app.route('/results_graph', methods = ['POST', 'GET'])
def results_graph():
    #Read in information
    ratios = request.args.get('ratios')
    ratios = ast.literal_eval(ratios)


    ratios6 = request.args.get('ratios6')
    ratios6 = ast.literal_eval(ratios6)


    ratios7 = request.args.get('ratios7')
    ratios7 = ast.literal_eval(ratios7)


    ratios_other = request.args.get('ratios_other')
    ratios_other = ast.literal_eval(ratios_other)



    x = request.args.get('x_axis')
    x = ast.literal_eval(x)



    y = request.args.get('y_axis')
    y = ast.literal_eval(y)



    OTU = request.args.get('OTU')
    OTU = ast.literal_eval(OTU)



    taxa_name = request.args.get('taxa_name')
    taxa_name = ast.literal_eval(taxa_name.replace("&#39;", "'"))
    taxa_name = ast.literal_eval(taxa_name)

    dbScanResults = request.args.get('dbScanResults')
    dbScanResults = ast.literal_eval(dbScanResults)

    colorBy = request.args.get('color_by')
    
    #Run get plot function
    return getPlot(x, y, colorBy, ratios, ratios6, ratios7, ratios_other, OTU, taxa_name, dbScanResults)

