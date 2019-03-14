import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import mpld3
import json

from flask import render_template
from itertools import compress



#Monkeypatch broken dependencies
class  NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
            np.int16, np.int32, np.int64, np.uint8,
            np.uint16,np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32, 
            np.float64)):
            return float(obj)
        elif isinstance(obj,(np.ndarray,)): #### This is the fix
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

mpld3._display.NumpyEncoder = NumpyEncoder

def getPlot(x, y, colorBy, ratios, ratios6, ratios7, ratios_other, OTU, taxa_name, dbScanResults):

    Names = []
    Marker = []
    #names for the different taxa
    for i in range(len(taxa_name)):
        Names.append(
                'OTU#: ' + str(OTU[i]) +
                ' Taxa: ' + str(taxa_name[i]) +
                ' Cluter: ' + str(dbScanResults[i]) +
                ', Metabolic Genes: ' + str(ratios[i] * 100)[0:4] + '%')
        if taxa_name[i][0] == 'o' or taxa_name[i][0] == 'c' or taxa_name[i][0] == 'p' or taxa_name[i][0] == 'k':
            Marker.append("o")
        else:
            Marker.append("s")

    #Organizes ratios so that red to blue scale is always from most red to most blue
    maxRatio = max(ratios)
    minRatio = min(ratios)
    #color by metabolism
    colors = [0] * len(ratios)
    if(colorBy == "Metabolic Genes"):
        for i in range(len(colors)):
            adjustedRatio = (ratios[i] - minRatio) / (maxRatio - minRatio)
            colors[i] = matplotlib.colors.to_hex([1 - adjustedRatio, 0, adjustedRatio])
 
 #non-metabolic coloring
    else:
        groups = []
        colorMap = []
        colors = []
        count = 0
        for item in dbScanResults:
            if item in groups:
                colors.append(colorMap[groups.index(item)])
            else:
                groups.append(item)
                tmpColor = np.random.rand(3,)
                colorMap.append(tmpColor)
                colors.append(tmpColor)
                count += 1







    #plot results. Different markers require different plots
    resultsPlot, ax= plt.subplots(figsize=(12, 6.75) )
    plt.plot()
    
    #marks if family is known
    FamilyList = []
    for i in Marker:
        if i == "s":
            FamilyList.append(True)
        else:
            FamilyList.append(False)
    
    
    
    #plot results for family or lower being lowest category known
    xFamily = list(compress(x, FamilyList))
    yFamily = list(compress(y, FamilyList))
    colorsFamily = list(compress(colors, FamilyList))
    labelsFamily = list(compress(Names, FamilyList))
    labels = Names
    points = plt.scatter(xFamily, yFamily, alpha = .5, s = 100, c = colorsFamily, marker = "o")
    mpld3.plugins.connect(resultsPlot, mpld3.plugins.PointLabelTooltip(points, labels = labelsFamily))

    #Plots those for order (or higher) is lowest category known
    orderList = []
    for i in Marker:
        if i == "o":
            orderList.append(True)
        else:
            orderList.append(False)
    xOrder = list(compress(x, orderList))
    yOrder = list(compress(y, orderList))
    colorsOrder = list(compress(colors, orderList))
    labelsOrder = list(compress(Names, orderList))
    labels = Names
    points = plt.scatter(xOrder, yOrder, alpha = .5, s = 100, c = colorsOrder, marker = "v")
    mpld3.plugins.connect(resultsPlot, mpld3.plugins.PointLabelTooltip(points, labels = labelsOrder))



    ax.axis('off')

    return render_template('results_graph.html', figure=mpld3.fig_to_html(resultsPlot))

