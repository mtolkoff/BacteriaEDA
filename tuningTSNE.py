import pandas as pd
import numpy as np
import numpy
import sklearn.metrics as metrics
import sklearn.manifold as manifold
import sklearn.decomposition as decomp
import sklearn.cluster as cluster
import sys
import json
from flask import Flask, render_template, request
import mpld3
import matplotlib.pyplot as plt
import matplotlib
import colorsys

#Monkeypatch broken dependencies
class  NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, (numpy.int_, numpy.intc, numpy.intp, numpy.int8,
            numpy.int16, numpy.int32, numpy.int64, numpy.uint8,
            numpy.uint16,numpy.uint32, numpy.uint64)):
            return int(obj)
        elif isinstance(obj, (numpy.float_, numpy.float16, numpy.float32, 
            numpy.float64)):
            return float(obj)
        elif isinstance(obj,(numpy.ndarray,)): #### This is the fix
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

mpld3._display.NumpyEncoder = NumpyEncoder


#Open KO
with open('koOrthologs.json') as f:
    categories = json.load(f)
     
#Post orer traversal to find path to KO ortholog to gene.
#Returns an array, [0] is the path used to find the ortholog value:
#example path: categories['children][0]['children'][1]['children'][1]['children'][13]
#[1] returns the name of the ortholog
def findPath(KO_num, path, node):
    if(path != ''):
        token = ', '
    else:
        token = ''
    if KO_num in node['name']: 
        path = list(map(int, path.split(",")))
        return([node['name'], path])
    elif 'children' in node.keys():
        for i in range(len(node['children'])):
            nodeBelow = findPath(KO_num, path + token + str(i), node['children'][i])
            if not nodeBelow == False:
                return nodeBelow
    return False

data = pd.read_table('picrustDB_Named', sep = ',')


#Builds the graph
def getHTML(analyzer, colorBy):
    
    #removes features which are all 0
    analyzerCount = analyzer.sum()        
    analyzer = analyzer.loc[:, analyzerCount != 0]


#        LDA = decomp.LatentDirichletAllocation(n_components = 3)
    tsne = manifold.TSNE()
    #get the names
    Names = list(np.transpose(analyzer.iloc[:, 0:1].values)[0])
    #removes the names so it's not part of the analysis
    namelessData = analyzer.iloc[:, 1:]
    Features = list(namelessData)
    #fit TSNE
    tsne.fit(namelessData)
    

    colors = [] 
    groupList = []
    FullNames = []
    
    #get names of orthologs 
    for name in list(namelessData):
        foundPath = findPath(name, '', categories)
        if foundPath == False:
            groupList.append(-1)
#            colors.append('white')
        else:
            groupList.append(foundPath[1][0])
#            if foundPath[1][0:1] in groupList:
#                colors.append(colorList[groupList.index(foundPath[1][0:1])])
#            else:
#                groupList.append(foundPath[1][0:1])
#                colorList.append(list(matplotlib.colors.CSS4_COLORS.items())[count][1])
#                colors.append(list(matplotlib.colors.CSS4_COLORS.items())[count][1])
#                count = (count + 1) % len(list(matplotlib.colors.CSS4_COLORS.items()))
    
    groupFreq = []
    groupFreqKey = []
    for i in groupList:
        if i in groupFreqKey:
            groupFreq[groupFreqKey.index(i)] += 1
        else:
            groupFreq.append(1)
            groupFreqKey.append(i)
    
    #compute ratios for each group
    totalGenes = [0] * namelessData.shape[0]
    groupGenes = [0] * namelessData.shape[0]
    ratios = [0] * namelessData.shape[0]
    
    
    for i in range(namelessData.shape[0]):
        for j in range(namelessData.shape[1]):
            if namelessData.iloc[i, j] == 1:
                totalGenes[i] += 1
                if groupList[j] == 0:
                    groupGenes[i] += 1

    for i in range(len(groupGenes)):
        ratios[i] = groupGenes[i] / totalGenes[i]
    #print(totalGenes)
    #print(groupGenes)
    #print(ratios)
    
    #group 2
    totalGenes7 = [0] * namelessData.shape[0]
    groupGenes7 = [0] * namelessData.shape[0]
    ratios7 = [0] * namelessData.shape[0]

    for i in range(namelessData.shape[0]):
        for j in range(namelessData.shape[1]):
            if namelessData.iloc[i, j] == 1:
                totalGenes7[i] += 1
                if groupList[j] == 7:
                    groupGenes7[i] += 1

    for i in range(len(groupGenes)):
        ratios7[i] = groupGenes7[i] / totalGenes7[i]
    
    #group 3
    totalGenes6 = [0] * namelessData.shape[0]
    groupGenes6 = [0] * namelessData.shape[0]
    ratios6 = [0] * namelessData.shape[0]

    for i in range(namelessData.shape[0]):
        for j in range(namelessData.shape[1]):
            if namelessData.iloc[i, j] == 1:
                totalGenes6[i] += 1
                if groupList[j] == 6:
                    groupGenes6[i] += 1

    for i in range(len(groupGenes)):
        ratios6[i] = groupGenes6[i] / totalGenes6[i]
    
    #other
    totalGenes_other = [0] * namelessData.shape[0]
    groupGenes_other = [0] * namelessData.shape[0]
    ratios_other = [0] * namelessData.shape[0]

    for i in range(namelessData.shape[0]):
        for j in range(namelessData.shape[1]):
            if namelessData.iloc[i, j] == 1:
                totalGenes_other[i] += 1
                if groupList[j] not in [0, 6, 7]:
                    groupGenes_other[i] += 1

    for i in range(len(groupGenes)):
        ratios_other[i] = groupGenes_other[i] / totalGenes_other[i]
    
    
    
    #print(ratios6)
    #print(ratios7)
    #print(ratios_other)
    
    #DBSCAN
    x = np.transpose(tsne.embedding_)[0]
    y = np.transpose(tsne.embedding_)[1] 


    dbscan = cluster.DBSCAN(eps = 4.5, min_samples = 10)
    dbScanResults = dbscan.fit(tsne.embedding_)

    #Labels
    for i in range(len(Names)):
        Names[i] = str(Names[i]) + " DBSCAN label: " + str(dbScanResults.labels_[i]) + ", Group 1:  " + str(ratios[i] * 100)[0:4] + "%, Group 2: " + str(ratios6[i] * 100)[0:4] + "%, Group 3: " + str(ratios7[i] * 100)[0:4] + "%, Other: " + str(ratios_other[i] * 100)[0:4] + "%"


    #Silhouette metrics if you're interested
    #print(metrics.silhouette_score(tsne.embedding_, dbScanResults.labels_, metric = 'euclidean'))
    #print(metrics.silhouette_score(namelessData, dbScanResults.labels_, metric = 'euclidean'))
    #print(metrics.silhouette_score(np.array(ratios).reshape(-1, 1), dbScanResults.labels_, metric = 'euclidean'))
    #print(metrics.silhouette_score(np.array(ratios6).reshape(-1, 1), dbScanResults.labels_, metric = 'euclidean'))
    #print(metrics.silhouette_score(np.array(ratios7).reshape(-1, 1), dbScanResults.labels_, metric = 'euclidean'))
    #print(metrics.silhouette_score(np.array(ratios_other).reshape(-1, 1), dbScanResults.labels_, metric = 'euclidean'))


    #metabolic coloring
    maxRatio = max(ratios)
    minRatio = min(ratios)
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
        for item in dbScanResults.labels_:
            if item in groups:
                colors.append(colorMap[groups.index(item)])
            else:
                groups.append(item)
                tmpColor = np.random.rand(3,)
                colorMap.append(tmpColor)
                colors.append(tmpColor)
                count += 1

   
    g1 = categories['children'][0]['name']
    g2 = categories['children'][7]['name']
    g3 = categories['children'][6]['name']
  
    #plot results
    resultsPlot, ax= plt.subplots(figsize=(12, 6.75) )
    ax.plot()
    points = ax.scatter(x, y, alpha = .5, s = 100, c = colors)
    labels = Names
    mpld3.plugins.connect(resultsPlot, mpld3.plugins.PointLabelTooltip(points, labels = labels))


    ax.axis('off')
    return render_template('resultMPL.html', figure=mpld3.fig_to_html(resultsPlot)) 



#analyzer1 = data.sample(500, random_state = 1234)
#analyzer2 = data.sample(500, random_state = 12345)

#htmlPreselect1_metabolism = getHTML(analyzer1, "Metabolic Genes")
#htmlPreselect2_metabolism = getHTML(analyzer2, "Metabolic Genes")
#htmlPreselect1_Cluster =  getHTML(analyzer1, "Metabolrc Genes")
#htmlPreselect2_Cluster =  getHTML(analyzer2, "Metabolrc Genes")




app = Flask(__name__)

@app.route('/')
def main():
    speciesList = data.iloc[:, 0:1].values
    speciesListFormatted = [item[0] for item in speciesList]
    return render_template('index.html', speciesList = speciesListFormatted)



#handles results and different requests
@app.route('/result', methods = ['POST', 'GET'])
def result():
    colorBy = request.form["colorBy"]
    if request.form["submit"] == 'Preselected 1':
        #if(colorBy == "Metabolic Genes"):
        #    return htmlPreselect1_metabolism
        #else:
        #    return htmlPreselect1_Cluster
        analyzer = data.sample(500, random_state = 1234)
    elif request.form["submit"] == "Preselected 2":
        #if(colorBy == "Metabolic Genes"):
        #    return htmlPreselect2_metabolism
        #else:
        #    return htmlPreselect2_Cluster
        analyzer = data.sample(500, random_state = 12345)
    elif request.form["submit"] == "Random Selection":
        analyzer = data.sample(int(request.form["random"]))
    elif request.form["submit"] == "Upload":
        analyzeList = []
        lines = request.files["uploader"].readlines()
        for line in lines:
            analyzeList.append(line.decode("utf-8").rstrip('\n'))
        #result = request.form 
        #analyzeList = [item[0] for item in result.items()]
        #print(analyzeList)
        #analyzeListInt = list(map(int, analyzeList))
        #print(analyzeList)
        #print(data.iloc[:, 0:1].isin(analyzeList).sum())
        analyzer = data[np.array(data.iloc[:, 0].isin(analyzeList), dtype=bool)]
    return getHTML(analyzer, colorBy)


if __name__ == '__main__':
    app.run(debug = True)
    
