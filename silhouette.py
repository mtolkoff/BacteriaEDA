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

with open('koOrthologs.json') as f:
    categories = json.load(f)
     
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


def getSilhouette(analyzer):

    analyzerCount = analyzer.sum()        
    analyzer = analyzer.loc[:, analyzerCount != 0]


#        LDA = decomp.LatentDirichletAllocation(n_components = 3)
    tsne = manifold.TSNE()
    Names = list(np.transpose(analyzer.iloc[:, 0:1].values)[0])
    namelessData = analyzer.iloc[:, 1:]
    Features = list(namelessData)
    tsne.fit(namelessData)
    

    colors = [] 
    groupList = []
    FullNames = []
    
   
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
    
    
    x = np.transpose(tsne.embedding_)[0]
    y = np.transpose(tsne.embedding_)[1] 


    dbscan = cluster.DBSCAN(eps = 4.5, min_samples = 10)
    dbScanResults = dbscan.fit(tsne.embedding_)
    #dbScanResults.labels_
    for i in range(len(Names)):
        Names[i] = str(Names[i]) + " DBSCAN label: " + str(dbScanResults.labels_[i]) + ", Group 1:  " + str(ratios[i] * 100)[0:4] + "%, Group 2: " + str(ratios6[i] * 100)[0:4] + "%, Group 3: " + str(ratios7[i] * 100)[0:4] + "%, Other: " + str(ratios_other[i] * 100)[0:4] + "%"


    silhouette = []
    silhouette.append(metrics.silhouette_score(tsne.embedding_, dbScanResults.labels_, metric = 'euclidean'))
    silhouette.append(metrics.silhouette_score(namelessData, dbScanResults.labels_, metric = 'euclidean'))
    silhouette.append(metrics.silhouette_score(np.array(ratios).reshape(-1, 1), dbScanResults.labels_, metric = 'euclidean'))
    silhouette.append(metrics.silhouette_score(np.array(ratios6).reshape(-1, 1), dbScanResults.labels_, metric = 'euclidean'))
    silhouette.append(metrics.silhouette_score(np.array(ratios7).reshape(-1, 1), dbScanResults.labels_, metric = 'euclidean'))
    silhouette.append(metrics.silhouette_score(np.array(ratios_other).reshape(-1, 1), dbScanResults.labels_, metric = 'euclidean'))

    return silhouette 


count = 0
silFinal = [0, 0, 0, 0, 0, 0]

for i in range(100):
    print(count)
    count += 1
    sil = getSilhouette(data.sample(500))
    for i in range(6):
        silFinal[i] += sil[i]

for i in range(6):
    print(silFinal[i] / 100)



