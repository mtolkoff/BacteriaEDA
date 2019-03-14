import numpy as np
import json
from flask import render_template

import sklearn.metrics as metrics
import sklearn.manifold as manifold
import sklearn.decomposition as decomp
import sklearn.cluster as cluster

from COloNY_functions.helper_functions.findPath import *

#Open KO orthologs
with open('koOrthologs.json') as f:
    categories = json.load(f)

#Assemble information for the results page
def getHTML(analyzer):
    
    #removes features which are all 0
    analyzerCount = analyzer.sum()         
    taxa_name = list(analyzer['taxa_name'])
    analyzer = analyzer.drop('taxa_name', axis = 1)
    analyzer = analyzer.loc[:, analyzerCount != 0]


    tsne = manifold.TSNE()
    #get the names
    OTU = analyzer.index
    analyzer = analyzer
    Features = list(analyzer)
    #fit TSNE
    tsne.fit(analyzer)
    

    colors = [] 
    groupList = []
    FullNames = []
    
    #get names of orthologs 
    for name in list(analyzer):
        foundPath = findPath(name, '', categories)
        if foundPath == False:
            groupList.append(-1)
        else:
            groupList.append(foundPath[1][0])
    
    #compute percentage of genes in each group for each species
    totalGenes = [0] * analyzer.shape[0]
    groupGenes = [0] * analyzer.shape[0]
    ratios = [0] * analyzer.shape[0] 
    
    for i in range(analyzer.shape[0]):
        for j in range(analyzer.shape[1]):
            if analyzer.iloc[i, j] == 1:
                totalGenes[i] += 1
                if groupList[j] == 0:
                    groupGenes[i] += 1

    for i in range(len(groupGenes)):
        ratios[i] = groupGenes[i] / totalGenes[i]


    #group 2
    totalGenes7 = [0] * analyzer.shape[0]
    groupGenes7 = [0] * analyzer.shape[0]
    ratios7 = [0] * analyzer.shape[0]

    for i in range(analyzer.shape[0]):
        for j in range(analyzer.shape[1]):
            if analyzer.iloc[i, j] == 1:
                totalGenes7[i] += 1
                if groupList[j] == 7:
                    groupGenes7[i] += 1

    for i in range(len(groupGenes)):
        ratios7[i] = groupGenes7[i] / totalGenes7[i]
    
    #group 3
    totalGenes6 = [0] * analyzer.shape[0]
    groupGenes6 = [0] * analyzer.shape[0]
    ratios6 = [0] * analyzer.shape[0]

    for i in range(analyzer.shape[0]):
        for j in range(analyzer.shape[1]):
            if analyzer.iloc[i, j] == 1:
                totalGenes6[i] += 1
                if groupList[j] == 6:
                    groupGenes6[i] += 1

    for i in range(len(groupGenes)):
        ratios6[i] = groupGenes6[i] / totalGenes6[i]
    
    #other
    totalGenes_other = [0] * analyzer.shape[0]
    groupGenes_other = [0] * analyzer.shape[0]
    ratios_other = [0] * analyzer.shape[0]

    for i in range(analyzer.shape[0]):
        for j in range(analyzer.shape[1]):
            if analyzer.iloc[i, j] == 1:
                totalGenes_other[i] += 1
                if groupList[j] not in [0, 6, 7]:
                    groupGenes_other[i] += 1

    for i in range(len(groupGenes)):
        ratios_other[i] = groupGenes_other[i] / totalGenes_other[i]
    
    
    
    #DBSCAN
    x = np.transpose(tsne.embedding_)[0]
    y = np.transpose(tsne.embedding_)[1] 


    dbscan = cluster.DBSCAN(eps = 4.5, min_samples = 10)
    dbScanResults = dbscan.fit(tsne.embedding_)


    
    TaxaLabel = np.sort(OTU)

    return render_template('resultMPL.html', taxa = TaxaLabel, x_axis = list(x), y_axis = list(y), ratios = ratios, ratios6 = ratios6, ratios7 = ratios7, ratios_other = ratios_other, OTU = list(OTU), taxa_name = list(taxa_name), dbScanResults = list(dbScanResults.labels_))
