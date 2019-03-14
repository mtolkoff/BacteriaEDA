import pandas as pd

otuToNames_file = open('../Downloads/99_otu_taxonomy.txt', 'r')
data_test = pd.read_csv('picrust_midSize', index_col = 0)
otuToName = otuToNames_file.read()

otuToNameSplit = []
for i in range(len(otuToName) - 1):
     otuToNameSplit.append(otuToName[i].split('\t'))

#otuToNameSplit = np.transpose(otuToNameSplit)

dataNames = []
for i in range(5000):
    dataNames.append(np.transpose(otuToNameSplit)[1][list(np.transpose(otuToNameSplit)[0]).index(str(data_test.index[i]))])

finalName = []
for i in range(5000):
    count = 6
    isFound = False
    treeLine = dataNames[i].split("; ")
    while not isFound:
        if len(treeLine[count]) <= 4:
                count-= 1
        else:
                finalName.append(treeLine[count].rstrip('\n')
                isFound = True

data_test['taxa_name'] = finalName
