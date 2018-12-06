#read file
otuToNames_file = open('99_otu_taxonomy.txt', 'r')
otuToName = otuToNames_file.read()
otuToName = otuToName.split('\n')

#splits taxonomic groupings from OTU numbers
otuToNameSplit = []
for i in range(len(otuToName) - 1):
     otuToNameSplit.append(otuToName[i].split('\t'))

#otuToNameSplit = np.transpose(otuToNameSplit)

#Replaces OTU number with full taxonomic information
dataNames = []
for i in range(5000):
    data_test.iloc[i, 0] = np.transpose(otuToNameSplit)[1][list(np.transpose(otuToNameSplit)[0]).index(str(data_test.iloc[i, 0]))]

#replaces full taxonomic information with just name of lowest level group
for i in range(5000):
    count = 6
    isFound = False
    treeLine = data_test.iloc[i, 0].split("; ")
    while not isFound:
        if len(treeLine[count]) <= 4:
                count-= 1
        else:
                data_test.iloc[i, 0] = treeLine[count].rstrip('\n')
                isFound = True

