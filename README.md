# COloNY Description
This was a project I created with Insight Data Science in consultation with Indigo Agriculture. Indigo is a company which makes "seed coatings," or bacteria and fungi which they coat plant seeds with to help agricultural products grow better. In order to facilitate research and development Indigo sequenced many new fungi and bacteria. They asked me to create and EDA tool which would allow them to understand what species are similar to each other. Additionally, they were concerned that different subsets of microbes would produce different groupings, and so I needed to be able to perform my analysis online. I used TSNE and DBSCAN to visualize and then cluster these bacteria.

## Data

### PicrustDB_Named
Because of the sensitive nature of Indigo's data, I instead used a placeholder dataset from a project called PICRUSt. This data consists of copy number repeats of KO orthologs which I binarized since the number of repeats isn't necessarily reflective of expression. I subset the number of species down to 5000 and the number of orthologs down to 250 for speed purposes (the original purpose of this app was to be part of a 5 minute demonstration. Speed for Indigo is less critical, but speed for the demonstration purposes were important. Since the dataset is a placeholder, we can still use subsets as a POC). The genes were selected such that they were not very sparse within the dataset. While the original dataset used unique OTU identfiers, as per Indigo's request, I replaced those identifiers with names of the lowest level named classificiation group, which are not necessarily unique. Parsed in dataprocessMPL.py

### 99_otu_taxonomy.txt
A list of names and corresponding OTU numbers. Parsed in addNamesLoop.py

### koOrthologs.json
A list of ko orthologs in JSON format with the corresponding names. Parsed in dataprocessMPL.py.

### Ensparsinator.py
Simple code which gets one out of every 10 lines from a very large file.

### addNamesLoop.py
Code which takes the original database format and replace the OTU numbers with names. I do not believe it can be run as-is however.


## Model

### dataprocess.MPL
This runs the flask app and performs the TSNE and DBSCAN analyses. This takes in data from PicrustDB_Named and uses the koOrthologs.json data to formulate information about the gene composition of each species and puts that information on the graph.

## Testing

### tuning.py
The DBSCAN algorithm has a number of parameters which need tuning. This is the code I used to tune the EPS parameter. I use the metric of a silhouette score on the percentage of genes related to metabolism. I chose the bounds to test since lower created too many clusters or too many outliers, and higher produced too few cluters (roughly a 90-10 split, for 2 clusters, even thought this technically improved my silhouette measure).


### silhouette.py
I have no idea what this is.

### tuningTSNE.py
Pretty sure this is garbage and I should throw it out

### consistency.py
I discuss briefly the concerns about speed in PicrustDB_Named section, and because of those concerns I cannot set the stopping rules of my TSNE to achieve maximum consistency and get a fast answer. What I can do is measure how consistent my answers are. I take a number of subsamples and for each subsample I measure how consistently a species is grouped with another species. If they are grouped together >50% of the time then I say that "truth" is that they are grouped together. I then measure the chance that if you belong together you wind up grouped together as well as the probability that you are not grouped together if you don't belong together. Through multiple runs I find this first number winds up being between about 86%-90%, and the second number is roughly 0%.
