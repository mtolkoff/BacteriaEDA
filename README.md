# COloNY Description
This was a project I created with Insight Data Science in consultation with Indigo Agriculture. Indigo is a company which makes "seed coatings," or bacteria and fungi which they coat plant seeds with to help agricultural products grow better. In order to facilitate research and development Indigo sequenced many new fungi and bacteria. They asked me to create and EDA tool which would allow them to understand what species are similar to each other. Additionally, they were concerned that different subsets of microbes would produce different groupings, and so I needed to be able to perform my analysis online. I used TSNE and DBSCAN to visualize and then cluster these bacteria.

## Data

### PicrustDB_Named
Because of the sensitive nature of Indigo's data, I instead used a placeholder dataset from a project called PICRUSt. This data consists of copy number repeats of KO orthologs which I binarized since the number of repeats isn't necessarily reflective of expression. I subset the number of species down to 5000 and the number of orthologs down to 250 for speed purposes (the original purpose of this app was to be part of a 5 minute demonstration. Speed for Indigo is less critical, but speed for the demonstration purposes were important. Since the dataset is a placeholder, we can still use subsets as a POC). The genes were selected such that they were not very sparse within the dataset. While the original dataset used unique OTU identfiers, as per Indigo's request, I replaced those identifiers with names of the lowest level named classificiation group, which are not necessarily unique. Parsed in dataprocessMPL.py

## 99_otu_taxonomy.txt
A list of names and corresponding OTU numbers. Parsed in addNamesLoop.py

## koOrthologs.json
A list of ko orthologs in JSON format with the corresponding names. Parsed in dataprocessMPL.py.

### Ensparsinator.py

### addNamesLoop.py

## Model

### dataprocess.MPL

## Testing

### tuning.py

### silhouette.py

### tuningTSNE.py

### consistency.py
