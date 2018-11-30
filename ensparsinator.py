count = 0

results = open('picrustDB_thinned2.csv', 'w')

with open('picrustDB_massaged.csv') as f:
    for line in f:
        if count % 10 == 0:
            results.write(line)
        count = (count + 1) % 10

results.close()
f.close()
