import pandas as pd
import csv

def getRefseqs(metaDataFile,databasePriority):
    data = pd.read_csv(metaDataFile, sep='\t', header=None)
    metaData = data.values.tolist()
    output = []
    #Find column for refseqs
    for i,col in enumerate(metaData[0]):
        if col == 'database(RefSeq)':
            refseqCol = i
    for i,entry in enumerate(metaData[1:]):
        name = ''
        found = 0
        try:
            refseqs = entry[refseqCol].split(';')
            if len(refseqs) == 1:
                name = refseqs[0]
            else:
                for database in databasePriority:
                    for refseq in refseqs:
                        if refseq[0:1] == database:
                            name = refseq
                            found = 1
                            break
                if found == 0:
                    name = refseqs[0]
        except:
                name = 'nan'
        if len(name) > 5:
            output.append(name)
    return(output)

path = 'rubisco/'
metaDataFile = path + 'meta_data_curated.tab'
databasePriority = ['WP','YP','NP','XP','AP']
refseqList = getRefseqs(metaDataFile,databasePriority)
outputfile = path + 'refseqList.txt'
with open(outputfile, 'w+', newline='') as o:
    for refseq in refseqList:
        o.write(refseq + '\n')
print('Wrote {0} refseqs to {1}'.format(len(refseqList),outputfile))