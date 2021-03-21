import pandas as pd
import csv

def readEntryTypes(filename):
    data = pd.read_csv(filename,sep='\t', header=1)
    data = data.values.tolist()
    entries = {}
    for entry_ac,entry_type,entry_name in data:
        entries[entry_ac] = {'type': entry_type, 'name': entry_name}
    return(entries)

def defLevel(line):
    level = 0
    for char in line:
        if char == '-':
            level += 1
        else:
            break
    return int(level/2)


interpro = readEntryTypes('./interpro/Types.tsv')
treefile = './interpro/ParentChildTreeFile.txt'
with open(treefile) as t:
    lines = t.read().split('\n')

topname = ''
for row in lines:
    try:
        accession,name = row.split('::')
        level = defLevel(accession)
        accession = accession[level*2:]
        if level == 0:
            topaccession = accession
            topname = name
        else:
            interpro[accession]['level'] = level
            interpro[accession]['topaccession'] = topaccession
            interpro[accession]['topname'] = topname
    except:
        pass

with open('./interpro/DomainAssociations.tsv', 'w+', newline='') as o:
    for key in interpro:
        try:
            level = interpro[key]['level']
            topaccession = interpro[key]['topaccession']
            topname = interpro[key]['topname']
        except:
            level = ''
            topaccession = ''
            topname = '' 
        line = '{}\t{}\t{}\t{}\t{}\t{}\n'.format(key,interpro[key]['type'],interpro[key]['name'],level,topaccession,topname)
        o.write(line)