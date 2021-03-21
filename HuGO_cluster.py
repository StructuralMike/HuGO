import os
import csv
import time
from pandas import read_csv
import subprocess
import csv
from HuGO import curate_fasta
from HuGO import source_key
from HuGO import HuGOvars

def cdhit(f_input, f_output, cutoff):
    cmd_cdhit = [
            'cd-hit',
            '-i', f_input,
            '-o', f_output,
            '-c', '0.' + cutoff,
            '-M', '8000',
            '-T', '0',
            '-sc', '1'
            ]
    output,error = subprocess.Popen(cmd_cdhit, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return(output,error)

def analyseCdhit(f_analyse):
    hitdict = {}
    accessionCluster = {}
    columnNames = [i for i in range(0, 6)]
    data = read_csv(f_analyse, sep=' ', names=columnNames, header=None, low_memory=False)
    inputdata = data.values.tolist()

    for n in inputdata:
        if n[0][0] == '>':
            idn = n[1]
            hitdict[idn] = []
        else:
            val1 = n[1][1:-3]
            if n[2] != '*':
                val2 = n[3][0:-1]
            else:
                val2 = '100.00'
            hitdict[idn].append((val1,val2))
            accessionCluster[val1] = (idn,val2)
    return(hitdict,accessionCluster)

def userSet(f_userlist):
    try:
        userset = set(open(f_userlist).read().split())
    except:
        userset = set()
    return userset

def criteria(row,userset):
    try:
        len(row[24])    # Does a 3D structure exist
        return True
    except:
        if row[4] == 'reviewed':
            return True
        elif row[2] in userset:
            return True
    return False

def clusterMetaData(f_meta_data_curated,cutoff,f_userlist,hitdict,accessionCluster,accessionColumn=2):
    trackDict = {}
    for i in hitdict:
        trackDict[i] = []

    data = read_csv(f_meta_data_curated, sep='\t', header=None,low_memory=False)          # Pandas automatically sets the correct types for each column
    metaData = data.values.tolist()
    metaDataRed = [metaData[0]]
    userset = userSet(f_userlist)

    #Include obligatory entries
    for row in metaData:
        if row[0] != 'Name':
            if criteria(row,userset):
                metaDataRed.append(row)
                try:
                    idn = accessionCluster[row[accessionColumn]][0]
                    del trackDict[idn]
                except:
                    pass
            else:
                try:
                    idn = accessionCluster[row[accessionColumn]][0]
                    if idn in trackDict:
                        trackDict[idn].append(row)
                    else:
                        pass
                except:
                    pass

    #Include representative from clusters not yet accounted for
    for cluster in hitdict:
        if cluster in trackDict:            
            for entry in trackDict[cluster]:
                if hitdict[cluster][0][0] == entry[accessionColumn]:
                    metaDataRed.append(trackDict[cluster][0])

    f_cdhit_meta_data = os.path.splitext(f_meta_data_curated)[0] + '_clustered_' + cutoff + '.tsv'
    with open(f_cdhit_meta_data, 'w+', newline='') as o:
        wr = csv.writer(o,dialect='excel',delimiter='\t')
        wr.writerows(metaDataRed)
    return(metaDataRed)

def main(HuGO=HuGOvars(full=False)):
    if HuGO['full']:
        f_input = HuGO['files']['curated_fasta_full']
    else:
        f_input = HuGO['files']['curated_fasta']

    cutoff = HuGO['pid']

    print('Running HuGO Cluster {0}% using {1}'.format(cutoff, f_input))

    f_output = HuGO['paths']['cdhit'] + os.path.splitext(os.path.basename(f_input))[0] + '_clustered_' + cutoff + '.fasta'
    output,error = cdhit(f_input, f_output, cutoff)
    f_clusterlog = HuGO['paths']['cdhit'] + 'cdhit.log'
    with open(f_clusterlog, 'w+', newline='') as o:
        o.write(output)

    f_analyse = f_output + '.clstr'
    hitdict,accessionCluster = analyseCdhit(f_analyse)
    f_clustered = os.path.splitext(f_output)[0] + '.tsv'

    with open(f_clustered, 'w+', newline='') as o:
        for n in hitdict:
            for acc,pid in hitdict[n]:
                o.write('{}\t{}\t{}\n'.format(n,acc,pid))

    f_userlist = HuGO['paths']['project'] + 'user_supplied_uniprot.txt'
    if HuGO['source'] == 'pfam':
        accessionColumn = 0
    else:
        accessionColumn = 2
    metaDataRed = clusterMetaData(HuGO['files']['meta_data_curated'],cutoff,f_userlist,hitdict,accessionCluster,accessionColumn)

    f_curate_fasta = os.path.splitext(HuGO['files']['curated_fasta'])[0] + '_clustered_' + cutoff + '.fasta'
    f_curate_fasta_full = os.path.splitext(HuGO['files']['curated_fasta_full'])[0] + '_clustered_' + cutoff + '.fasta'
    curate_fasta(metaDataRed,f_curate_fasta,f_curate_fasta_full,source=HuGO['source'],full=HuGO['full'])
    print('HuGO Cluster done.')

if __name__ == "__main__":
    main()
