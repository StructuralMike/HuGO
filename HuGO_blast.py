import os
import csv
import time
import pandas as pd
import subprocess
import numpy as np
import math
from HuGO import HuGOvars

def dataCleanup(inputdata):
    unique_rows = {}
    for n in inputdata:
        # dict.get return None if key doesn't exist
        if n[0] != n[1]:
            row_key = create_key(n)
            if unique_rows.get(row_key) is None:
                newE = (n[-1]/(2**(n[3]/2)))       # A new E-value is by dividing the bitscore by an arbitrary value to get around the round-down-to-0-problem that occurs at E-180 in the blastp algorithm and at E-324 in python
                unique_rows[row_key] = n[0],n[1],newE
    return [row for key, row in unique_rows.items()]

def create_key(row):
#    key_vals = row[0], row[1], row[9], row[8], row[11], row[10]
    key_vals = row[0], row[1]
    key_vals_sorted = sorted([str(val) for val in key_vals])
    return "_".join(key_vals_sorted)

def processBlastOutput(f_blastpout,f_blastpreduced):
    data = pd.read_csv(f_blastpout, sep=',', header=None)
    inputdata = data.values.tolist()

    print('\tSorting blast results...')
    inputdata.sort(key=lambda x: x[2])

    print('\tRemoving redundancy...')
    outputdata = dataCleanup(inputdata)

    print('\tSorting output data...')
    outputdata.sort(key=lambda x: x[2], reverse=True)

    print('\tWriting to file...')
    with open(f_blastpreduced, 'w+', newline='') as o:
        wr = csv.writer(o,dialect='excel')
        wr.writerows(outputdata)

def allvsall_blastp(blast_db_inputfile,blast_db_outputfile,blast_inputfile,blast_outputfile,threads):
    command_line_mkblastdb = [
            'makeblastdb',
            '-dbtype','prot',
            '-in', blast_db_inputfile,
            '-out', blast_db_outputfile,
            '-input_type', 'fasta'
            ]
    subprocess.call(command_line_mkblastdb)

    print('\tRunning all vs all blastp...')
    command_line_blastp = [
            'blastp',
            '-matrix', 'BLOSUM62',
            '-num_threads', threads,
            '-out', blast_outputfile,
            '-query', blast_inputfile,
            '-db', blast_db_outputfile,
            '-outfmt', "10 qseqid sseqid evalue bitscore pident nident mismatch positive gaps qstart qend sstart send length"
        ]
    output,error = subprocess.Popen(command_line_blastp, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return(output,error)

def cutOffs(evalues):
    evalue,efreq = np.unique(evalues,return_counts=True)
    esum = np.sum(efreq)
    ecmlfreq = [n/esum for n in np.cumsum(efreq)]
    lim = 0.1
    cutoffs = []
    for i,y in enumerate(list(ecmlfreq)):
        if y > lim:
            #Check which e-cutoff is closest
            current = abs(y - lim)
            previous = abs(ecmlfreq[i-1] - lim)
            if current <= previous:
                cutoffs.append([round(lim*100),evalue[i]])
            else:
                cutoffs.append([round(lim*100),evalue[i-1]])
            lim += 0.1
    return cutoffs

def magnitude(value):
    if value == 0:
        return np.PINF
    else:
        return -math.floor(math.log10(value))

def fileCutoff(f_blastpreduced):
    print('\tWriting cytoscape input files...')

    data = pd.read_csv(f_blastpreduced, sep=',', header=None)
    inputdata = data.values.tolist()
    evalues = np.array([magnitude(n[2]) for n in inputdata])
    cutoffs = cutOffs(evalues)

    cutoffdict = {}
    for cutoff,e in cutoffs:
        cutoffdict[cutoff] = [e,[]]

    for i,n in enumerate(list(evalues)):
        for key in cutoffdict:
            if n >= cutoffdict[key][0] or np.isposinf(evalues[i]):
                cutoffdict[key][1].append(inputdata[i])

    for key in cutoffdict:
        try:
            cutoff = '-percent_E' + str(-int(cutoffdict[key][0]))
        except:
            cutoff = '-percent_0'
        outputfile = (os.path.splitext(f_blastpreduced)[0] + '_cutoff-' + str(key) + cutoff + '.csv')
        with open(outputfile, 'w+', newline='') as o:
            wr = csv.writer(o,dialect='excel')
            wr.writerows(cutoffdict[key][1])

def main(HuGO=HuGOvars()):
    files = HuGO['files']
    threads = HuGO['threads']
    pid = HuGO['pid']
    if HuGO['full']:
        f_input = files['curated_fasta_full_clustered']
    else:
        f_input = files['curated_fasta_clustered']
    f_db = HuGO['paths']['blastp'] + os.path.splitext(os.path.basename(f_input))[0]
    f_query = f_input
    f_blastout = HuGO['paths']['blastp'] + 'blastp_allvsall_' + pid + '.csv'

    print('Running HuGO blast using {}'.format(f_input))

    output,error = allvsall_blastp(f_input,f_db,f_query,f_blastout,threads)
    f_blastplog = HuGO['paths']['blastp'] + 'blastp.log'
    with open(f_blastplog, 'w+', newline='') as o:
        o.write(output)

    f_blastreduced = HuGO['paths']['cytoscape'] + os.path.basename(f_blastout)
    processBlastOutput(f_blastout,f_blastreduced)
    fileCutoff(f_blastreduced)

    print('HuGO blast done.')

if __name__ == "__main__":
    main()
