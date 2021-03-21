import os
import csv
import time
import pandas as pd
import subprocess
from HuGO import HuGOvars

def clustalo(path, f_input):
    cmd_clustalo = [
            'clustalo',
            '--infile',f_input,
            '--MAC-RAM', '8000',
            '--verbose',
            '--guidetree-out', path + os.path.splitext(os.path.basename(f_input))[0] + '_tree.tree',
            '--outfmt', 'clustal',
            '--resno',
            '--outfile', path + os.path.splitext(os.path.basename(f_input))[0] + '_alignment.clustal',
            '--output-order', 'tree-order',
            '--seqtype', 'protein',
            '--force',
            '--distmat-out=' + path + os.path.splitext(os.path.basename(f_input))[0] + '_matrix.pim',
            '--percent-id',
            '--full',
#            '--hmm-in=HuGO/gbpc/GbpC.hmm',
#            '--max-hmm-iterations=-1'
            ]

    output,error = subprocess.Popen(cmd_clustalo, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return(output,error)

def main(HuGO = HuGOvars()):
    files = HuGO['files']
    if HuGO['full']:
        f_input = files['curated_fasta_full_clustered']
    else:
        f_input = files['curated_fasta_clustered']

    print('Running HuGO tree using {}'.format(f_input))
    output,error = clustalo(HuGO['paths']['clustalo'], f_input)
    print(error)
    f_treelog = HuGO['paths']['clustalo'] + 'clustalo.log'
    with open(f_treelog, 'w+', newline='') as o:
        o.write(output)
    
    print('HuGO tree done.')

if __name__ == "__main__":
    main()