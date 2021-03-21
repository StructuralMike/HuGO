import os
import csv
import time
import pandas as pd
import subprocess
from HuGO import HuGOvars

def interproscan(path, f_input, threads):
    cmd_interproscan = [
            'interproscan.sh',
            '--input',f_input,
            '--formats', 'tsv',
            '--output-file-base', path + 'iprscan5_' + os.path.splitext(os.path.basename(f_input))[0],
            '--cpu', threads,
            '--seqtype', 'p',
            '--tempdir', path + 'temp',
            '--verbose'
            ]

    output,error = subprocess.Popen(cmd_interproscan, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    return(output,error)

def main(HuGO = HuGOvars()):
    files = HuGO['files']
    if HuGO['full']:
        f_input = files['curated_fasta_full_clustered']
    else:
        f_input = files['curated_fasta_clustered']

    print('Running HuGO interproscan using {}'.format(f_input))

    output,error = interproscan(HuGO['paths']['interproscan'], f_input, HuGO['threads'])
    f_interprolog = HuGO['paths']['interproscan'] + 'interproscan.log'
    with open(f_interprolog, 'w+', newline='') as o:
        o.write(output)

    print('HuGO interproscan done.')

if __name__ == "__main__":
    main()