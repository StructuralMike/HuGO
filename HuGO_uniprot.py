import urllib.parse
import urllib.request
from Bio import SeqIO
from pandas import read_csv
import csv
import platform
from HuGO import curate_fasta
from HuGO import source_key
from HuGO import HuGOvars


# Take a fasta file and return two dicts, names/seqs.
# names dict is given the accession code as key which is extracted from the fasta id based on f_key,
# and each accession code is associated with a list of fasta ids
# seqs dict have each fasta id as a key and one corresponding sequence
def parse_fasta(f_input,f_key):
    names = {}
    seqs = {}
    fasta_sequences = SeqIO.parse(open(f_input),'fasta')
    for entry in fasta_sequences:
        name, sequence = entry.id, str(entry.seq).upper()
        if name.split(f_key)[0] in names:
            names[name.split(f_key)[0]].append(name)
        else:
            names[name.split(f_key)[0]] = [name]
        seqs[name] = sequence
    return(names,seqs)

# Split a list into iterable chunks of size n
def chunks(lst, n):
    size = len(lst)
    for i in range(0,size,n):
        c = i+n
        if c > size:
            c = size
        print('\tFetching meta data, {0}/{1}'.format(c,size))
        yield lst[i:i + n]

# Submit a REST request to UNIPROT with:
# names: a list of uniprot accession codes
# columns: what meta data columns to retrieve
# meta_data_fname: which file to write the output to
def uniprot_rest(names,columns,f_meta_data):
    currentOS = platform.system()
    if currentOS == 'Darwin':
        from HuGO import sslContext
        ssl_context = sslContext()
    all_entries = [*names]    
    cols = ','.join(columns)
    with open(f_meta_data, 'w+', newline='') as o:
        for chunk in chunks(all_entries,200):
            url = 'https://www.uniprot.org/uniprot/'
            query = 'id:('
            query += str(' OR '.join(chunk))
            query += ')'
            params = {
                    'query': query,
                    'format': 'tab',
                    'columns': cols
                    }
            data = urllib.parse.urlencode(params)
            data = data.encode('utf-8')
            req = urllib.request.Request(url, data)
            if currentOS == 'Darwin':
                with urllib.request.urlopen(req,context=ssl_context) as f:
                    response = f.read()
            else:
                with urllib.request.urlopen(req) as f:
                    response = f.read()                
            o.write(response.decode('utf-8'))

# For each fasta id, check the uniprot accession code
# then append the meta data that was fetched from UNIPROT using uniprot_rest()
# and write the output to a new file 
def curateData(f_input,seqs,columns,f_key,f_meta_data_curated):
    print('\tCurating meta data')
    outputdata = [['Name','Sequence']]
    outputdata[0].extend(columns)
    data = read_csv(f_input, sep='\t', header=None)          # Pandas automatically sets the correct types for each column
    inputdata = data.values.tolist()
    for entry in seqs:
        next_append = []
        next_append.extend([entry,seqs[entry]])
        acc = entry.split(f_key)[0]
        for n in inputdata:
            if n[0] == acc and n[8] != 'fragment' and not 'X' in str(n[18]).upper():
                if n[2] == 'reviewed' or n[2] == 'unreviewed':
                    next_append.extend(n)
                    outputdata.append(next_append)
                break
    with open(f_meta_data_curated, 'w+', newline='') as o:
        wr = csv.writer(o,dialect='excel',delimiter='\t')
        wr.writerows(outputdata)
    return(outputdata)

def main(HuGO=HuGOvars()):
    #Input variables 
    files = HuGO['files']
    f_key = source_key(HuGO['source'])                          

    print('Running HuGO Uniprot using {}'.format(files['input_fasta']))
    
    # The meta data columns to fetch from uniprot
    columns = [                                             
        'id',
        'entry name',
        'reviewed',
        'protein names',
        'genes',
        'organism',
        'organism-id',
        'length',
        'fragment',
        'encodedon',
        'database(InterPro)',
        'database(Pfam)',
        'lineage(CLASS)',
        'lineage(FAMILY)',
        'lineage(GENUS)',
        'lineage(KINGDOM)',
        'lineage(ORDER)',
        'lineage(PHYLUM)',
        'sequence',
        'citation',
        'citationmapping',
        'database(RefSeq)',
        '3d',
        'ec',
        'go-id'
        ]

    #Program execution order
    #Parse fasta file, fetch column meta data from uniprot, then remove non-relevant entries and ouput curated meta data, and curated fasta.

    names,seqs = parse_fasta(files['input_fasta'],f_key)
    uniprot_rest(names,columns,files['meta_data'])
    curatedData = curateData(files['meta_data'],seqs,columns,f_key,files['meta_data_curated'])
    curate_fasta(curatedData,files['curated_fasta'],files['curated_fasta_full'],source=HuGO['source'],full=True)
    print('HuGO Uniprot done.')

if __name__ == "__main__":
    main()
