# Mac OS needs ssl certificates specified
def sslContext():
    import certifi
    import ssl
    return ssl.create_default_context(cafile=certifi.where())

# The key separating accession code from the rest of the fasta id
def source_key(source):
    if source == 'interpro':
        f_key = '|'
    elif source == 'panther':
        f_key = '|'
    elif source == 'pfam':
        f_key = '/'
    elif source == 'pfam2':
        f_key = '_'
    return(f_key)

# Generate fasta files from a meta data file
def curate_fasta(curatedData,f_curated_fasta,f_curated_fasta_full,source,full=True):
    curated_fasta = []
    sequenceColumn = 1
    if source == 'pfam':
        accessionColumn = 0
    else:
        accessionColumn = 2
    for i,n in enumerate(curatedData):
        if i > 0:
            curated_fasta.append('>' + str(n[accessionColumn]))
            curated_fasta.append(str(n[sequenceColumn]))
    with open(f_curated_fasta, 'w+') as o:
        for item in curated_fasta:
            o.write(str(item) + '\n')

    if not full:
        curated_fasta = {}
        accessionColumn = 2
        sequenceColumn = 20
        for i,n in enumerate(curatedData):
            if i > 0:
                curated_fasta[n[accessionColumn]] = n[sequenceColumn]

        with open(f_curated_fasta_full, 'w+') as o:
            for key in curated_fasta:
                o.write('>{}\n{}\n'.format(key,curated_fasta[key]))
    print('\tWrote curated fasta to: {}'.format(f_curated_fasta_full))

def HuGOvars(filename='IPR032358.fasta',project='DUF4867',fastaformat='interpro',maxthreads='14',percent_identity='99red',full=True):
    import os
    HuGO = {}

    source = fastaformat
    threads = maxthreads
    pid = percent_identity

    p_project = 'HuGO/' + project.replace(' ','_') + '/'
    if not os.path.exists(p_project[0:-1]):
        try:
            os.makedirs(p_project[0:-1])
        except:
            raise 'Invalid characters in project name'

    workdirs = ['metadata','fasta','blastp','cdhit','cytoscape','clustalo','interproscan','itol']
    paths = {'project': p_project}
    for d in workdirs:
        paths[d] = paths['project'] + d + '/'
        if not os.path.exists(p_project + d):
            try:
                os.makedirs(paths['project']+d)
            except:
                print('\tCould not create directory {}'.format(paths['project']+d))


    files = {
        'input_fasta': p_project + filename,
        'meta_data': paths['metadata'] + 'meta_data.tsv',
        'meta_data_curated': paths['metadata'] + 'meta_data_curated.tsv',
        'curated_fasta': paths['fasta'] + project + '_curated.fasta',
        'curated_fasta_full': paths['fasta'] + project + '_curated_full.fasta',
        'itoldomains': paths['itol'] + project + 'iToldomains.txt'
    }
    files['curated_fasta_clustered'] = os.path.splitext(files['curated_fasta'])[0] + '_clustered_' + pid + '.fasta'
    files['curated_fasta_full_clustered'] = os.path.splitext(files['curated_fasta_full'])[0] + '_clustered_' + pid + '.fasta'

    HuGO['source'] = source
    HuGO['threads'] = threads
    HuGO['pid'] = pid
    HuGO['files'] = files
    HuGO['paths'] = paths
    HuGO['full'] = full
    HuGO['project'] = str(project)
    return HuGO

def main():
    import HuGO_getFamily
#    HuGO_getFamily.main()

    import HuGO_uniprot
#    HuGO_uniprot.main()

    import HuGO_cluster
#    HuGO_cluster.main()

    import HuGO_blast
    HuGO_blast.main()

    import HuGO_interproscan
#    HuGO_interproscan.main(HuGO = HuGOvars(full=True))

    import HuGO_itoldomains
#    HuGO_itoldomains.main(HuGO = HuGOvars(full=True))

    import HuGO_tree
#    HuGO_tree.main()

if __name__ == "__main__":
    main()