import urllib.parse
import urllib.request
import platform
from HuGO import HuGOvars

def getFamily(family,f_key,f_input_fasta):
    currentOS = platform.system()
    if currentOS == 'Darwin':
        from HuGO import sslContext
        ssl_context = sslContext()
    with open(f_input_fasta, 'w+', newline='') as o:
        if f_key == 'interpro':
            url = 'https://www.ebi.ac.uk/interpro/legacy/entry/{}/proteins-matched'.format(family)
        elif f_key == 'panther':
            url = 'https://www.ebi.ac.uk/interpro/legacy/signature/{}/proteins-matched'.format(family.replace(':','%3A'))
        params = {'export': 'fasta'}
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

def main(HuGO=HuGOvars()):
    source = 'interpro'
    family = 'IPR032358'
    f_input_fasta = './HuGO/DUF4867/{}.fasta'.format(family.replace(':','_'))
    print('Fetching family {} from {}'.format(family,source))
    getFamily(family,source,f_input_fasta)

if __name__ == "__main__":
    main()