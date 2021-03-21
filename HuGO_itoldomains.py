import pandas as pd
import os
from HuGO import HuGOvars

class Colors:
    def __init__(self,colors):
        self.colors = list(colors.keys())
        self.size = len(self.colors)
        self.used = 0
 
    def __next__(self):
        if self.used >= self.size:
            raise StopIteration
        else:
            self.used += 1
            return self.colors[self.used-1]

def constructOutput(name,length,foundDomains,domainAttributes,autoColors,toplevels,customAssociations):
    exclusions = {'IPR018247','IPR012706','IPR004116','IPR027267','G3DSA:3.10.20.320','G3DSA:3.10.20.470','G3DSA:3.10.20.460','G3DSA:1.20.5.250','IPR032300','IPR026345','IPR002048'}
    construct = f'{name},{length}'

    for start,end,domainID,domainName in foundDomains:
        if domainID not in exclusions:
            if domainID in customAssociations:
                domainID = customAssociations[domainID]
                domainName = toplevels[domainID]
            try:
                entry_type = toplevels[domainID]['entry_type']
                level = toplevels[domainID]['level']
            except:
                entry_type = 'domain'
                level = 0

            if level > 0:
                key = toplevels[domainID]['top_key']
                name = toplevels[domainID]['top_name']
            else:
                key = domainID
                name = domainName

            try:
                shape = domainAttributes[key]['SHAPE']
                label = domainAttributes[key]['LABEL']
                color = domainAttributes[key]['COLOR']
            except:
                if entry_type == 'active_site' or entry_type == 'conserved_site' or entry_type == 'binding_site':
                    shape = 'DI'
                elif entry_type == 'domain' or entry_type == 'family' or entry_type == 'homologous_superfamily':
                    shape = 'RE'
                elif entry_type == 'repeat':
                    shape = 'EL'
                else:
                    shape ='RE'
                label = name[0:3]
                try:
                    color = next(autoColors)
                except:
                    color = '#666666'
                domainAttributes[key] = {'NAME': name, 'LABEL': label, 'SHAPE': shape, 'COLOR': color}

            construct += f',{shape}|{start}|{end}|{color}|{label}'

    return construct

def domainAnalysis(foundDomains):
    pass

def proteinDomains(filename):
    columnNames = [i for i in range(0, 15)]
    data = pd.read_csv(filename, names=columnNames, sep='\t', header=None)
    return(data.values.tolist())

def splitIntoProteins(proteinDomains):
    proteins = {}
    for entry in proteinDomains:
        if proteins.get(entry[0]) is None: 
            proteins[entry[0]] = [entry]
        else:
            proteins[entry[0]].append(entry)
    return proteins

def customAss(tuples):
    customAssociations = {}
    for child,parent in tuples:
        customAssociations[child] = parent
    return customAssociations

def associateDomains(proteins,toplevels,domainAttributes,autoColors,customAssociations):
    domainCount = {}
    for protein in proteins:
        for domain in proteins[protein]:
            domainID = domain[11]
            domainName = domain[12]
            try:
                len(domainID)
                if domainID in customAssociations:
                    domainID = customAssociations[domainID]
                    domainName = toplevels[domainID]
                entry_type = toplevels[domainID]['entry_type']
                level = toplevels[domainID]['level']
                if level > 0:
                    key = toplevels[domainID]['top_key']
                    name = toplevels[domainID]['top_name']
                else:
                    key = domainID
                    name = domainName

                if key in domainCount:
                    domainCount[key] += 1
                else:
                    domainCount[key] = 1
            except:
                pass
    listofTuples = sorted(domainCount.items() , reverse=True, key=lambda x: x[1])
    for key,count in listofTuples:
        if key in domainAttributes:
            pass
        else:
            entry_type = toplevels[key]['entry_type']
            name = toplevels[key]['entry_name']
            if entry_type == 'active_site' or entry_type == 'conserved_site' or entry_type == 'binding_site':
                shape = 'DI'
            elif entry_type == 'domain' or entry_type == 'family' or entry_type == 'homologous_superfamily':
                shape = 'RE'
            elif entry_type == 'repeat':
                shape = 'EL'
            else:
                shape ='RE'
            label = name[0:3]
            try:
                color = next(autoColors)
            except:
                color = '#666666'
            domainAttributes[key] = {'NAME': name, 'LABEL': label, 'SHAPE': shape, 'COLOR': color}
        domainAttributes[key]['OCCURANCES'] = count

    return domainAttributes

def iTolDomains(proteinDomains, domainAttributes, autoColors, toplevels, customAssociations):

    proteinNameCol = 0 
    protLenCol = 2
    dbCol = 3
    idCol = 4
    nameCol = 5
    startCol = 6
    endCol = 7
    interproCol = 11
    interproNameCol = 12
    name = proteinDomains[0][proteinNameCol]
    length = proteinDomains[0][protLenCol]
    proteinDomains.sort(key=lambda x: x[startCol])
#    validDomains = ['MobiDBLite','ProSiteProfiles','Coils']
    validDomains = []

    foundDomains = []
    previousEnd = 0
    for domain in proteinDomains:
        try:
            len(domain[interproNameCol])
            if domain[startCol] >= previousEnd:
                foundDomains.append((domain[startCol],domain[endCol],domain[interproCol],domain[interproNameCol]))
                previousEnd = domain[endCol]
        except:
            if domain[dbCol] in validDomains:
                foundDomains.append((domain[startCol],domain[endCol],domain[idCol],domain[nameCol]))

    return(constructOutput(name,length,foundDomains,domainAttributes,autoColors,toplevels,customAssociations))

def uniqueColors(domainAttributes,colors):
    for entry in domainAttributes:
        if colors.get(domainAttributes[entry]['COLOR']) is None:
            pass
        else:
            del colors[domainAttributes[entry]['COLOR']]
    return colors

def domainAssociations(filename='HuGO/HuGOresources/interpro/DomainAssociations.tsv'):
    data = pd.read_csv(filename,sep='\t', header=None, engine='python')
    data = data.values.tolist()
    toplevels = {}
    for key,entry_type,entry_name,level,top_key,top_name in data:
        toplevels[key] = {'entry_type': entry_type, 'entry_name': entry_name, 'level': level, 'top_key': top_key, 'top_name': top_name}
    return toplevels

def main(HuGO = HuGOvars()):
    
    toplevels = domainAssociations()
    domainAttributes = {
    #            'IPR022263': {'NAME': 'KxYKxGKxW signal peptide', 'LABEL': 'KxY', 'SHAPE': 'TR', 'COLOR': '#6a3d9a'},
    #            'IPR021197': {'NAME': 'Cross-wall-targeting lipoprotein motif', 'LABEL': 'CWT', 'SHAPE': 'TR', 'COLOR': '#a6cee3'},
    #            'IPR005877': {'NAME': 'YSIRK Gram-positive signal peptide', 'LABEL': 'YRK', 'SHAPE': 'TR', 'COLOR': '#cab2d6'},
    #            'IPR041324': {'NAME': 'Antigen I/II N-terminal', 'LABEL': 'ANT', 'SHAPE': 'RE', 'COLOR': '#e31a1c'},
    #            'IPR009578': {'NAME': 'Cell surface antigen repeat', 'LABEL': 'CAR', 'SHAPE': 'RE', 'COLOR': '#1f78b4'},
    #            'IPR036234': {'NAME': 'Glucan-binding protein C/Surface antigen I/II, V-domain', 'LABEL': 'GBP', 'SHAPE': 'RE', 'COLOR': '#33a02c'},
    #            'IPR013574': {'NAME': 'Glucan-binding protein C/Surface antigen I/II, V-domain', 'LABEL': 'GBP', 'SHAPE': 'RE', 'COLOR': '#33a02c'},
    #            'IPR026345': {'NAME': 'Adhesin isopeptide-forming adherence domain', 'LABEL': 'CSA', 'SHAPE': 'RE', 'COLOR': '#fdbf6f'},
    #            'IPR032300': {'NAME': 'Cell surface antigen, C-terminal', 'LABEL': 'CSC', 'SHAPE': 'RE', 'COLOR': '#ff7f00'},
    #            'IPR019948': {'NAME': 'Gram-positive LPXTG cell wall anchor', 'LABEL': 'LPX', 'SHAPE': 'TL', 'COLOR': '#fb9a99'},
    #            'IPR019931': {'NAME': 'Gram-positive LPXTG cell wall anchor', 'LABEL': 'LPX', 'SHAPE': 'TL', 'COLOR': '#fb9a99'},
    #            'IPR027578': {'NAME': 'Proline-rich tail region repeat', 'LABEL': 'PRR', 'SHAPE': 'RE', 'COLOR': '#b15928'},

    #            'IPR009459': {'NAME': 'MucBP domain', 'LABEL': 'MUC', 'SHAPE': 'RE', 'COLOR': '#b2df8a'},
    #            'IPR041558': {'NAME': 'Mucin binding domain', 'LABEL': 'MUC', 'SHAPE': 'RE', 'COLOR': '#b2df8a'},
    #            'IPR041495': {'NAME': 'Mub B2-like domain', 'LABEL': 'MUB2', 'SHAPE': 'RE', 'COLOR': '#b2df8a'},
    #            'IPR026466': {'NAME': 'Fimbrial isopeptide formation D2 domain', 'LABEL': 'FIM', 'SHAPE': 'RE', 'COLOR': '#ffff99'},
    #            'IPR012569': {'NAME': 'Leucine-rich repeat-containing adjacent domain', 'LABEL': 'LRA', 'SHAPE': 'RE', 'COLOR': '#f1b6da'},
    #            'IPR032675': {'NAME': 'Leucine-rich repeat domain superfamily', 'LABEL': 'LRR', 'SHAPE': 'RE', 'COLOR': '#c51b7d'},
    #            'IPR008160': {'NAME': 'Collagen triple helix repeat', 'LABEL': 'CTH', 'SHAPE': 'RE', 'COLOR': '#276419'},
    #            'IPR019734': {'NAME': 'Tetratricopeptide repeat', 'LABEL': 'TTP', 'SHAPE': 'RE', 'COLOR': '#dfc27d'}
                }

    workdir = HuGO['paths']['interproscan']
    files = []
    for file in os.listdir(HuGO['paths']['interproscan']):
        if file.startswith("iprscan5") and (file.endswith(".txt") or file.endswith(".tsv")):
            files.append(os.path.join(workdir,file))

    print('Running HuGO iToldomains using\n{}\nand {} predefined domains.'.format(files,len(domainAttributes)))

    colors = {
        '#a6cee3': 0,'#1f78b4': 0,'#b2df8a': 0,'#33a02c': 0,'#fb9a99': 0,'#e31a1c': 0,'#fdbf6f': 0,'#ff7f00': 0,'#cab2d6': 0,'#6a3d9a': 0,
        '#ffff99': 0,'#b15928': 0,'#8dd3c7': 0,'#ffffb3': 0,'#bebada': 0,'#fb8072': 0,'#80b1d3': 0,'#fdb462': 0,'#b3de69': 0,'#fccde5': 0,
        '#d9d9d9': 0,'#bc80bd': 0,'#ccebc5': 0,'#ffed6f': 0,'#543005': 0,'#8c510a': 0,'#bf812d': 0,'#dfc27d': 0,'#f6e8c3': 0,'#f5f5f5': 0,
        '#c7eae5': 0,'#80cdc1': 0,'#35978f': 0,'#01665e': 0,'#003c30': 0,'#8e0152': 0,'#c51b7d': 0,'#de77ae': 0,'#f1b6da': 0,'#fde0ef': 0,
        '#f7f7f7': 0,'#e6f5d0': 0,'#b8e186': 0,'#7fbc41': 0,'#4d9221': 0,'#276419': 0,'#9e0142': 0,'#d53e4f': 0,'#f46d43': 0,'#fdae61': 0,
        '#fee08b': 0,'#ffffbf': 0,'#e6f598': 0,'#abdda4': 0,'#66c2a5': 0,'#3288bd': 0,'#5e4fa2': 0,'#ffffd9': 0,'#edf8b1': 0,'#c7e9b4': 0,
        '#7fcdbb': 0,'#41b6c4': 0,'#1d91c0': 0,'#225ea8': 0,'#253494': 0,'#081d58': 0
        }
    colors = uniqueColors(domainAttributes,colors)
    autoColors = Colors(colors)
    itolRow = []
    proteins = {}
    for f in files:
        proteins = {**proteins, **splitIntoProteins(proteinDomains(f))}

#    tuples = []
    tuples = [('IPR013574','IPR036234'),('IPR019948','IPR019931'),('IPR026345','IPR032300'),('IPR002048','IPR032300'),('PR01217','IPR027578')]
    customAssociations = customAss(tuples)
    domainAttributes = associateDomains(proteins,toplevels,domainAttributes,autoColors,customAssociations)
    
    for protein in proteins:
        itolRow.append(iTolDomains(proteins[protein],domainAttributes,autoColors,toplevels, customAssociations) + '\n')

    legendShapes = 'LEGEND_SHAPES'
    legendColors = 'LEGEND_COLORS'
    legendLabels = 'LEGEND_LABELS'
    for entry in domainAttributes:
        try:
            NAME,SHAPE,COLOR,OCCURANCES = str.replace(domainAttributes[entry]['NAME'].strip(),',',';'),domainAttributes[entry]['SHAPE'],domainAttributes[entry]['COLOR'],domainAttributes[entry]['OCCURANCES']
        except:
            NAME,SHAPE,COLOR,OCCURANCES = str.replace(domainAttributes[entry]['NAME'].strip(),',',';'),domainAttributes[entry]['SHAPE'],domainAttributes[entry]['COLOR'],'nd'
        legendShapes += f',{SHAPE}'
        legendColors += f',{COLOR}'
        legendLabels += f',[{entry}] {NAME} (#{OCCURANCES})'

    outputfile = HuGO['files']['itoldomains']
    with open(outputfile, 'w+', newline='') as o:
        o.write('DATASET_DOMAINS\n')
        o.write('SEPARATOR COMMA\n')
        o.write('DATASET_LABEL,Domains\n')
        o.write('COLOR,#000000\n')
        o.write('LEGEND_TITLE,' + HuGO['project'] + '\n')
        o.write(legendShapes + '\n')
        o.write(legendColors + '\n')
        o.write(legendLabels + '\n')
        o.write('SHOW_DOMAIN_LABELS,0\n')
        o.write('LABELS_ON_TOP,0\n')
        o.write('BACKBONE_COLOR,#000000\n')
        o.write('BORDER_WIDTH,1\n')
        o.write('BORDER_COLOR,#000000\n')
        o.write('DATA\n')
        for row in itolRow:
            o.write(row)
    
    print('HuGO iToldomains done.')

if __name__ == "__main__":
    main()
