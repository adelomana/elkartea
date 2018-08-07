###
### This script (1) builds a table of all available annotation, including GO from BLAST orthologs and (2) performs hypergeometric tests on different functional categories for each comparison
###

import sys

def annotationFileReader():

    '''
    This function reads the transcript names of C5 and builds a dictionary of different IDs.
    '''

    transcriptNames=[]; annotationStone={}; transcriptDescriptions={}
    
    with open(annotationFile,'r') as f:
        next(f)
        next(f)
        for line in f:
            v=line.split('\t')
            field=v[8].replace('\n','')
            transcriptName=field.split('transcript_id=')[1].split(';')[0]
            ID=field.split('ID=')[1].split(';')[0]
            try:
                product=field.split('product=')[1].split(';')[0]
            except:
                product=''
            
            transcriptNames.append(transcriptName)
            annotationStone[ID]=transcriptName

            transcriptDescriptions[transcriptName]=[ID,product]

    return transcriptNames,annotationStone,transcriptDescriptions

def COGreader():

    '''
    This function reads the provided COG annotation.
    '''

    COGstone={}
    with open(COGfile,'r') as f:
        next(f)
        for line in f:
            v=line.split()

            ID=annotationStone[v[0]]
            COGterm=v[9]
            COGstone[ID]=COGterm

    return COGstone

def GOreader():

    '''
    This function gets GO annotation through MicrobesOnline KIP7 and BLAST orthology.
    '''

    # build dictionary of KIP7 and GO terms
    KIP7GOstone={}
    with open(GOfile,'r') as f:
        next(f)
        for line in f:
            v=line.split('\t')
            ID=v[0]
            GOterm=v[-3]
            
            KIP7GOstone[ID]=GOterm
            
    # build orthology
    orthologyStone={}
    with open(orthologyFile,'r') as f:
        next(f)
        for line in f:
            v=line.split('\t')
            VIMSS=v[0].replace('VIMSS','')
            ga=v[1]
            orthologyStone[ga]=VIMSS

    # map KIP7 to C5 through orthology
    Gaterms=list(set(list(orthologyStone.keys())))
    
    GOstone={}
    for term in Gaterms:
        GOstone[term]=KIP7GOstone[orthologyStone[term]]

    return GOstone,orthologyStone

def KEGGreader():

    '''
    This function reads the provided KEGG annotation.
    '''

    KEGGstone={}
    with open(KEGGfile,'r') as f:
        next(f)
        for line in f:
            v=line.split()

            ID=annotationStone[v[0]]
            KEGGterm=v[9]
            KEGGstone[ID]=KEGGterm

    return KEGGstone

def PFAMreader():

    '''
    This function reads the provided PFAM annotation.
    '''

    PFAMstone={}
    with open(PFAMfile,'r') as f:
        next(f)
        for line in f:
            v=line.split()

            ID=annotationStone[v[0]]
            PFAMterm=v[8]
            PFAMstone[ID]=PFAMterm

    return PFAMstone

def TIGRreader():

    '''
    This function reads the provided TIGR annotation.
    '''

    TIGRstone={}
    with open(TIGRfile,'r') as f:
        next(f)
        for line in f:
            v=line.split()

            ID=annotationStone[v[0]]
            TIGRterm=v[6]
            TIGRstone[ID]=TIGRterm

    return TIGRstone

# 0. user defined variables
print('defining variables...')

genomeDir='/Volumes/omics4tb/alomana/projects/ENIGMA/data/annotations/C5/genome/'
annotationFile=genomeDir+'2767802313.gff'
KEGGfile=genomeDir+'2767802313.ko.tab.txt'
COGfile=genomeDir+'2767802313.cog.tab.txt'
PFAMfile=genomeDir+'2767802313.pfam.tab.txt'
TIGRfile=genomeDir+'2767802313.tigrfam.tab.txt'
GOfile='/Volumes/omics4tb/alomana/projects/ENIGMA/data/annotations/microbesOnline/genomeInfo.txt'
orthologyFile='/Volumes/omics4tb/alomana/projects/ENIGMA/data/annotations/microbesOnline/reciprocal.blast.microbesOnline.C5.output.txt'

functionalAnnotationFile='/Volumes/omics4tb/alomana/projects/ENIGMA/data/annotations/C5/functional/C5.annotation.txt'

# 1. read files
print('reading files...')

# 1.1. read annotation
print('\t reading gene C5 annotation...')
transcriptNames,annotationStone,transcriptDescriptions=annotationFileReader()

# 1.2. read GO annotation
print('\t reading GO annotation...')
GOstone,orthologyStone=GOreader()

# 1.3. read KEGG annotation
print('\t reading KEGG annottion...')
KEGGstone=KEGGreader()

# 1.4. read PFAM annotation
print('\t reading PFAM annottion...')
PFAMstone=PFAMreader()

# 1.5. read COG annotation
print('\t reading COG annottion...')
COGstone=COGreader()

# 1.6. read TIGR annotation
print('\t reading TIGR annottion...')
TIGRstone=TIGRreader()

# 2. build a table of annotation
print('writing functional annotation file...')

g=open(functionalAnnotationFile,'w')
g.write('#C5.transcript.name\tC5.transcript.ID\tC5.transcript.description\tKIP7.ortholog\tGO\tKEGG\tPFAM\tCOG\tTIGR\n')

for transcriptName in transcriptNames:

    # write name
    line2Write='{}\t{}\t{}'.format(transcriptName,transcriptDescriptions[transcriptName][0],transcriptDescriptions[transcriptName][1])
    g.write(line2Write)

    # write GO
    if transcriptName in GOstone:
        g.write('\t{}'.format(orthologyStone[transcriptName]))
        g.write('\t{}'.format(GOstone[transcriptName]))
    else:
        g.write('\t ')
        g.write('\t ')

    # write KEGG
    if transcriptName in KEGGstone:
        g.write('\t{}'.format(KEGGstone[transcriptName]))
    else:
        g.write('\t ')

    # write PFAM
    if transcriptName in PFAMstone:
        g.write('\t{}'.format(PFAMstone[transcriptName]))
    else:
        g.write('\t ')

    # write COG
    if transcriptName in COGstone:
        g.write('\t{}'.format(COGstone[transcriptName]))
    else:
        g.write('\t ')

    # write TIGR
    if transcriptName in TIGRstone:
        g.write('\t{}'.format(TIGRstone[transcriptName]))
    else:
        g.write('\t ')

    # close the line
    g.write('\n')

g.close()

sys.exit()

# 3. perform functional enrichment on comparisons
print('performing hypergeometric test...')

functionalCategories=['GO','KEGG','PFAM','COG','TIGR']

for functionalCategory in functionalCategories:
    numberOfTests=0
    for group in comparisons:
        for term in terms:
            # compute background for that particular term
            # run hypergeometric test
            pval = hypergeom.sf(k-1, M, n, N) # https://blog.alexlenail.me/understanding-and-implementing-the-hypergeometric-test-in-python-a7db688a7458, make sure I get the same results as in R
            numberOfTests=numberOfTests+1
    # correct by Bonferroni

    # report

# dont forget to Bonferroni correct for number of tests: number of groups times category.
# restrict for at least 3 genes.
