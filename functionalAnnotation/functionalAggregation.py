###
### This script (1) builds a table of all available annotation, including GO from BLAST orthologs.
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

def mo2refseqReader():

    '''
    This function reads the BLAST orthologs from KIP7 Microbes Online and RefSeq.
    '''

    # first get the map from WP to INTCA
    stone={}
    with open(refSeqMultiFASTA,'r') as f:
        for line in f:
            if line[0] == '>':
                v=line.split()
                wp=v[0].replace('>','')
                intcaID=line.split('locus_tag=')[1].split(']')[0]
                old_intcaID=line                
                stone[wp]=intcaID

    # now we are ready for ortholog mapping
    intcaStone={}
    with open(mo2refseqOrthologyFile,'r') as f:
        next(f)
        for line in f:
            v=line.split('\t')
            intcaStone[v[0]]=stone[v[1]]

    return intcaStone
    
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

def refseqReader():

    '''
    This function maps old and new INTCA IDs.
    '''

    newoldINTCAStone={}

    with open(refSeqGFF3,'r') as f:
       next(f)
       next(f)
       next(f)
       for line in f:
           v=line.split('\t')
           if v[2] == 'gene':
               field=v[8]
               if 'old_locus_tag=' in field:
                   new=field.split('locus_tag=')[1].split(';')[0]
                   old=field.split('old_locus_tag=')[1].split(';')[0].replace('\n','')
                   newoldINTCAStone[new]=old

    return newoldINTCAStone

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
refSeqMultiFASTA='/Volumes/omics4tb/alomana/projects/ENIGMA/data/annotations/refseq/transcriptome/icalvum.NC_014830.1.cs.fasta'
mo2refseqOrthologyFile='/Volumes/omics4tb/alomana/projects/ENIGMA/data/annotations/refseq/mo2refseq.txt'
refSeqGFF3='/Volumes/omics4tb/alomana/projects/ENIGMA/data/annotations/refseq/genome/icalvum.NC_014830.1.2018.05.22.gff3'

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

# 1.7. read KIP7 Microbes-online to Refseq BLAST for annotation mapping
intcaStone=mo2refseqReader()

# 1.8. read RefSeq annotation for old Intca IDs
newoldINTCAStone=refseqReader()

# 2. build a table of annotation
print('writing functional annotation file...')

g=open(functionalAnnotationFile,'w')
g.write('#C5.transcript.name\tC5.transcript.ID\tC5.transcript.description\tKIP7.ortholog.MO\tKIP7.ortholog.RefSeq\tKIP7.ortholog.RefSeq.old.ID\tGO\tKEGG\tPFAM\tCOG\tTIGR\n')

for transcriptName in transcriptNames:

    # write name
    line2Write='{}\t{}\t{}'.format(transcriptName,transcriptDescriptions[transcriptName][0],transcriptDescriptions[transcriptName][1])
    g.write(line2Write)

    # write GO
    if transcriptName in GOstone:
        g.write('\t{}'.format(orthologyStone[transcriptName]))
        # intca
        try:
            VIMSSID='VIMSS'+orthologyStone[transcriptName]
            intcaID=intcaStone[VIMSSID]
            old=newoldINTCAStone[intcaID]
        except:
            intcaID=' '
            old=' '
        g.write('\t{}\t{}'.format(intcaID,old))
        # GO
        g.write('\t{}'.format(GOstone[transcriptName]))
    else:
        g.write('\t ')
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

