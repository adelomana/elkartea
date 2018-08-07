#!/bin/bash

#time python blast_rbh.py /Volumes/omics4tb/alomana/projects/ENIGMA/data/annotations/microbesOnline/transcriptome.txt /Volumes/omics4tb/alomana/projects/ENIGMA/data/annotations/C5/transcriptome/2767802313.multi.fasta -o /Volumes/omics4tb/alomana/projects/ENIGMA/data/annotations/microbesOnline/reciprocal.blast.microbesOnline.C5.output.txt -t blastn -a nucl --threads 8

time python blast_rbh.py /Volumes/omics4tb/alomana/projects/ENIGMA/data/annotations/microbesOnline/transcriptome.txt /Volumes/omics4tb/alomana/projects/ENIGMA/data/annotations/refseq/transcriptome/icalvum.NC_014830.1.cs.fasta -o /Volumes/omics4tb/alomana/projects/ENIGMA/data/annotations/refseq/mo2refseq.txt -t blastn -a nucl --threads 8
