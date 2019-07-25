#WORK IN PROGRESS
#PROCEED WITH CAUTION BC STUFF MIGHT BE BROKEN WHOO

#imports
from random import *
import matplotlib.pyplot as plt
import numpy as np
import re
import statistics
import os
import json
import operator
import argparse
import grimoire.genome as genome
import grimoire.decode as decode
import grimoire.hmm as hmm
import grimoire.toolbox as toolbox
import grimoire.sequence as sequence

#Command Line
parser = argparse.ArgumentParser(description='Evaluate 3 Different Models')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<path>', help='path to fasta file (%(type)s)')
parser.add_argument('--gff3', required=True, type=str,
	metavar='<path>', help='path to GFF file (%(type)s)')
parser.add_argument('--work', required=True, type=str,
	metavar='<path>', help='path to the working directory (%(type)s)')
arg = parser.parse_args()

os.system(
    'python3 bin/forge '+
    '--fasta '+arg.fasta+' '+
    '--gff3 '+arg.gff3+' '+
    '--model '+'mRNA1'+' '+
    '--hmm '+arg.work+'C.elegans.1percent.mRNA1.hmm'+' '+
    '--source '+arg.work+'C.elegans.1percent.mRNA1.source'
    )

os.system(
    'python3 bin/forge '+
    '--fasta '+arg.fasta+' '+
    '--gff3 '+arg.gff3+' '+
    '--model '+'mRNA2'+' '+
    '--hmm '+arg.work+'C.elegans.1percent.mRNA2.hmm'+' '+
    '--source '+arg.work+'C.elegans.1percent.mRNA2.source'
    )

'''
os.system(
    'python3 bin/chromify '+
    '--fasta '+'data/C.elegans.1percent.fasta'+' '+
    '--gff3 '+'data/C.elegans.1percent.gff3'+' '+
    '--out '+'data/C.elegans.source'
    )
'''

fasta=toolbox.FASTA_stream(filename = arg.work+'C.elegans.1percent.mRNA1.source.fasta')
gff =toolbox.GFF_file(filename = arg.work+'C.elegans.1percent.mRNA1.source.gff')

#Extract from GFF3
GFF3_genes = {}
dnas = []
for entry in fasta:
    dna = sequence.DNA(seq=entry.seq, name=entry.id)
    dnas.append(dna)
    gffs = gff.get(chrom=dna.name)
    for g in gffs:
        if g.type == 'ATG':
            GFF3_genes.update({entry.id:g.beg})

#mRNA1 model
mRNA1_genes = []
mRNA1_hmm = hmm.HMM.read(arg.work+'C.elegans.1percent.mRNA1.hmm')
mRNA1_hmm.convert2log()
for dna in dnas:
    mRNA1 = decode.Viterbi(model=mRNA1_hmm,dna=dna)
    path = mRNA1.generate_path()
    features = path.features()
    for f in features:
        if f.type == 'ATG':
            mRNA1_genes.append(f.beg)

#mRNA2 model
mRNA2_genes = []
mRNA2_hmm = hmm.HMM.read(arg.work+'C.elegans.1percent.mRNA2.hmm')
mRNA2_hmm.convert2log()
for dna in dnas:
    mRNA2 = decode.Viterbi(model=mRNA2_hmm,dna=dna)
    path = mRNA2.generate_path()
    features = path.features()
    for f in features:
        if f.type == 'ATG':
            mRNA2_genes.append(f.beg)

#Report
print('Id'+'\t \t \t \t'+'LORF'+'\t'+'mRNA1'+'\t'+'mRNA2')
for i, key in enumerate(GFF3_genes.keys()):
    print(key + '\t' + str(GFF3_genes[key])+'\t'+str(mRNA1_genes[i])+'\t'+str(mRNA2_genes[i]))
