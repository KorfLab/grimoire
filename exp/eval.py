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
import grimoire.toolbox as toolbox
import grimoire.sequence as sequence
import grimoire.genome as genome
import grimoire.decode as decode
import grimoire.hmm as hmm

#Command Line
parser = argparse.ArgumentParser(description='Evaluate 3 Different Models')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<path>', help='path to fasta file (%(type)s)')
parser.add_argument('--gff3', required=True, type=str,
	metavar='<path>', help='path to GFF file (%(type)s)')
parser.add_argument('--mRNA1', required=True, type=str,
	metavar='<path>', help = 'path to the mRNA 1 HMM model. if HMM path not found it'+
	' will be generated for you (%(type)s)')
parser.add_argument('--mRNA2', required=True, type=str,
	metavar='<path>', help = 'path to the mRNA 2 HMM model. if HMM path not found it'+
	' will be generated for you (%(type)s)')
parser.add_argument('--overwrite', action ='store_true', required=False,
	help='overwrite previous HMM')
parser.add_argument('--work', required=True, type=str,
	metavar='<path>', help='path to the working directory (%(type)s)')
arg = parser.parse_args()

if (os.path.exists(arg.mRNA1)==False or os.path.exists(arg.mRNA1)==False) or arg.overwrite==True:
	print('Overwriting/Creating')
	os.system(
		'python3 bin/forge '+
		'--fasta '+arg.fasta+' '+
		'--gff3 '+arg.gff3+' '+
		'--model '+'mRNA1'+' '+
		'--hmm '+arg.mRNA1+' '+
		'--source '+arg.work+'C.elegans.mRNA1.source'+' '+
		'--u5_ctx '+'5 '+
		'--u3_ctx '+'5 '+
		'--koz_ctx '+'1 '+
		'--atg_ctx '+'2 '+
		'--stop_ctx '+'2 '+
		'--cds_ctx '+'2'
		)
	os.system(
		'python3 bin/forge '+
		'--fasta '+arg.fasta+' '+
		'--gff3 '+arg.gff3+' '+
		'--model '+'mRNA2'+' '+
		'--hmm '+arg.mRNA2+' '+
		'--source '+arg.work+'C.elegans.mRNA2.source'+' '+
		'--u5_ctx '+'5 '+
		'--u3_ctx '+'5 '+
		'--koz_ctx '+'1 '+
		'--atg_ctx '+'2 '+
		'--stop_ctx '+'2 '+
		'--cds_ctx '+'2'
		)
	print('Forge Completed')

fasta=toolbox.FASTA_stream(filename = arg.work+'C.elegans.mRNA1.source.fasta')
gff =toolbox.GFF_file(filename = arg.work+'C.elegans.mRNA1.source.gff')

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
mRNA1_genes = {}
mRNA1_hmm = hmm.HMM.read(arg.work+'C.elegans.mRNA1.hmm')
mRNA1_hmm.convert2log()
for dna in dnas:
	mRNA1 = decode.Viterbi(model=mRNA1_hmm,dna=dna)
	path = mRNA1.generate_path()
	features = path.features()
	for f in features:
		if f.type == 'ATG':
			mRNA1_genes.update({f.dna.name:f.beg})

#mRNA2 model
mRNA2_genes = {}
mRNA2_hmm = hmm.HMM.read(arg.work+'C.elegans.mRNA2.hmm')
mRNA2_hmm.convert2log()
for dna in dnas:
	mRNA2 = decode.Viterbi(model=mRNA2_hmm,dna=dna)
	path = mRNA2.generate_path()
	features = path.features()
	for f in features:
		if f.type == 'ATG':
			mRNA2_genes.update({f.dna.name:f.beg})

#Report
print('Id'+'\t \t \t \t'+'LORF'+'\t'+'mRNA1'+'\t'+'mRNA2')
for key in GFF3_genes.keys():
	print(key + '\t' + str(GFF3_genes[key])+'\t'+str(mRNA1_genes[key])+'\t'+str(mRNA2_genes[key]))
