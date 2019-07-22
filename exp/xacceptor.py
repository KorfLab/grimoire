#!/usr/bin/env python3

import argparse
import copy
import os
import re
import subprocess
import statistics
import sys

import grimoire.hmm as hmm
import grimoire.decode as decode
import grimoire.toolbox as toolbox
import grimoire.sequence as sequence
from grimoire.genome import Feature

parser = argparse.ArgumentParser(description='Experiment runner: acceptor')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<path>', help='path to fasta file (%(type)s)')
parser.add_argument('--gff3', required=True, type=str,
	metavar='<path>', help='path to GFF file (%(type)s)')
parser.add_argument('--work', required=True, type=str,
	metavar='<path>', help='path to the working directory (%(type)s)')
parser.add_argument('--xvalid', required=False, type=int, default=3,
	metavar='<int>', help='cross-validation depth [%(default)d]')
parser.add_argument('--overwrite', action='store_true',
	help='overwrite contents of work directory')
arg = parser.parse_args()

try:
	os.mkdir(arg.work)
except FileExistsError:
	if arg.overwrite:
		sys.stderr.write('Overwriting working directory: ' + arg.work + '\n') 
	else:
		sys.stderr.write('Working directory (' + arg.work
			+ ') already exists. Will not overwrite without --overwite flag.\n')
		sys.exit(0)


## Split data set into subsets ##

p = subprocess.run([
	'python3', 'bin/chromify',
	'--fasta',  'data/C.elegans.1percent.fasta.gz',
	'--gff3',  'data/C.elegans.1percent.gff3.gz',
	'--out', arg.work + '/set',
	'--split', str(arg.xvalid)])
p.check_returncode()

## Run cross-validation experiments ##
accuracy = []
for i in range(arg.xvalid):
	test_dna = arg.work + '/set-' + str(i) + '.fasta'
	test_gff = arg.work + '/set-' + str(i) + '.gff3'
	train_dna = arg.work + '/train-' + str(i) + '.fasta'
	train_gff = arg.work + '/train-' + str(i) + '.gff3'
	dna = ''
	gff = ''
	for j in range(arg.xvalid):
		if i != j:
			dna += arg.work + '/set-' + str(j) + '.fasta '
			gff += arg.work + '/set-' + str(j) + '.gff3 '
	os.system('cat ' + dna + '> ' + train_dna)
	os.system('cat ' + gff + '> ' + train_gff)
		
	# create HMM
	hmm_file = arg.work + '/' + str(i) + '.hmm'
	subprocess.run([
		'bin/forge',
		'--fasta',  train_dna,
		'--gff3',   train_gff,
		'--model',  'acc',
		'--hmm',    hmm_file])
	model = hmm.HMM.read(hmm_file)
	model.convert2log()
		
	# create test sets
	source = arg.work + '/' + 'test' + '-' + str(i)
	subprocess.run([
		'bin/forge',
		'--fasta',  train_dna,
		'--gff3',   train_gff,
		'--model',  'acc',
		'--source', source])
		
	# decode and evalute HMM performance
	perf = decode.Performance(model)
	fasta = toolbox.FASTA_stream(filename=source + '.fasta')
	gff = toolbox.GFF_file(filename=source + '.gff')
	for entry in fasta:
		dna = sequence.DNA(seq=entry.seq, name=entry.id)
		gffs = gff.get(chrom=dna.name)
		truth = []
		for g in gffs:
			truth.append(Feature(dna, g.beg, g.end, g.strand, g.type))
		v = decode.Viterbi(model=model, dna=dna)
		p = v.generate_path()
		predicted = p.features()
		perf.compare(source=truth, prediction=predicted)

	ac = perf.nt_same / (perf.nt_same + perf.nt_diff)
	print(i, ac)
	accuracy.append(ac)

print('mean', statistics.mean(accuracy))