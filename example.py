#!/usr/bin/python3

"""test harness for grimoire"""

import sys
import json
import toolbox.dna

### FASTA ###

from toolbox.fasta import FastaFile
# ff = Fasta('data/TAIR10_1.fasta')
# for id in ff.ids:
# 	entry = ff.get(id)
# 	print(entry.id, entry.desc, entry.seq[0:50])

### GFF ###

# from toolbox.gff import Gff
# gf = Gff('data/TAIR10_1.gff3')
# print(gf.types)
# print(gf.chroms)
# for f in gf.get(chrom='Chr1', beg=7500, end=8500):
# 	print(f.beg, f.end, f.type, f.strand)

## DNA ##

# seq = 'ACCCCGAGGAGAGGACCCATAGGC'
# rev = toolbox.dna.revcomp(seq)
# print(seq, rev)

## HMM ##
from hmm import HMM
import viterbi
hmm = HMM.read('toy.hmm')
ff = FastaFile('toy.fasta')
for id in ff.ids:
	seq = ff.get(id).seq
	print(seq)
# 	path, score = viterbi.decode(model=hmm, seq=seq)
# 	print(path)
# 	print(score)
	paths, scores = viterbi.stochastic(model=hmm, seq=seq, n=10)
	for i in range(len(paths)):
		print(paths[i])
		print(scores[i])
	break