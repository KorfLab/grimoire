import sys
import json

import toolbox.dna
from toolbox.gff import Gff
from toolbox.fasta import FastaFile
from hmm import HMM
import viterbi

gf = Gff('internal.gff')
ff = FastaFile('data/TAIR10_1.fasta')

seqs = []

for id in ff.ids:
	seq = ff.get(id).seq
	for feature in gf.get(chrom='Chr1', type='model'):
		if feature.strand == '-': continue
		seqs.append(seq[feature.beg : feature.end])

hmm = HMM.read('internal.hmm')

for seq in seqs:
	print(seq)
	path, score = viterbi.decode(model=hmm, seq=seq)
	print(path)
	print(score)
# 	paths, scores = viterbi.stochastic(model=hmm, seq=seq, n=100)
# 	for i in range(len(paths)):
# 		print(paths[i])
# 		print(scores[i])