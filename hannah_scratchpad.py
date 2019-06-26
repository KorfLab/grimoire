import sys
import json

import toolbox
from toolbox import GFF_file as GFF
from toolbox import FASTA_file as FASTA
from hmm import HMM

gf = GFF('internal.gff')
ff = FASTA('data/TAIR10_1.fasta')

seqs = []

for id in ff.ids:
	seq = ff.get(id).seq
	for feature in gf.get(chrom='Chr1', type='model'):
		if feature.strand == '-': continue
		seqs.append(seq[feature.beg : feature.end])

hmm = HMM.read('internal.hmm')

for seq in seqs:
	print(seq)
	path, score = hmm.decode(model=hmm, seq=seq)
	print(path)
	print(score)
# 	paths, scores = viterbi.stochastic(model=hmm, seq=seq, n=100)
# 	for i in range(len(paths)):
# 		print(paths[i])
# 		print(scores[i])